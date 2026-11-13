"""Chemistry primitives. One place where a molecule becomes an identity.

Two decisions live here and nowhere else, because having them in two places is how
a leakage bug gets written:

IDENTITY
    A metabolite is identified by the first block of its InChIKey (14 chars, the
    skeleton block). This collapses stereochemistry, protonation and charge state.
    That is deliberate: Rhea writes ATP as a 4-minus anion and KEGG writes it
    neutral, and they are the same molecule for every question this project asks.
    Where stereochemistry matters (it does for the verifier's features) the full
    key is available via ikey_full.

CURRENCY
    PLAN.md 6.3, frozen: a participant is currency if it is in the curated name
    list OR its reaction-degree is above the 99.5th percentile. Degree is counted
    over LR-direction reactions only - Rhea stores each reaction twice (LR and RL),
    and counting both inflates every degree by exactly 2x while changing no ranking,
    but it also makes the printed numbers lies, so we do not do it.

    Currency metabolites are removed as graph NODES and retained as verifier
    FEATURES. This module only decides membership; removal happens at graph build.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, RDLogger

from . import config

# RDKit is loud about the wildcard-containing SMILES Rhea uses for macromolecular
# participants. We handle those explicitly below; silence the per-molecule noise.
RDLogger.DisableLog("rdApp.*")

WILDCARD = "*"


# ------------------------------------------------------------------ identity
@lru_cache(maxsize=1 << 20)
def parse(smi: str) -> Chem.Mol | None:
    """SMILES -> sanitised Mol, or None. Cached: the same metabolite recurs
    thousands of times across the reaction set."""
    if not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    return m


@lru_cache(maxsize=1 << 20)
def ikey_full(smi: str) -> str | None:
    m = parse(smi)
    if m is None:
        return None
    try:
        k = Chem.MolToInchiKey(m)
    except Exception:  # noqa: BLE001 - InChI fails on wildcards and odd valences
        return None
    return k or None


def ikey14(smi: str) -> str | None:
    """Skeleton block of the InChIKey. The project's metabolite identity."""
    k = ikey_full(smi)
    return k[:14] if k else None


def is_generic(smi: str) -> bool:
    """True for Rhea's macromolecular/polymer participants, written with a '*'
    attachment point (e.g. a lysine residue inside a protein). These have no
    well-defined InChIKey and no RetroRules template can be scored against them,
    so every experiment reports them as a separate bucket rather than a failure."""
    return WILDCARD in smi


# ------------------------------------------------------------------ reactions
def split_reaction(rxn_smiles: str) -> tuple[list[str], list[str]]:
    """'A.B>>C.D' -> (['A','B'], ['C','D']). Rhea uses no agent field."""
    if ">>" not in rxn_smiles:
        return [], []
    left, right = rxn_smiles.split(">>", 1)
    return (
        [s for s in left.split(".") if s],
        [s for s in right.split(".") if s],
    )


# ------------------------------------------------------------------ currency
CURRENCY_CACHE = config.DATA_PROCESSED / "currency_set.json"


def reaction_participants(reaction_smiles: dict[str, str]) -> dict[str, set[str]]:
    """reaction id -> set of participant InChIKey-14. Generics and unparsable
    participants are dropped; counts of what was dropped come back via the
    companion function below so no loss is silent."""
    out: dict[str, set[str]] = {}
    for rid, rxn in reaction_smiles.items():
        subs, prods = split_reaction(rxn)
        seen: set[str] = set()
        for smi in subs + prods:
            if is_generic(smi):
                continue
            k = ikey14(smi)
            if k:
                seen.add(k)
        out[rid] = seen
    return out


def participant_losses(reaction_smiles: dict[str, str]) -> tuple[int, int]:
    """(generic participants skipped, participants with no derivable InChIKey)."""
    generic = unparsed = 0
    for rxn in reaction_smiles.values():
        subs, prods = split_reaction(rxn)
        for smi in subs + prods:
            if is_generic(smi):
                generic += 1
            elif ikey14(smi) is None:
                unparsed += 1
    return generic, unparsed


def degree_table(reaction_smiles: dict[str, str]) -> Counter:
    """How many distinct reactions each metabolite appears in.

    Superseded as a currency criterion (docs/DECISIONS.md D1) but still reported:
    it is the quantity the plan's original rule used, and the revision is only
    auditable against what it replaced.
    """
    degree: Counter = Counter()
    for seen in reaction_participants(reaction_smiles).values():
        for k in seen:
            degree[k] += 1
    return degree


def breadth_table(
    candidates: dict[str, dict], participants: dict[str, set[str]]
) -> dict[str, set[str]]:
    """metabolite -> set of ordinary KEGG maps it appears in.

    candidates: data.candidate_set(), i.e. kegg_rid -> {rhea_master, pathways, ...}
    participants: keyed by Rhea MASTER id.
    """
    out: dict[str, set[str]] = defaultdict(set)
    for rec in candidates.values():
        for k in participants.get(rec["rhea_master"], ()):  # type: ignore[arg-type]
            out[k].update(rec["pathways"])
    return dict(out)


def verify_curated_formulae(example_smiles: dict[str, str]) -> list[str]:
    """Check each curated cofactor's InChIKey-14 still names the stated molecule.

    A curated list keyed by InChIKey with the formula written only in a comment is
    a list that silently rots. Returns the problems found; the caller decides
    whether to raise. An entry absent from the data is reported, not an error -
    GTP need not occur in every subset of reactions.
    """
    from rdkit.Chem import rdMolDescriptors

    problems: list[str] = []
    for key, (formula, _why) in config.CURATED_COFACTORS.items():
        smi = example_smiles.get(key)
        if smi is None:
            problems.append(f"{key} ({formula}): not present in this reaction set")
            continue
        m = parse(smi)
        got = rdMolDescriptors.CalcMolFormula(m) if m else "?"
        if got != formula:
            problems.append(
                f"{key}: curated as {formula} but the data says {got} - "
                "the key no longer names the molecule it claims to"
            )
    return problems


def currency_from_breadth(
    breadth: dict[str, set[str]],
    *,
    threshold: int = config.CURRENCY_BREADTH_MIN,
) -> set[str]:
    """The revised rule (PLAN.md 6.3, docs/DECISIONS.md D1):

        currency = {breadth >= threshold} OR curated cofactor list

    Threshold chosen for zero false positives on the labelled references, not for
    best accuracy: a false positive deletes a real metabolite from the graph, a
    false negative only leaves a cofactor hub in it.
    """
    auto = {k for k, maps in breadth.items() if len(maps) >= threshold}
    return auto | set(config.CURATED_COFACTORS)


def load_currency() -> set[str]:
    """The frozen currency set as InChIKey-14. Fails loudly if not yet computed:
    a silently-empty currency set would quietly change every downstream number."""
    if not CURRENCY_CACHE.exists():
        raise FileNotFoundError(
            f"{CURRENCY_CACHE} missing - run experiments/03_currency_audit.py first"
        )
    return set(json.loads(CURRENCY_CACHE.read_text())["currency_ikey14"])


def non_currency(smiles_list: list[str], currency: set[str]) -> list[str]:
    """Drop currency participants and generics; keep order, keep duplicates out."""
    out: list[str] = []
    seen: set[str] = set()
    for smi in smiles_list:
        if is_generic(smi):
            continue
        k = ikey14(smi)
        if k is None or k in currency or k in seen:
            continue
        seen.add(k)
        out.append(smi)
    return out
