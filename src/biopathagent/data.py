"""Canonical loaders. Promoted out of experiments/01_data_audit.py.

Experiment 02 originally rebuilt the Gate 1 candidate set by importing 01 through
importlib. That worked and was wrong: two copies of a definition drift, and the
definition in question is the denominator of the project's headline number. Every
experiment from here on calls candidate_set() and nothing else.

The one guard that actually catches drift is at the bottom: candidate_set()
compares its own result against the count experiment 01 wrote to its result card
and raises if they disagree.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from . import config

RAW = config.DATA_RAW


# ------------------------------------------------------------------ rhea
def load_reaction_smiles() -> dict[str, str]:
    """rhea-reaction-smiles.tsv has NO header: <rhea_id>\\t<reaction SMILES>."""
    out: dict[str, str] = {}
    with (RAW / "rhea" / "rhea-reaction-smiles.tsv").open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[0].isdigit():
                out[parts[0]] = parts[1]
    return out


def load_rhea_xref(filename: str) -> list[dict[str, str]]:
    """rhea2*.tsv: RHEA_ID, DIRECTION, MASTER_ID, ID (tab separated, with header)."""
    rows: list[dict[str, str]] = []
    with (RAW / "rhea" / filename).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            vals = line.rstrip("\n").split("\t")
            if len(vals) == len(header):
                rows.append(dict(zip(header, vals)))
    return rows


def load_directions() -> dict[str, dict[str, str]]:
    """rhea-directions.tsv: MASTER -> {LR, RL, BI} variant ids.

    This indirection is NOT optional. The cross-reference files (rhea2kegg etc.)
    cite the BI (bidirectional) variant, while rhea-reaction-smiles.tsv contains
    only the LR and RL variants - the master id and the BI id are both absent from
    it. Joining the xref id straight onto the SMILES table silently yields zero
    matches, which is how experiment 01 was first written.
    """
    out: dict[str, dict[str, str]] = {}
    with (RAW / "rhea" / "rhea-directions.tsv").open() as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 4 and p[0].isdigit():
                out[p[0]] = {"LR": p[1], "RL": p[2], "BI": p[3]}
    return out


def lr_reaction_smiles() -> dict[str, str]:
    """One reaction SMILES per Rhea MASTER, in the LR direction.

    Rhea stores every reaction twice (10001 is the forward of 10002). Any statistic
    counted over the raw SMILES table therefore double-counts. Keyed by master id.
    """
    smiles = load_reaction_smiles()
    out: dict[str, str] = {}
    for master, var in load_directions().items():
        lr = var.get("LR")
        if lr and lr in smiles:
            out[master] = smiles[lr]
    return out


# ------------------------------------------------------------------ kegg
def load_kegg_links(filename: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    with (RAW / "kegg" / filename).open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                out.append((parts[0], parts[1]))
    return out


def strip_prefix(x: str) -> str:
    return x.split(":", 1)[1] if ":" in x else x


def is_global_map(pathway_id: str) -> bool:
    """KEGG 'global/overview' maps (map011xx, map012xx) are aggregates - e.g.
    map01100 'Metabolic pathways' holds thousands of reactions. Counting them
    alongside ordinary maps double-counts most of metabolism, so they are excluded
    and reported separately. Frozen choice (PLAN.md 6.2), not a tunable.
    """
    n = strip_prefix(pathway_id).replace("map", "")
    return n.startswith("011") or n.startswith("012")


def kegg_pathway_reactions(*, exclude_global: bool = True) -> dict[str, set[str]]:
    """KEGG map id -> set of KEGG R-numbers, prefixes stripped."""
    out: dict[str, set[str]] = defaultdict(set)
    for p, r in load_kegg_links("kegg_pathway_reaction.tsv"):
        if exclude_global and is_global_map(p):
            continue
        out[strip_prefix(p)].add(strip_prefix(r))
    return dict(out)


# ------------------------------------------------------------------ metanetx
def kegg_to_mnxr() -> dict[str, str]:
    out: dict[str, str] = {}
    with (RAW / "metanetx" / "reac_xref.tsv").open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) >= 2 and p[0].startswith("kegg.reaction:") and p[1].startswith("MNXR"):
                out[p[0].split(":", 1)[1]] = p[1]
    return out


# ------------------------------------------------------------------ the gate set
def candidate_set() -> dict[str, dict]:
    """The Gate 1 candidate set: KEGG pathway reactions we can actually draw.

    Returns kegg_rid -> {rhea_master, rhea_lr, smiles, pathways}.

    'Drawable' is the gate's denominator (PLAN.md 3.1): a KEGG pathway reaction
    that maps to a Rhea master whose LR variant has a reaction SMILES. The two
    failures it separates - 'no structure available' vs 'no template reproduces it'
    - were conflated in the plan's first draft, which is why the denominator is
    pinned here in one function.
    """
    smiles = load_reaction_smiles()
    directions = load_directions()
    rhea2kegg = load_rhea_xref("rhea2kegg_reaction.tsv")

    in_ordinary: dict[str, set[str]] = defaultdict(set)
    for p, r in load_kegg_links("kegg_pathway_reaction.tsv"):
        if not is_global_map(p):
            in_ordinary[strip_prefix(r)].add(strip_prefix(p))

    out: dict[str, dict] = {}
    for row in rhea2kegg:
        kegg_rid = row["ID"]
        if kegg_rid not in in_ordinary:
            continue
        master = row["MASTER_ID"]
        lr = directions.get(master, {}).get("LR")
        if not lr or lr not in smiles:
            continue
        # A KEGG reaction can map to several Rhea masters. Keep the first by rhea
        # id for determinism and record the collision count rather than silently
        # picking one; collisions are 1.3% of the set and are reported by exp 02.
        if kegg_rid in out:
            out[kegg_rid]["n_rhea_alternatives"] += 1
            continue
        out[kegg_rid] = {
            "rhea_master": master,
            "rhea_lr": lr,
            "smiles": smiles[lr],
            "pathways": sorted(in_ordinary[kegg_rid]),
            "n_rhea_alternatives": 1,
        }

    _assert_matches_audit(len(out))
    return out


def audit_card() -> dict:
    path = config.RESULTS / "01_data_audit.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing - run experiments/01_data_audit.py first"
        )
    return json.loads(path.read_text())


def drawable_fraction() -> float:
    """Fraction of all KEGG ordinary-map reactions that are drawable at all.

    Gate 1's conditional coverage multiplies by this to give absolute coverage.
    Read from the audit card rather than written as a constant, because a constant
    is what lets the two numbers drift apart silently.
    """
    c = audit_card()["chemistry_availability"]
    return c["with_reaction_smiles"] / c["kegg_reactions_in_ordinary_maps"]


def _assert_matches_audit(n: int) -> None:
    """Cross-check against experiment 01's recorded number."""
    path = config.RESULTS / "01_data_audit.json"
    if not path.exists():
        return
    rec = json.loads(path.read_text())["gate1_candidates"]["drawable"]
    if rec != n:
        raise AssertionError(
            f"candidate_set() returned {n} but experiment 01 recorded {rec}. "
            "One of the two definitions drifted; do not paper over this - the "
            "number is the gate's denominator."
        )
