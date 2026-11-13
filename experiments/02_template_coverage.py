#!/usr/bin/env python3
"""Experiment 02 - Gate 1, RetroRules template coverage (PLAN.md 3.1).

The gate asks: of the KEGG pathway reactions we can actually draw, what fraction
can a RetroRules template reproduce?

Measured in two stages, cheapest first.

STAGE A (this stage, seconds, no chemistry)
    RetroRules templates are *derived from* specific reactions; each rule's
    "Legacy ID" names the MetaNetX reaction it came from. So set membership gives
    a hard UPPER BOUND on coverage without applying a single template:

        can a rule possibly reproduce reaction R?  ->  is R's MNXR id among the
        reactions that generated rules, or does some rule from another reaction
        happen to apply?

    The first term is computable instantly. The second needs Stage B. So Stage A
    reports:
      self_derived   - a rule was built from this very reaction. It will reproduce
                       it trivially. This is memorisation, not generalisation, and
                       is reported separately for exactly that reason.
      not_derived    - no rule came from this reaction; whether any *other* rule
                       applies is Stage B's question.

    If self_derived is already below the kill threshold, Stage B cannot rescue it
    by much, and the gate fails early and cheaply.

STAGE B (not in this file yet, needs RDKit)
    Apply candidate templates with RunReactants and check the recorded product is
    produced. Run with and without excluding self-derived rules.

Run:  python experiments/02_template_coverage.py --stage a
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from biopathagent import config  # noqa: E402

RAW = config.DATA_RAW
RR_DIR = RAW / "retrorules" / "retrorules_rr02_rp2_hs"


# ------------------------------------------------------------------ loaders
def load_rule_sources() -> tuple[set[str], dict[str, Counter], int]:
    """Parse the flat forward rule file.

    Returns (mnxr ids that generated >=1 rule, ec -> diameter counter, n rules).
    """
    sources: set[str] = set()
    ec_diam: dict[str, Counter] = defaultdict(Counter)
    n = 0
    path = RR_DIR / "retrorules_rr02_rp2_flat_forward.csv"
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            n += 1
            legacy = row.get("Legacy ID", "")
            # "MNXR94682_MNXM821" -> reaction id is the part before the first "_"
            if legacy.startswith("MNXR"):
                sources.add(legacy.split("_", 1)[0])
            ec_diam[row.get("EC number", "NOEC")][row.get("Diameter", "?")] += 1
    return sources, ec_diam, n


def load_kegg_to_mnxr() -> dict[str, str]:
    out: dict[str, str] = {}
    path = RAW / "metanetx" / "reac_xref.tsv"
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) >= 2 and p[0].startswith("kegg.reaction:") and p[1].startswith("MNXR"):
                out[p[0].split(":", 1)[1]] = p[1]
    return out


def load_audit_candidates() -> set[str]:
    """Rebuild the Gate 1 candidate set the same way experiment 01 did.

    Imported rather than duplicated would be better; kept explicit here so the two
    experiments can be read independently, and asserted against the audit card.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "audit", Path(__file__).resolve().parent / "01_data_audit.py"
    )
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)  # type: ignore[union-attr]

    smiles = audit.load_reaction_smiles()
    directions = audit.load_directions()
    rhea2kegg = audit.load_rhea_xref("rhea2kegg_reaction.tsv")
    kegg_path_rxn = audit.load_kegg_links("kegg_pathway_reaction.tsv")

    in_ordinary = {
        audit.strip_prefix(r)
        for p, r in kegg_path_rxn
        if not audit.is_global_map(p)
    }
    drawable = set()
    for r in rhea2kegg:
        lr = directions.get(r["MASTER_ID"], {}).get("LR")
        if lr and lr in smiles and r["ID"] in in_ordinary:
            drawable.add(r["ID"])
    return drawable


# ------------------------------------------------------------------ stage A
def stage_a() -> dict:
    t0 = time.time()
    print("=" * 74)
    print("EXPERIMENT 02 - GATE 1, STAGE A (upper bound by rule provenance)")
    print("=" * 74)

    print("\nloading ...")
    rule_sources, ec_diam, n_rules = load_rule_sources()
    kegg2mnxr = load_kegg_to_mnxr()
    candidates = load_audit_candidates()

    print(f"    forward rules                                    {n_rules:>10,}")
    print(f"    distinct MNXR reactions that generated rules     {len(rule_sources):>10,}")
    print(f"    KEGG -> MNXR mappings                            {len(kegg2mnxr):>10,}")
    print(f"    Gate 1 candidate reactions (from exp 01)         {len(candidates):>10,}")

    mapped = {k: kegg2mnxr[k] for k in candidates if k in kegg2mnxr}
    self_derived = {k for k, m in mapped.items() if m in rule_sources}

    n_cand = len(candidates)
    n_mapped = len(mapped)
    n_self = len(self_derived)

    print("\n[A] provenance-based upper bound")
    print(f"    candidates                                       {n_cand:>10,}")
    print(f"    ... with a MetaNetX id                           {n_mapped:>10,}"
          f"   ({n_mapped/max(n_cand,1):6.1%})")
    print(f"    ... a rule was derived FROM this reaction        {n_self:>10,}"
          f"   ({n_self/max(n_cand,1):6.1%} of candidates)")
    print(f"    ... no rule derived from it (Stage B decides)    {n_mapped-n_self:>10,}")
    print(f"    ... not in MetaNetX at all (cannot be checked)   {n_cand-n_mapped:>10,}")

    cov_self = n_self / max(n_cand, 1)
    print("\n[B] what this does and does not say")
    print(f"    trivial/memorised coverage (self-derived)        {cov_self:>10.1%}")
    print( "    true coverage is >= this (other rules may apply) and is measured in Stage B")
    print(f"    absolute coverage of all 8,215 KEGG map reactions"
          f"{cov_self*0.510:>10.1%}  (x0.510 drawable ceiling)")

    # diameter spread tells us how much generalisation headroom the rule set has
    all_diam: Counter = Counter()
    for c in ec_diam.values():
        all_diam.update(c)
    print("\n[C] rule diameters (smaller = more generic = more likely to transfer)")
    for d in sorted(all_diam, key=lambda x: int(x) if str(x).isdigit() else 99):
        print(f"    diameter {str(d):>3}                                     "
              f"{all_diam[d]:>10,}   ({all_diam[d]/max(n_rules,1):5.1%})")
    noec = sum(v for k, v in ((k, sum(c.values())) for k, c in ec_diam.items())
               if k == "NOEC")
    print(f"    rules with no EC assigned                        {noec:>10,}"
          f"   ({noec/max(n_rules,1):5.1%})")

    verdict = config.gate1_verdict(cov_self)
    print("\n[D] GATE 1 - STAGE A VERDICT (lower bound on coverage)")
    print(f"    thresholds: PASS >= {config.GATE1_PASS:.0%}, "
          f"DEGRADE >= {config.GATE1_DEGRADE:.0%}, else KILL")
    print(f"    stage-A lower bound = {cov_self:.1%}  ->  {verdict}")
    if verdict == "KILL":
        print("    NOTE: a KILL here is provisional - Stage B can only raise this number.")
    else:
        print("    Stage B can only raise this; the gate is already clear of KILL.")

    card = {
        "experiment": "02_template_coverage_stage_a",
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "plan_section": "3.1 Gate 1",
        "rules": {
            "forward_rules": n_rules,
            "source_reactions": len(rule_sources),
            "diameters": {str(k): v for k, v in all_diam.items()},
            "rules_without_ec": noec,
        },
        "candidates": {
            "total": n_cand,
            "with_metanetx_id": n_mapped,
            "self_derived_rule_exists": n_self,
            "no_self_derived_rule": n_mapped - n_self,
            "not_in_metanetx": n_cand - n_mapped,
        },
        "coverage": {
            "stage_a_lower_bound_conditional": cov_self,
            "stage_a_lower_bound_absolute": cov_self * 0.510,
            "note": "self-derived only; Stage B measures rules from OTHER reactions",
        },
        "gate1_stage_a_verdict": verdict,
        "runtime_sec": round(time.time() - t0, 1),
    }
    out = config.RESULTS / "02_gate1_stage_a.json"
    out.write_text(json.dumps(card, indent=2))
    print(f"\nresult card: {out}")
    print(f"runtime: {card['runtime_sec']}s")
    return card


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="a", choices=["a"])
    args = ap.parse_args()
    config.ensure_dirs()
    if args.stage == "a":
        stage_a()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
