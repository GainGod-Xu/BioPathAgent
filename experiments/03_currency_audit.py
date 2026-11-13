#!/usr/bin/env python3
"""Experiment 03 - freeze the currency-metabolite set (PLAN.md 6.3, Tier A).

Why this runs before Gate 1 Stage B
-----------------------------------
"Does a template reproduce this reaction?" is meaningless without saying which
product counts. For ATP + glucose -> ADP + glucose-6-P, a template applied to
glucose should be judged on whether it makes glucose-6-P. If ADP counts, then a
phosphate-transfer template "succeeds" on every kinase in the database by
producing the cofactor half, and coverage becomes a measurement of how common ATP
is. So the currency set is a dependency of the gate, not a later refinement.

The rule, as revised (docs/DECISIONS.md D1)
-------------------------------------------
    currency(m) = breadth(m) >= 59  OR  m in config.CURATED_COFACTORS

breadth(m) = number of distinct ordinary KEGG maps containing a reaction m takes
part in. 59 is the smallest threshold with no labelled pathway intermediate above
it. The threshold is set for zero false positives rather than best accuracy,
because the two errors are not symmetric: deleting pyruvate from the graph
corrupts every pathway result, while leaving FAD in it adds one hub node.

This file reports three things it would be dishonest to omit:
  [B]  the superseded degree-percentile rule's output, so the revision can be
       checked against what it replaced
  [E]  the threshold sensitivity over 40..100, so "59" can be seen for what it is
       - a value fitted to 27 hand-labelled references, not a cross-validated one
  [F]  the curated list's formula assertions, so a stale InChIKey fails loudly

Run:  python experiments/03_currency_audit.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rdkit.Chem import rdMolDescriptors  # noqa: E402

from biopathagent import chem, config, data  # noqa: E402

# Labelled references. Same set as experiment 04; imported from there would be
# tidier, but these labels are the calibration data for a frozen threshold and are
# written out in full on purpose.
REF_COFACTOR = {
    "GPRLSGONYQIRFK", "XLYOFNOQVPJJNP", "ZKHQWZAMYRWXGA", "XTWYTFMLZFPYCI",
    "UDMBCSSLTHHNCD", "BAWFJGJZGIEFAR", "BOPGDPNILDQYTO", "XJLXINKUBYWONI",
    "ACFIXJIJDZMPPO", "RGJOEKWQDUBAIZ", "NBIIXXVUZAFLBC", "XPPKVPWEQAFLFU",
    "MEFKEPWMEQBLKI", "ZJUKTBDSGOFHSH", "CURLTUGMZLYLDI", "MYMOFIZGZYHOMD",
    "IMGVNJNCCGXBHD", "ANKZYBDXHMZBDK",
}
REF_INTERMEDIATE = {
    "WQZGKKKJIJFFOK", "LCTONWCANYUPML", "KPGXRSRHYNQIFN", "KDYFGRWQOYBRFD",
    "QTBSBXVTEAMEQO", "QNAYBMKLOCPYGJ", "WHUUTDBJXJRKMK", "ZSLZBFCDCINBPY",
    "ZQPPMHVWECSIRJ",
}


def main() -> int:
    config.ensure_dirs()
    t0 = time.time()
    print("=" * 84)
    print("EXPERIMENT 03 - CURRENCY METABOLITES (revised rule, PLAN.md 6.3 / D1)")
    print("=" * 84)

    rxns = data.lr_reaction_smiles()
    parts = chem.reaction_participants(rxns)
    cands = data.candidate_set()
    breadth = chem.breadth_table(cands, parts)
    degree = chem.degree_table(rxns)
    n_generic, n_unparsed = chem.participant_losses(rxns)
    n_maps = len({p for r in cands.values() for p in r["pathways"]})

    example: dict[str, str] = {}
    for rxn in rxns.values():
        for smi in sum(chem.split_reaction(rxn), []):
            if not chem.is_generic(smi):
                k = chem.ikey14(smi)
                if k and k not in example:
                    example[k] = smi

    def formula(k: str) -> str:
        m = chem.parse(example.get(k, ""))
        return rdMolDescriptors.CalcMolFormula(m) if m else "?"

    print(f"\nreactions counted (one LR variant per Rhea master)   {len(rxns):>10,}")
    print(f"drawable KEGG pathway reactions (gate denominator)   {len(cands):>10,}")
    print(f"ordinary KEGG maps spanned                           {n_maps:>10,}")
    print(f"distinct metabolites (InChIKey-14)                   {len(degree):>10,}")
    print(f"generic (wildcard '*') participants skipped          {n_generic:>10,}")
    print(f"participants with no derivable InChIKey              {n_unparsed:>10,}")

    # -------------------------------------------------------------- the set
    currency = chem.currency_from_breadth(breadth)
    auto = {k for k, m in breadth.items() if len(m) >= config.CURRENCY_BREADTH_MIN}
    curated_only = set(config.CURATED_COFACTORS) - auto

    print(f"\n[A] the frozen currency set: {len(currency)} metabolites")
    print(f"    reached by breadth >= {config.CURRENCY_BREADTH_MIN:<3}                           "
          f"{len(auto):>10,}")
    print(f"    added by the curated list only                   {len(curated_only):>10,}")
    print(f"    curated entries the threshold also reaches       "
          f"{len(set(config.CURATED_COFACTORS) & auto):>10,}")

    print("\n    by breadth:")
    print(f"    {'breadth':>7} {'degree':>7}  {'formula':<22} {'ikey14':<16} source")
    for k in sorted(currency, key=lambda x: -len(breadth.get(x, ()))):
        src = "breadth" if k in auto else "curated"
        if k in config.CURATED_COFACTORS and k in auto:
            src = "both"
        print(f"    {len(breadth.get(k, ())):>7} {degree.get(k, 0):>7,}  "
              f"{formula(k):<22} {k:<16} {src}")

    # ------------------------------------------- the superseded rule, for audit
    counts = sorted(degree.values(), reverse=True)
    n = len(counts)

    def pct_cut(p: float) -> int:
        idx = min(max(int(round((1 - p / 100) * n)), 0), n - 1)
        return counts[idx]

    print("\n[B] the SUPERSEDED degree-percentile rule, for comparison")
    for p in config.CURRENCY_PERCENTILE_SWEEP:
        thr = pct_cut(p)
        sel = {k for k, v in degree.items() if v >= thr}
        bad = sel & REF_INTERMEDIATE
        print(f"    degree >= {thr:<6,} (p{p})  -> {len(sel):>4} members, "
              f"{len(bad)} labelled intermediates wrongly included")
    old_thr = pct_cut(config.CURRENCY_DEGREE_PERCENTILE)
    old_bad = sorted({k for k, v in degree.items() if v >= old_thr} & REF_INTERMEDIATE)
    print(f"    the plan's original p99.5 rule would have removed these as nodes:")
    for k in old_bad:
        print(f"      {formula(k):<22} {k}  degree {degree[k]:>5,}  "
              f"breadth {len(breadth.get(k, ())):>3}")

    # ------------------------------------------------- reference performance
    def errs(sel: set[str]) -> tuple[int, int]:
        return len(sel & REF_INTERMEDIATE), len(REF_COFACTOR - sel)

    fp, fn = errs(currency)
    fp_auto, fn_auto = errs(auto)
    print("\n[C] performance on the 27 labelled references")
    print(f"    breadth threshold alone : false positives {fp_auto}, false negatives {fn_auto}")
    print(f"    threshold + curated list: false positives {fp}, false negatives {fn}")
    if fp:
        print("    FALSE POSITIVES (real metabolites wrongly called currency):")
        for k in sorted(currency & REF_INTERMEDIATE):
            print(f"      {formula(k):<22} {k}")
    if fn:
        print("    false negatives (cofactors left in the graph as hub nodes):")
        for k in sorted(REF_COFACTOR - currency):
            print(f"      {formula(k):<22} {k}")

    # ------------------------------------------------------ settled questions
    print("\n[D] PLAN.md appendix B question 1, settled by this measurement")
    for k, (f, why) in config.SETTLED_NOT_CURRENCY.items():
        in_set = k in currency
        flag = "  *** CONTRADICTED - it IS in the set ***" if in_set else ""
        print(f"    NOT currency: {f:<22} breadth {len(breadth.get(k, ())):>3}  {why}{flag}")

    # ----------------------------------------------------------- sensitivity
    print("\n[E] threshold sensitivity (59 is fitted to the references, not validated)")
    print(f"    {'thr':>4} {'n_auto':>7} {'FP':>3} {'FN':>3}   metabolites entering vs thr=59")
    sens = {}
    for thr in (40, 45, 50, 55, 59, 65, 70, 80, 90, 100):
        sel = {k for k, m in breadth.items() if len(m) >= thr}
        f_p, f_n = errs(sel | set(config.CURATED_COFACTORS))
        delta = sorted(sel - auto, key=lambda x: -len(breadth[x]))
        sens[thr] = {"n_auto": len(sel), "fp": f_p, "fn": f_n,
                     "extra_vs_59": [formula(k) for k in delta]}
        extra = ", ".join(formula(k) for k in delta[:6]) or "-"
        print(f"    {thr:>4} {len(sel):>7} {f_p:>3} {f_n:>3}   {extra}")

    # -------------------------------------------------------- list integrity
    print("\n[F] curated-list integrity (InChIKey still names the stated formula)")
    problems = chem.verify_curated_formulae(example)
    if problems:
        for p in problems:
            print(f"    {p}")
    else:
        print("    all curated entries verified against the data")
    hard = [p for p in problems if "the data says" in p]
    if hard:
        raise AssertionError(
            f"{len(hard)} curated cofactor key(s) no longer name the stated "
            "molecule; fix config.CURATED_COFACTORS before using the set"
        )

    out = {
        "experiment": "03_currency_audit",
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "plan_section": "6.3",
        "rule": f"breadth >= {config.CURRENCY_BREADTH_MIN} OR curated cofactor list",
        "decision_record": "docs/DECISIONS.md D1",
        "currency_ikey14": sorted(currency),
        "members": [
            {"ikey14": k, "formula": formula(k), "breadth": len(breadth.get(k, ())),
             "degree": degree.get(k, 0),
             "source": "both" if (k in auto and k in config.CURATED_COFACTORS)
                       else ("breadth" if k in auto else "curated")}
            for k in sorted(currency, key=lambda x: -len(breadth.get(x, ())))
        ],
        "n_currency": len(currency),
        "n_from_breadth": len(auto),
        "n_curated_only": len(curated_only),
        "reference_false_positives": fp,
        "reference_false_negatives": fn,
        "superseded_rule_false_positives": [formula(k) for k in old_bad],
        "threshold_sensitivity": {str(k): v for k, v in sens.items()},
        "curated_list_problems": problems,
        "n_reactions_counted": len(rxns),
        "n_distinct_metabolites": len(degree),
        "generic_participants_skipped": n_generic,
        "unparsed_participants": n_unparsed,
        "runtime_sec": round(time.time() - t0, 1),
    }
    chem.CURRENCY_CACHE.parent.mkdir(parents=True, exist_ok=True)
    chem.CURRENCY_CACHE.write_text(json.dumps(out, indent=2))
    card = config.RESULTS / "03_currency_audit.json"
    card.write_text(json.dumps(out, indent=2))
    print(f"\nfrozen set: {chem.CURRENCY_CACHE}")
    print(f"result card: {card}")
    print(f"runtime: {out['runtime_sec']}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
