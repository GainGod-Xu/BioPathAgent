#!/usr/bin/env python3
"""Experiment 04 - is reaction-degree the right currency criterion? (PLAN.md 6.3)

Why this experiment exists
--------------------------
Experiment 03 applied the frozen rule - currency = curated list OR reaction-degree
above the 99.5th percentile - and the degree half recovered the cofactors
correctly at the top (H+, H2O, O2, CoA, ATP, NAD(P)(H), Pi, PPi, CO2, SAM, ADP,
AMP, FMN, SAH, NH4+, H2O2 are ranks 1-21). Below that it also selected:

    C6H12O6    glucose          degree 239
    C3H3O3-    pyruvate         degree 221
    C5H4O5-2   2-oxoglutarate   degree 474
    C4H4O4-2   fumarate         degree 374
    C2H3O2-    acetate          degree 191
    C3H7NO2    L-alanine        degree 156

Removing those as graph nodes deletes glycolysis and the TCA cycle. So the rule as
written is wrong, and it is wrong in a way that would have corrupted every pathway
result downstream rather than failing loudly. This is the measurement that says so.

The hypothesis under test
-------------------------
Degree conflates two different things. A cofactor is used in many reactions
*spread across the whole of metabolism*. A central intermediate is used in many
reactions *concentrated in a few pathways*. So the discriminator is not how many
reactions a metabolite appears in, but how many distinct pathways:

    breadth(m) = number of distinct ordinary KEGG maps containing a reaction in
                 which m participates

Prediction: ATP and glucose have comparable degree but very different breadth. If
that separation exists, breadth replaces degree in the rule. If it does not, the
rule needs hand curation instead, and we will say so.

Either outcome is recorded in docs/DECISIONS.md. This is not a free parameter to
tune later - it is measured once, decided once, and frozen.

Run:  python experiments/04_currency_discriminator.py
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rdkit.Chem import rdMolDescriptors  # noqa: E402

from biopathagent import chem, config, data  # noqa: E402

# Reference metabolites with known status, used to read the separation off the
# measurement instead of asserting it. Identified by InChIKey-14 as printed by
# experiment 03, with the formula kept alongside as a check that the key still
# refers to what we think it does.
KNOWN_COFACTOR = {
    "GPRLSGONYQIRFK": "H+",
    "XLYOFNOQVPJJNP": "H2O",
    "ZKHQWZAMYRWXGA": "ATP C10H12N5O13P3-4",
    "XTWYTFMLZFPYCI": "ADP C10H12N5O10P2-3",
    "UDMBCSSLTHHNCD": "AMP C10H12N5O7P-2",
    "BAWFJGJZGIEFAR": "NAD+ C21H26N7O14P2-",
    "BOPGDPNILDQYTO": "NADH C21H27N7O14P2-2",
    "XJLXINKUBYWONI": "NADP+ C21H25N7O17P3-3",
    "ACFIXJIJDZMPPO": "NADPH C21H26N7O17P3-4",
    "RGJOEKWQDUBAIZ": "CoA C21H32N7O16P3S-4",
    "NBIIXXVUZAFLBC": "phosphate HO4P-2",
    "XPPKVPWEQAFLFU": "diphosphate HO7P2-3",
    "MEFKEPWMEQBLKI": "SAM C15H23N6O5S+",
    "ZJUKTBDSGOFHSH": "SAH C14H20N6O5S",
    "CURLTUGMZLYLDI": "CO2",
    "MYMOFIZGZYHOMD": "O2",
    "IMGVNJNCCGXBHD": "FAD C27H30N9O15P2-3",
    "ANKZYBDXHMZBDK": "FMN C17H18N4O9P-3",
}
KNOWN_INTERMEDIATE = {
    "WQZGKKKJIJFFOK": "glucose C6H12O6",
    "LCTONWCANYUPML": "pyruvate C3H3O3-",
    "KPGXRSRHYNQIFN": "2-oxoglutarate C5H4O5-2",
    "KDYFGRWQOYBRFD": "fumarate C4H4O4-2",
    "QTBSBXVTEAMEQO": "acetate C2H3O2-",
    "QNAYBMKLOCPYGJ": "L-alanine C3H7NO2",
    "WHUUTDBJXJRKMK": "L-glutamate C5H8NO4-",
    "ZSLZBFCDCINBPY": "acetyl-CoA C23H34N7O17P3S-4",
    "ZQPPMHVWECSIRJ": "oleate C18H33O2-",
}


def main() -> int:
    config.ensure_dirs()
    t0 = time.time()
    print("=" * 86)
    print("EXPERIMENT 04 - CURRENCY DISCRIMINATOR: degree vs pathway breadth")
    print("=" * 86)

    # ---------------------------------------------------------------- degree
    rxns = data.lr_reaction_smiles()
    degree: Counter = Counter()
    participants: dict[str, set[str]] = defaultdict(set)  # master -> ikeys
    for master, rxn in rxns.items():
        subs, prods = chem.split_reaction(rxn)
        seen = set()
        for smi in subs + prods:
            if chem.is_generic(smi):
                continue
            k = chem.ikey14(smi)
            if k:
                seen.add(k)
        participants[master] = seen
        for k in seen:
            degree[k] += 1

    # ------------------------------------------------------------- breadth
    # KEGG pathway -> reactions -> rhea master -> participants.
    cands = data.candidate_set()
    breadth: dict[str, set[str]] = defaultdict(set)
    for kegg_rid, rec in cands.items():
        parts = participants.get(rec["rhea_master"], set())
        for p in rec["pathways"]:
            for k in parts:
                breadth[k].add(p)
    n_maps = len({p for r in cands.values() for p in r["pathways"]})

    print(f"\nreactions (one LR variant per master)                {len(rxns):>8,}")
    print(f"drawable KEGG pathway reactions (gate denominator)   {len(cands):>8,}")
    print(f"ordinary KEGG maps they span                         {n_maps:>8,}")
    print(f"distinct metabolites with a degree                   {len(degree):>8,}")
    print(f"distinct metabolites with a pathway breadth          {len(breadth):>8,}")

    # ---------------------------------------------------- read the separation
    def row(k: str, label: str) -> dict:
        return {
            "ikey14": k,
            "label": label,
            "degree": degree.get(k, 0),
            "breadth": len(breadth.get(k, ())),
            "breadth_frac": len(breadth.get(k, ())) / max(n_maps, 1),
        }

    cof = sorted((row(k, v) for k, v in KNOWN_COFACTOR.items()),
                 key=lambda r: -r["breadth"])
    inter = sorted((row(k, v) for k, v in KNOWN_INTERMEDIATE.items()),
                   key=lambda r: -r["breadth"])

    print("\n[A] known COFACTORS (should be currency)")
    print(f"    {'degree':>7}  {'breadth':>7}  {'% maps':>7}   label")
    for r in cof:
        print(f"    {r['degree']:>7,}  {r['breadth']:>7}  {r['breadth_frac']:>6.1%}   {r['label']}")

    print("\n[B] known PATHWAY INTERMEDIATES (must NOT be currency)")
    print(f"    {'degree':>7}  {'breadth':>7}  {'% maps':>7}   label")
    for r in inter:
        print(f"    {r['degree']:>7,}  {r['breadth']:>7}  {r['breadth_frac']:>6.1%}   {r['label']}")

    cof_b = [r["breadth"] for r in cof]
    int_b = [r["breadth"] for r in inter]
    cof_d = [r["degree"] for r in cof]
    int_d = [r["degree"] for r in inter]
    print("\n[C] does either criterion separate the two reference groups?")
    print(f"    degree : cofactors min {min(cof_d):>6,}   intermediates max {max(int_d):>6,}"
          f"   {'SEPARATES' if min(cof_d) > max(int_d) else 'OVERLAPS'}")
    print(f"    breadth: cofactors min {min(cof_b):>6,}   intermediates max {max(int_b):>6,}"
          f"   {'SEPARATES' if min(cof_b) > max(int_b) else 'OVERLAPS'}")

    # A single threshold that is right on the references, if one exists.
    sep_breadth = min(cof_b) > max(int_b)
    threshold = None
    if sep_breadth:
        threshold = (max(int_b) + min(cof_b)) // 2 or min(cof_b)
        print(f"    -> a breadth threshold of {threshold} classifies all "
              f"{len(cof)+len(inter)} references correctly")

    # ------------------------------------------- what the candidate set looks like
    print("\n[D] metabolites ranked by pathway breadth (top 40)")
    ranked = sorted(breadth.items(), key=lambda kv: -len(kv[1]))
    example: dict[str, str] = {}
    for rxn in rxns.values():
        for smi in sum(chem.split_reaction(rxn), []):
            if not chem.is_generic(smi):
                k = chem.ikey14(smi)
                if k and k not in example:
                    example[k] = smi
    print(f"    {'rank':>4} {'breadth':>7} {'degree':>7}  {'formula':<22} {'ikey14':<16} note")
    rows_out = []
    for i, (k, maps) in enumerate(ranked[:200], 1):
        m = chem.parse(example.get(k, ""))
        formula = rdMolDescriptors.CalcMolFormula(m) if m else "?"
        note = KNOWN_COFACTOR.get(k, "") or KNOWN_INTERMEDIATE.get(k, "")
        tag = ""
        if k in KNOWN_COFACTOR:
            tag = "  <cofactor ref>"
        elif k in KNOWN_INTERMEDIATE:
            tag = "  <INTERMEDIATE ref>"
        rows_out.append({"rank": i, "ikey14": k, "formula": formula,
                         "breadth": len(maps), "degree": degree.get(k, 0),
                         "ref": note})
        if i <= 40:
            print(f"    {i:>4} {len(maps):>7} {degree.get(k,0):>7,}  {formula:<22} "
                  f"{k:<16}{tag}")

    out = {
        "experiment": "04_currency_discriminator",
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "plan_section": "6.3",
        "question": "does pathway breadth separate cofactors from central intermediates "
                    "where reaction degree does not?",
        "n_reactions": len(rxns),
        "n_candidate_reactions": len(cands),
        "n_ordinary_maps": n_maps,
        "reference_cofactors": cof,
        "reference_intermediates": inter,
        "degree_separates": bool(min(cof_d) > max(int_d)),
        "breadth_separates": bool(sep_breadth),
        "breadth_threshold_if_separating": threshold,
        "top200_by_breadth": rows_out,
        "runtime_sec": round(time.time() - t0, 1),
    }
    card = config.RESULTS / "04_currency_discriminator.json"
    card.write_text(json.dumps(out, indent=2))
    print(f"\nresult card: {card}")
    print(f"runtime: {out['runtime_sec']}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
