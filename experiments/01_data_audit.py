#!/usr/bin/env python3
"""Experiment 01 - data audit (PLAN.md phase 1).

Answers, with numbers rather than assertions:
  1. How big is each source, after parsing?
  2. How much of the KEGG pathway graph can we attach real chemistry to?
     KEGG pathway -> KEGG reaction -> Rhea reaction -> reaction SMILES
  3. What is the cross-database overlap (Rhea / KEGG / MetaCyc)?
  4. How many pathway edges survive to be candidates for Gate 1?

This is a precondition for Gate 1: template coverage can only be measured over
reactions we can actually draw. Writes a JSON result card to RESULTS/.

Run:  python experiments/01_data_audit.py
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from biopathagent import config  # noqa: E402

RAW = config.DATA_RAW


# ------------------------------------------------------------------ loaders
def load_reaction_smiles() -> dict[str, str]:
    """rhea-reaction-smiles.tsv has NO header: <rhea_id>\\t<reaction SMILES>."""
    out: dict[str, str] = {}
    path = RAW / "rhea" / "rhea-reaction-smiles.tsv"
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[0].isdigit():
                out[parts[0]] = parts[1]
    return out


def load_rhea_xref(filename: str) -> list[dict[str, str]]:
    """rhea2*.tsv: RHEA_ID, DIRECTION, MASTER_ID, ID (tab separated, with header)."""
    rows: list[dict[str, str]] = []
    path = RAW / "rhea" / filename
    with path.open() as fh:
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
    matches, which is how this was first written.
    """
    out: dict[str, dict[str, str]] = {}
    path = RAW / "rhea" / "rhea-directions.tsv"
    with path.open() as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 4 and p[0].isdigit():
                out[p[0]] = {"LR": p[1], "RL": p[2], "BI": p[3]}
    return out


def load_kegg_links(filename: str) -> list[tuple[str, str]]:
    path = RAW / "kegg" / filename
    out: list[tuple[str, str]] = []
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                out.append((parts[0], parts[1]))
    return out


# ------------------------------------------------------------------ helpers
def strip_prefix(x: str) -> str:
    return x.split(":", 1)[1] if ":" in x else x


def is_global_map(pathway_id: str) -> bool:
    """KEGG 'global/overview' maps (map011xx, map012xx) are aggregates of many
    pathways - e.g. map01100 'Metabolic pathways' holds thousands of reactions.

    Counting them alongside ordinary maps double-counts most of metabolism and
    inflates every per-pathway statistic, so they are excluded and reported
    separately. This is a frozen choice (PLAN.md 6.2), not a tunable.
    """
    n = strip_prefix(pathway_id).replace("map", "")
    return n.startswith("011") or n.startswith("012")


def split_reaction_smiles(rxn: str) -> tuple[list[str], list[str]] | None:
    if ">>" not in rxn:
        return None
    lhs, rhs = rxn.split(">>", 1)
    return lhs.split("."), rhs.split(".")


# ------------------------------------------------------------------ audit
def main() -> int:
    config.ensure_dirs()
    t0 = time.time()
    card: dict = {
        "experiment": "01_data_audit",
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "plan_section": "phase 1",
    }

    print("=" * 74)
    print("EXPERIMENT 01 - DATA AUDIT")
    print("=" * 74)

    # --- 1. source sizes -------------------------------------------------
    smiles = load_reaction_smiles()
    rhea2kegg = load_rhea_xref("rhea2kegg_reaction.tsv")
    rhea2ec = load_rhea_xref("rhea2ec.tsv")
    rhea2metacyc = load_rhea_xref("rhea2metacyc.tsv")
    rhea2uniprot = load_rhea_xref("rhea2uniprot_sprot.tsv")

    kegg_path_rxn = load_kegg_links("kegg_pathway_reaction.tsv")
    kegg_rxn_ec = load_kegg_links("kegg_reaction_ec.tsv")

    # Rhea ids in the SMILES file are DIRECTIONAL variants; MASTER_ID is the
    # undirected reaction. Both counts matter and are reported separately.
    masters_with_smiles = {r["MASTER_ID"] for r in rhea2kegg if r["RHEA_ID"] in smiles}

    sizes = {
        "rhea_directional_reactions_with_smiles": len(smiles),
        "rhea_master_reactions_xref_kegg": len({r["MASTER_ID"] for r in rhea2kegg}),
        "rhea_master_reactions_xref_ec": len({r["MASTER_ID"] for r in rhea2ec}),
        "rhea_master_reactions_xref_metacyc": len({r["MASTER_ID"] for r in rhea2metacyc}),
        "rhea_uniprot_pairs": len(rhea2uniprot),
        "rhea_reviewed_proteins": len({r["ID"] for r in rhea2uniprot}),
        "kegg_pathway_reaction_links": len(kegg_path_rxn),
        "kegg_reaction_ec_links": len(kegg_rxn_ec),
    }
    print("\n[1] source sizes")
    for k, v in sizes.items():
        print(f"    {k:<48} {v:>10,}")
    card["source_sizes"] = sizes

    # --- 2. KEGG pathway graph ------------------------------------------
    pathways: dict[str, set[str]] = defaultdict(set)
    global_maps: dict[str, set[str]] = defaultdict(set)
    for p, r in kegg_path_rxn:
        (global_maps if is_global_map(p) else pathways)[strip_prefix(p)].add(strip_prefix(r))

    rxn_in_normal = set().union(*pathways.values()) if pathways else set()
    rxn_in_global = set().union(*global_maps.values()) if global_maps else set()

    print("\n[2] KEGG pathway graph (global/overview maps excluded)")
    print(f"    ordinary pathway maps                            {len(pathways):>10,}")
    print(f"    global/overview maps excluded                    {len(global_maps):>10,}")
    print(f"    distinct reactions in ordinary maps              {len(rxn_in_normal):>10,}")
    print(f"    distinct reactions only in global maps           "
          f"{len(rxn_in_global - rxn_in_normal):>10,}")
    sizes_pw = sorted(len(v) for v in pathways.values())
    if sizes_pw:
        mid = sizes_pw[len(sizes_pw) // 2]
        print(f"    reactions per map  min/median/max                "
              f"{sizes_pw[0]:>4,} /{mid:>5,} /{sizes_pw[-1]:>5,}")
    card["kegg_graph"] = {
        "ordinary_maps": len(pathways),
        "global_maps_excluded": len(global_maps),
        "reactions_in_ordinary_maps": len(rxn_in_normal),
        "reactions_only_in_global_maps": len(rxn_in_global - rxn_in_normal),
        "reactions_per_map_median": sizes_pw[len(sizes_pw) // 2] if sizes_pw else 0,
        "reactions_per_map_max": sizes_pw[-1] if sizes_pw else 0,
    }

    # --- 3. the join: can we draw these reactions? ------------------------
    directions = load_directions()
    card["source_sizes"]["rhea_direction_masters"] = len(directions)

    kegg2rhea_master: dict[str, set[str]] = defaultdict(set)
    kegg2rhea_directional: dict[str, set[str]] = defaultdict(set)
    for r in rhea2kegg:
        master = r["MASTER_ID"]
        kegg2rhea_master[r["ID"]].add(master)
        # go master -> LR variant -> SMILES (see load_directions docstring)
        lr = directions.get(master, {}).get("LR")
        if lr and lr in smiles:
            kegg2rhea_directional[r["ID"]].add(lr)

    mapped = {k for k in rxn_in_normal if k in kegg2rhea_master}
    drawable = {k for k in rxn_in_normal if kegg2rhea_directional.get(k)}

    n_total = len(rxn_in_normal)
    print("\n[3] chemistry availability for KEGG pathway reactions")
    print(f"    reactions in ordinary maps                       {n_total:>10,}")
    print(f"    ... with a Rhea cross-reference                  {len(mapped):>10,}"
          f"   ({len(mapped)/max(n_total,1):6.1%})")
    print(f"    ... with an actual reaction SMILES               {len(drawable):>10,}"
          f"   ({len(drawable)/max(n_total,1):6.1%})")
    card["chemistry_availability"] = {
        "kegg_reactions_in_ordinary_maps": n_total,
        "with_rhea_xref": len(mapped),
        "with_reaction_smiles": len(drawable),
        "frac_with_rhea_xref": len(mapped) / max(n_total, 1),
        "frac_with_smiles": len(drawable) / max(n_total, 1),
    }

    # per-pathway coverage: a pathway is only usable if most of it is drawable
    per_pw = []
    for pid, rxns in pathways.items():
        d = len(rxns & drawable)
        per_pw.append((pid, len(rxns), d, d / max(len(rxns), 1)))
    fully = [p for p in per_pw if p[3] == 1.0]
    mostly = [p for p in per_pw if p[3] >= 0.8]
    half = [p for p in per_pw if p[3] >= 0.5]
    print(f"    maps 100% drawable                               {len(fully):>10,}"
          f"   ({len(fully)/max(len(per_pw),1):6.1%})")
    print(f"    maps >=80% drawable                              {len(mostly):>10,}"
          f"   ({len(mostly)/max(len(per_pw),1):6.1%})")
    print(f"    maps >=50% drawable                              {len(half):>10,}"
          f"   ({len(half)/max(len(per_pw),1):6.1%})")
    card["pathway_drawability"] = {
        "maps_total": len(per_pw),
        "maps_fully_drawable": len(fully),
        "maps_ge80pct": len(mostly),
        "maps_ge50pct": len(half),
    }

    # --- 4. cross-database overlap --------------------------------------
    m_kegg = {r["MASTER_ID"] for r in rhea2kegg}
    m_ec = {r["MASTER_ID"] for r in rhea2ec}
    m_mcyc = {r["MASTER_ID"] for r in rhea2metacyc}
    m_prot = {r["MASTER_ID"] for r in rhea2uniprot}
    print("\n[4] cross-database overlap, by Rhea MASTER id")
    print(f"    Rhea masters with KEGG xref                      {len(m_kegg):>10,}")
    print(f"    Rhea masters with MetaCyc xref                   {len(m_mcyc):>10,}")
    print(f"    KEGG and MetaCyc both                            {len(m_kegg & m_mcyc):>10,}")
    print(f"    KEGG only                                        {len(m_kegg - m_mcyc):>10,}")
    print(f"    MetaCyc only                                     {len(m_mcyc - m_kegg):>10,}")
    print(f"    with EC                                          {len(m_ec):>10,}")
    print(f"    with reviewed protein                            {len(m_prot):>10,}")
    print(f"    KEGG-linked AND EC AND protein                   "
          f"{len(m_kegg & m_ec & m_prot):>10,}")
    card["overlap"] = {
        "rhea_masters_kegg": len(m_kegg),
        "rhea_masters_metacyc": len(m_mcyc),
        "kegg_and_metacyc": len(m_kegg & m_mcyc),
        "kegg_only": len(m_kegg - m_mcyc),
        "metacyc_only": len(m_mcyc - m_kegg),
        "with_ec": len(m_ec),
        "with_reviewed_protein": len(m_prot),
        "kegg_ec_protein": len(m_kegg & m_ec & m_prot),
    }

    # --- 5. what Gate 1 will actually run on ------------------------------
    # A Gate 1 candidate needs: in an ordinary pathway map, has reaction SMILES,
    # and parses into reactants/products.
    parseable = 0
    multi = 0
    for k in drawable:
        rid = sorted(kegg2rhea_directional[k])[0]
        parts = split_reaction_smiles(smiles[rid])
        if parts:
            parseable += 1
            if len(parts[0]) > 1 or len(parts[1]) > 1:
                multi += 1
    print("\n[5] Gate 1 candidate set")
    print(f"    drawable pathway reactions                       {len(drawable):>10,}")
    print(f"    ... parse into reactants>>products               {parseable:>10,}")
    print(f"    ... multi-component (needs currency handling)    {multi:>10,}"
          f"   ({multi/max(parseable,1):6.1%})")
    card["gate1_candidates"] = {
        "drawable": len(drawable),
        "parseable": parseable,
        "multi_component": multi,
    }

    card["runtime_sec"] = round(time.time() - t0, 1)
    out = config.RESULTS / "01_data_audit.json"
    out.write_text(json.dumps(card, indent=2))
    print(f"\nresult card: {out}")
    print(f"runtime: {card['runtime_sec']}s")

    print("\nNOTE: this audit measures DATA availability only. It does not measure")
    print("      template coverage - that is Gate 1, experiment 02.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
