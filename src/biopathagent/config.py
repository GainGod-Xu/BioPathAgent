"""Central paths and constants. Import this; never hard-code a path elsewhere."""
from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------- locations
# Code lives in $HOME; data and results live on netscratch (group quota ~50 TiB).
CODE_ROOT = Path(__file__).resolve().parents[2]

PROJECT_ROOT = Path(
    os.environ.get(
        "BIOPATHAGENT_ROOT",
        "/n/netscratch/kempner_mzitnik_lab/Lab/hxu66/biopathagent",
    )
)

DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_INTERIM = PROJECT_ROOT / "data" / "interim"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DATA_SPLITS = PROJECT_ROOT / "data" / "splits"
RESULTS = PROJECT_ROOT / "results"
LOGS = PROJECT_ROOT / "logs"

_ALL_DIRS = [DATA_RAW, DATA_INTERIM, DATA_PROCESSED, DATA_SPLITS, RESULTS, LOGS]


def ensure_dirs() -> None:
    for d in _ALL_DIRS:
        d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------- sources
RHEA_TSV_BASE = "https://ftp.expasy.org/databases/rhea/tsv"

# Only files the project actually uses. Keep this list honest: if it is here, some
# module reads it; if nothing reads it, delete the line rather than download it.
RHEA_FILES = {
    "rhea-reaction-smiles.tsv": "RHEA:id -> reaction SMILES (the backbone)",
    "rhea-directions.tsv": "master id -> LR / RL / BI direction variants",
    "rhea2ec.tsv": "RHEA:id -> EC number",
    "rhea2kegg_reaction.tsv": "RHEA:id -> KEGG R-number (join to KEGG pathways)",
    "rhea2metacyc.tsv": "RHEA:id -> MetaCyc reaction id (no MetaCyc licence needed for ids)",
    "rhea2uniprot_sprot.tsv": "RHEA:id -> reviewed UniProt accession",
    "rhea-chebi-smiles.tsv": "ChEBI id -> SMILES, for participant-level work",
}

KEGG_REST = "https://rest.kegg.jp"

# KEGG REST is free for academic interactive use but is NOT a bulk-download endpoint.
# Every call here is a single whole-database LINK query (a handful of requests in total),
# which is what the service is designed to serve. Do not loop per-entry over it.
KEGG_LINKS = {
    "kegg_pathway_reaction.tsv": f"{KEGG_REST}/link/reaction/pathway",
    # NB: KEGG's database name for EC is "enzyme"; /link/ec/reaction returns HTTP 400.
    "kegg_reaction_ec.tsv": f"{KEGG_REST}/link/enzyme/reaction",
    "kegg_pathway_list.tsv": f"{KEGG_REST}/list/pathway",
    "kegg_reaction_list.tsv": f"{KEGG_REST}/list/reaction",
}

# RetroRules: released as a versioned dump. Resolved at fetch time rather than
# pinned blindly, because the project has moved hosts before.
RETRORULES_CANDIDATES = [
    "https://retrorules.org/dl/preparsed/rr02/rp2/hs",
    "https://zenodo.org/records/5827427/files/retrorules_rr02_rp2_hs.tar.gz",
]

# ---------------------------------------------------------------- chemistry
# Currency metabolites: frozen rule, PLAN.md 6.3. Removed as graph NODES everywhere;
# still available as verifier FEATURES. One list, one decision, not revisited per
# experiment. Identified by InChIKey-14 (skeleton block) so protonation and charge
# states collapse together.
#
# Each entry carries its justification. Adding a row is a design change: record it in
# docs/DECISIONS.md.
CURRENCY_NAMES = {
    "water": "universal solvent/participant",
    "proton": "charge balance only",
    "ATP": "phosphoryl donor, connects unrelated chemistry",
    "ADP": "counterpart of ATP",
    "AMP": "counterpart of ATP",
    "NAD(+)": "redox cofactor",
    "NADH": "redox cofactor",
    "NADP(+)": "redox cofactor",
    "NADPH": "redox cofactor",
    "CoA": "acyl carrier",
    "phosphate": "ubiquitous leaving group",
    "diphosphate": "ubiquitous leaving group",
    "CO2": "ubiquitous product",
    "O2": "ubiquitous oxidant",
    "ammonia": "ubiquitous N donor/acceptor",
    "ammonium": "ubiquitous N donor/acceptor",
    "GTP": "phosphoryl donor",
    "GDP": "counterpart of GTP",
    "FAD": "redox cofactor",
    "FADH2": "redox cofactor",
    "S-adenosyl-L-methionine": "methyl donor",
    "S-adenosyl-L-homocysteine": "counterpart of SAM",
    "hydrogen peroxide": "ubiquitous oxidant/product",
    "sulfate": "ubiquitous",
    "hydrogencarbonate": "ubiquitous",
}

# Degree-based rule that supplements the curated list (PLAN.md 6.3).
CURRENCY_DEGREE_PERCENTILE = 99.5
CURRENCY_PERCENTILE_SWEEP = (99.0, 99.5, 99.9)  # reported once, then frozen

# Open question, must be settled before Phase 3 (PLAN.md appendix B):
# acetyl-CoA is context-dependent and is deliberately NOT in the list above yet.
UNDECIDED_CURRENCY = {"acetyl-CoA", "L-glutamate", "2-oxoglutarate"}


# ---------------------------------------------------------------- gates
# Pre-registered in PLAN.md 3.1. Written down before the measurement exists so the
# project can fail honestly.
GATE1_PASS = 0.50
GATE1_DEGRADE = 0.20


def gate1_verdict(coverage: float) -> str:
    if coverage >= GATE1_PASS:
        return "PASS"
    if coverage >= GATE1_DEGRADE:
        return "DEGRADE-SCOPE"
    return "KILL"
