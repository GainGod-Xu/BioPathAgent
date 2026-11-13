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
# Verified 2026-11-13 via the Zenodo API. retrorules.org/dl/... and the record id
# 5827427 that circulate in older workflows both return 404 / the wrong release;
# resolve by record id instead of by the project's own host.
#   5828017  rr02-rp2-hs   70 MB  pre-parsed rules, multiple diameters  <- default
#   5827977  rr02-rp3-hs   61 MB
#   5827969  rr02-rp3-nohs 43 MB
#   18620886 RetroRules v1 SQLite, 2.4 GB, full rule set with metadata
RETRORULES_CANDIDATES = [
    "https://zenodo.org/api/records/5828017/files/retrorules_rr02_rp2_hs.tar.gz/content",
]
RETRORULES_SQLITE = (
    "https://zenodo.org/api/records/18620886/files/retrorules_2019_rr01_sqlite.zip/content"
)

# ---------------------------------------------------------------- chemistry
# Currency metabolites: PLAN.md 6.3. Removed as graph NODES everywhere; still
# available as verifier FEATURES. Identified by InChIKey-14 (skeleton block) so
# protonation and charge states collapse together.
#
# REVISED 2026-11-13 after experiments 03 and 04. See docs/DECISIONS.md D1.
# The original rule was "curated names OR reaction-degree above the 99.5th
# percentile". Experiment 03 showed the degree half selects glucose, pyruvate,
# 2-oxoglutarate, fumarate, acetate and L-alanine - removing those as nodes deletes
# glycolysis and the TCA cycle. Experiment 04 tested pathway breadth as a
# replacement; it does not separate the two classes either (reference cofactors
# reach down to breadth 21, reference intermediates up to 58).
#
# So the rule is now:
#     currency(m) = breadth(m) >= CURRENCY_BREADTH_MIN  OR  m in CURATED_COFACTORS
#
# with the threshold set for ZERO FALSE POSITIVES rather than best accuracy. The
# two errors are not symmetric: deleting pyruvate corrupts every pathway result
# downstream, while leaving FAD in the graph only adds a hub node. The curated list
# exists to catch the cofactors the threshold cannot reach, which is the job a
# curated list is actually for.
#
# The curated list is keyed by InChIKey-14 because that is what the code joins on.
# The formula is carried alongside and ASSERTED at load time - a key that no longer
# refers to the molecule named here is a silent correctness failure, so it is
# checked rather than trusted.

# Smallest threshold with no reference intermediate above it (exp 04: max was 58,
# 2-oxoglutarate). Calibrated on 27 hand-labelled references, NOT cross-validated.
# Sensitivity over 40..100 is reported by experiment 04.
CURRENCY_BREADTH_MIN = 59

# ikey14 -> (formula as it appears in Rhea, justification)
CURATED_COFACTORS = {
    # --- reached by the breadth threshold anyway; listed so the set is readable
    "GPRLSGONYQIRFK": ("H+", "charge balance only"),
    "XLYOFNOQVPJJNP": ("H2O", "universal solvent/participant"),
    "MYMOFIZGZYHOMD": ("O2", "ubiquitous oxidant"),
    "CURLTUGMZLYLDI": ("CO2", "ubiquitous one-carbon product"),
    "ZKHQWZAMYRWXGA": ("C10H12N5O13P3-4", "ATP, phosphoryl donor"),
    "XTWYTFMLZFPYCI": ("C10H12N5O10P2-3", "ADP, counterpart of ATP"),
    "BAWFJGJZGIEFAR": ("C21H26N7O14P2-", "NAD+, redox cofactor"),
    "BOPGDPNILDQYTO": ("C21H27N7O14P2-2", "NADH, redox cofactor"),
    "XJLXINKUBYWONI": ("C21H25N7O17P3-3", "NADP+, redox cofactor"),
    "ACFIXJIJDZMPPO": ("C21H26N7O17P3-4", "NADPH, redox cofactor"),
    "RGJOEKWQDUBAIZ": ("C21H32N7O16P3S-4", "free CoA, acyl carrier"),
    "NBIIXXVUZAFLBC": ("HO4P-2", "phosphate, ubiquitous leaving group"),
    "XPPKVPWEQAFLFU": ("HO7P2-3", "diphosphate, ubiquitous leaving group"),
    # --- below the threshold: the reason the curated list is needed at all
    "UDMBCSSLTHHNCD": ("C10H12N5O7P-2", "AMP, counterpart of ATP (breadth 56)"),
    "MEFKEPWMEQBLKI": ("C15H23N6O5S+", "SAM, methyl donor (breadth 53)"),
    "ZJUKTBDSGOFHSH": ("C14H20N6O5S", "SAH, counterpart of SAM (breadth 47)"),
    "ANKZYBDXHMZBDK": ("C17H18N4O9P-3", "FMN, redox cofactor (breadth 49)"),
    "YTNIXZGTHTVJBW": ("C17H21N4O9P-2", "FMNH2, redox cofactor (breadth 49)"),
    "IMGVNJNCCGXBHD": ("C27H30N9O15P2-3", "FAD, redox cofactor (breadth 21)"),
    "YPZRHBJKEMOYQH": ("C27H33N9O15P2-2", "FADH2, redox cofactor (breadth 21)"),
    "QGZKDVFQNNGYKY": ("H4N+", "ammonium, ubiquitous N donor/acceptor (breadth 57)"),
    "MHAJPDPJQMAIIY": ("H2O2", "ubiquitous oxidant/product (breadth 49)"),
    "XCCTYIAWTASOJW": ("C9H11N2O12P2-3", "UDP, glycosyl carrier counterpart (breadth 40)"),
    "XKMLYUALXHKNFT": ("C10H12N5O14P3-4", "GTP, phosphoryl donor"),
    "QGWNDRXFNXRZMB": ("C10H12N5O11P2-3", "GDP, counterpart of GTP"),
    "MAGIRAZQQVQNKP": ("Fe2S2+", "iron-sulfur cluster, protein-bound electron carrier"),
    "XSOVBBGAMBLACL": ("Fe2S2+2", "iron-sulfur cluster, oxidised form"),
    "CWYNVVGOOAEACU": ("Fe+2", "electron carrier"),
    "VTLYFUHAOXGGBS": ("Fe+3", "electron carrier"),
    "FKNQFGJONOIPTF": ("Na+", "counter-ion, carries no chemistry"),
    "VEXZGXHMUGYJMC": ("Cl-", "counter-ion, carries no chemistry"),
}

# SETTLED 2026-11-13 by experiments 03+04 (PLAN.md appendix B question 1).
# All three sit below the breadth threshold and all three carry the transferred
# carbon skeleton rather than a transferable group, so all three are NODES, not
# currency. Recorded here because the plan listed them as deliberately open and a
# later reader needs to see they were closed by measurement, not by taste.
SETTLED_NOT_CURRENCY = {
    "ZSLZBFCDCINBPY": ("C23H34N7O17P3S-4", "acetyl-CoA: breadth 55; the acetyl is "
                                           "the transferred skeleton, not a carrier"),
    "WHUUTDBJXJRKMK": ("C5H8NO4-", "L-glutamate: breadth 57; amino-acid metabolite"),
    "KPGXRSRHYNQIFN": ("C5H4O5-2", "2-oxoglutarate: breadth 58; TCA intermediate"),
    "LTYOQGRJFJAKNA": ("C24H33N7O19P3S-5", "malonyl-CoA: breadth 24; carbon donor"),
    "BDAGIHXWWSANSR": ("CHO2-", "formate: breadth 27; C1 metabolite"),
}

# Retained for the record: the superseded rule's parameters. Experiment 03 still
# reports the sweep so the revision is auditable against what it replaced.
CURRENCY_DEGREE_PERCENTILE = 99.5
CURRENCY_PERCENTILE_SWEEP = (99.0, 99.5, 99.9)


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
