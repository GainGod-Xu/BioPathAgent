# BioPathAgent — Executable Research Plan (Revision 2)

**Evidence-Grounded and Self-Verifying Biochemical Pathway Reasoning**

> Core thesis: *chemical validity is not biochemical validity.*
> Pipeline: **Generate → Ground → Verify → Calibrate → Plan**

Revision 2 date: 2026-11-13. Supersedes the v1 PDF (`docs/BioPathAgent_Complete_Research_Plan.pdf`).
This file is authoritative; the PDF is kept for provenance only.

---

## 0. What Revision 2 changes and why

v1 is strategically sound — benchmark first, agents last, LLM optional, no RL until
classical search is proven to be the bottleneck. Revision 2 does not relitigate any of
that. It fixes the places where v1 is **not yet falsifiable or not yet executable**.

| # | Gap in v1 | Fix in Rev 2 | Section |
|---|---|---|---|
| G1 | Every results table is empty; no pre-registered numbers. A plan with no stated expected values cannot fail. | Pre-registered success / kill thresholds for every work package. | §3 |
| G2 | No null controls. The plan can produce a high Top-k that is entirely carried by template enumeration, with grounding and verification contributing nothing. | Mandatory **shuffled-grounding** and **retrieval-only** controls; a result is not reportable without them. | §4 |
| G3 | "Pathway-disjoint" / "family-disjoint" named but not operationalized. Rhea, KEGG, MetaCyc and RetroRules overlap heavily; a KEGG test edge may be literally a Rhea training row. | Concrete cross-database dedup + split construction protocol, frozen before any model is trained. | §5 |
| G4 | "Known pathway M0→M1→…→M4" is undefined. KEGG maps are branching hypergraphs, not chains. | Explicit pathway-graph extraction algorithm, with reporting of yield. | §6 |
| G5 | Currency metabolites offered three alternative treatments ("down-weighted, excluded, or prohibited"). | One frozen rule, stated numerically. | §6.3 |
| G6 | Thermodynamics listed as *optional*, though it is the cheapest non-ML biochemical evidence and the main thing separating "biochemical" from "chemical". | Promoted to a **core** evidence stream. | §7, §9 |
| G7 | Hard negatives are generated from the positives, so the verifier can learn the construction artifact instead of biochemistry. | Adversarial artifact audit; trivial-feature probe must fail. | §8.3 |
| G8 | Calibration set is drawn from curated hard negatives, but deployment inputs are template-enumerated search candidates. This is a train/deploy distribution mismatch. | Calibration is fit on **search-generated** candidates, not benchmark negatives. | §9.2 |
| G9 | `C(path) = -Σ log p_i + λ·h(target distance) + …` folds the heuristic into the cost. That is not A*. | Separate g and h; state admissibility; report as weighted A* (ε-admissible). | §10.2 |
| G10 | "General LLM as biochemical judge" has no protocol, so the headline no-LLM claim is unfalsifiable. | Fixed LLM baseline protocol with identical metrics including calibration. | §11 |
| G11 | No statistics: no seeds, no intervals, no paired tests. | Wilson intervals for proportions; paired bootstrap over pathways for ablations. | §12 |
| G12 | Baselines unnamed. | Nine peer-reviewed, PubMed-verified comparators. | §13 |
| G13 | No compute budget, no timeline, no archival plan. | Added, CPU-first, within current cluster constraints. | §14, §15, §16 |

---

## 0.5 Build order: a working system first, upgrades second

**Governing rule (Hao, 2026-11-13): first stand up something that actually solves the
problem; upgrade after.** Rev 2 adds a lot of rigor machinery. Not all of it belongs in the
first build, and this section says which parts do.

The test is not "is this rigorous?" — it is **"does omitting it now force us to redo the
work later?"**

### Tier A — must be right in V0 (retrofitting invalidates every number already produced)

| Item | Why it cannot wait |
|---|---|
| Gate 1 template coverage (§3.1) | Decides the scope of everything. Cheap, days. |
| Cross-database dedup `sig(R)` (§5.1) | If train/test share reactions, every number ever measured is void. |
| Pathway-disjoint split S5 (§5.2) | Built once with the dataset; retrofitting means re-running all experiments. |
| Currency-metabolite rule (§6.3) | Changes the graph itself, hence every path and every metric. |
| Pathway extraction algorithm (§6.2) | Defines what a "pathway" even is. |
| Controls C2 + C5 (§4) | ~50 lines, run in the same script. Without them V0's Top-k means nothing. |

Tier A is roughly one extra week on top of a naive build, and it is the difference between
a number and an artifact.

### Tier B — V0 ships without these; add as upgrades

Verifier and the whole evidence-group apparatus; hard-negative benchmark and its artifact
audit (§8.3); enzyme-reaction grounding and C1; calibration and abstention (§9.2);
uncertainty-weighted A* (§10.2 — V0 uses plain Beam/A* on template priors); missing-knowledge
benchmark; agent layer; LLM ablations (§11.3); full statistics protocol (§12 — V0 reports
point estimates, V1 adds intervals and paired bootstrap); splits S2/S3/S4.

### V0 definition of done

```
Rhea / MetaCyc / KEGG
  → pathway graph (§6.2, currency removed per §6.3)
  → dedup + pathway-disjoint split (§5)
  → mask 1–3 steps
  → RetroRules expansion → RDKit validation
  → Beam / A* search
  → Top-k pathway recovery + search cost
  → C2 retrieval baseline + C5 frequency prior, same script
```

No LLM, no learned grounding, no verifier, no agents. **If V0 clears §3.3's 1- and 2-step
thresholds, the project is real and upgrades are worth building. If it does not, no amount
of Tier B saves it** — and that is exactly why V0 comes first.

Upgrade order after V0, each gated by whether it moves a V0 number:
`verifier → calibration → uncertainty-aware search → grounding → missing-knowledge → agents → case study`

Grounding sits late **because of scientific risk, not GPU cost** (Hao, 2026-11-13: GPU use is
fine, just keep it sane). It is the component most likely to turn out cosmetic, and control
C1 is what decides that — so it is not worth building before there is a V0 number for it to
move.

One exception, pulled early: **ESM-2 embedding extraction is a one-off, cacheable ~20–60
GPU-h**. Run it as soon as the protein table exists (Phase 1), so the verifier can use
protein features in V1 without waiting on the contrastive model. Extracting embeddings is not
the same commitment as training the grounding model.

---

## 1. Project vision (unchanged from v1)

The system answers whether a biochemical path is merely **chemically constructible** or
**genuinely supported** by enzyme, reaction, database and thermodynamic evidence. The
defining capability is not free-form generation; it is self-verification and principled
abstention when evidence is weak.

Positioning, verbatim from v1 and retained: do **not** lead with "we built five AI agents."
The scientific contribution is the rigorous separation of candidate generation from
biochemical verification, with calibrated uncertainty propagated from individual reaction
steps to pathway-level decisions. The multi-agent layer is an integration shell.

---

## 2. Central research questions

| Q | Objective | Decided by |
|---|---|---|
| Q1 | Can known metabolic pathways be reconstructed from incomplete information? | WP1 |
| Q2 | Can biochemical validity be distinguished from chemical plausibility? | WP4 |
| Q3 | Can the system recognize biochemical hallucinations (wrong enzyme / substrate / product / cofactor / direction)? | WP3+WP4 |
| Q4 | Can calibrated uncertainty improve pathway search by pruning unreliable branches without destroying recall? | WP6+WP7 |
| Q5 | Can the system recover deliberately removed biochemical knowledge rather than retrieve memorized database edges? | WP8 |

Q5 is the gate between a *reconstruction* paper and a *discovery* paper. It is the only
question whose answer changes the target venue, so it must not be answered by a split that
leaks.

---

## 3. Pre-registered thresholds (G1)

These are committed **before** data are seen. They are written so the project can fail
honestly. "Kill" means: stop, report the negative result, and do not build the next layer
on top.

### 3.1 Gate 1 — template coverage (WP2, experiment 1)

| Quantity | Pass | Degrade scope | Kill |
|---|---|---|---|
| Fraction of **drawable** KEGG pathway reactions reproducible by ≥1 RetroRules template (any radius), substrate→product | ≥ 50% | 20–50% | < 20% |

**Denominator, fixed after experiment 01 (2026-11-13).** The gate is measured over
*drawable* reactions, not all KEGG pathway edges, because two different failures were
otherwise being added together: "we have no structure for this reaction" and "no template
reproduces this reaction". Only the second is what Gate 1 is asking about.

Measured data ceiling (experiment 01, `results/01_data_audit.json`):

| | count | of ordinary-map reactions |
|---|---|---|
| KEGG reactions in 166 ordinary pathway maps | 8,215 | 100% |
| with a Rhea cross-reference | 4,263 | 51.9% |
| **with an actual reaction SMILES → Gate 1 candidate set** | **4,186** | **51.0%** |
| of those, multi-component (currency handling mandatory) | 3,965 | 94.7% of candidates |

So template coverage is reported two ways, always together:
`coverage_conditional` (of 4,186 drawable) is the gate; `coverage_absolute`
(of 8,215, i.e. conditional × 0.510) is the honest headline for "how much of KEGG
metabolism can this system touch at all". Neither number is quoted alone.

That 94.7% multi-component figure also settles an open design question early: currency
handling is not a refinement, it is load-bearing for almost every reaction in the set.

Degrade = restrict the biochemical scope (central carbon metabolism, or xenobiotic/drug
metabolism only) and re-measure before proceeding. Kill = the deterministic engine cannot
carry the paper; the project must be re-scoped around the verifier alone.

*Rationale:* this is the single load-bearing number in the whole plan. If coverage is in the
single digits, every downstream Top-k is a measurement of the enumerator's blind spots, not
of reasoning. We have been burned by exactly this failure mode before (a 5–6% coverage
ceiling silently capping an entire pipeline) — so coverage is experiment #1, not #4.

### 3.2 Single-step recovery (WP2)

Given true substrate + true EC, recover the true product.

| Metric | Pass | Kill |
|---|---|---|
| Top-1 (exact product match, canonical InChIKey) | ≥ 35% | < 15% |
| Top-10 | ≥ 60% | < 30% |
| Top-10 on **EC-family-disjoint** split | ≥ 40% | < 20% |

### 3.3 Masked pathway recovery (WP1)

Exact full-path recovery, Top-10, pathway-disjoint split:

| Hidden steps | Pass | Kill |
|---|---|---|
| 1 | ≥ 50% | < 25% |
| 2 | ≥ 30% | < 12% |
| 3 | ≥ 15% | < 5% |
| 4 | ≥ 8% | — (reported, not gated) |

Graceful degradation is itself a claim: recovery must fall roughly geometrically with hidden
step count. A cliff (e.g. 50% → 3%) indicates the 1-step result is retrieval, not reasoning.

### 3.4 Verifier (WP4)

| Metric | Pass | Kill |
|---|---|---|
| AUROC vs **family-matched hard** negatives | ≥ 0.80 | < 0.65 |
| AUROC on the *hardest* single error type | ≥ 0.70 | < 0.55 |
| Gain over chemistry-only features (ΔAUROC) | ≥ +0.05 | ≤ +0.01 |
| Trivial-feature probe AUROC (§8.3) | ≤ 0.60 | > 0.75 → negatives are artifactual, rebuild |

The ΔAUROC row is the central hypothesis of the paper. If biological evidence adds ≤ 0.01
over RDKit-level chemistry, the thesis "chemical validity ≠ biochemical validity" is not
supported *by this operationalization*, and that is the honest result to report.

### 3.5 Calibration (WP6)

| Metric | Pass | Kill |
|---|---|---|
| ECE (15 equal-mass bins) on deployment-distribution candidates | ≤ 0.05 | > 0.15 |
| Accuracy gain at 50% coverage vs full coverage | ≥ +10 points | ≤ +2 points |

### 3.6 Search utility (WP7)

| Metric | Pass | Kill |
|---|---|---|
| Nodes expanded, verifier-pruned vs unpruned, at matched Top-10 path recall | ≤ 0.5× | ≥ 0.9× |
| Top-10 path recall lost to pruning | ≤ 3 points | ≥ 10 points |

### 3.7 Missing knowledge (WP8)

Discovery framing is unlocked **only** if Level 4 (reaction template removed) and Level 5
(entire reaction/enzyme family removed) both exceed:

| Level | Pass for discovery framing |
|---|---|
| 4 | Top-10 recovery ≥ 20% |
| 5 | Top-10 recovery ≥ 10%, and > 3× a frequency-prior baseline |

Below this, retain the trustworthy-reconstruction framing (v1 decision rule, kept).

---

## 4. Mandatory null controls (G2)

**No number from WP1–WP8 is reportable without its paired control.** These are cheap and
they are the difference between a result and an artifact.

### C1 — Shuffled grounding
Permute the protein sequence / EC annotation across reactions within the test set, holding
everything else fixed. Re-run end to end.

- *Expectation:* performance collapses toward the enzyme-agnostic template baseline.
- *If it does not collapse:* the grounding module is not being read, and Contribution 3
  must be demoted. Report it.

### C2 — Retrieval-only baseline
For a masked pathway, retrieve the most similar pathway in the training set (metabolite
fingerprint over endpoints) and copy its intermediates.

- *Expectation:* weak on pathway-disjoint splits.
- *If it is strong:* the benchmark measures database lookup, and the split is leaking.
  Tighten §5 and rebuild.

### C3 — Verifier-blind planner
Replace verifier probabilities with a constant. Re-run search.

- Isolates how much of the planning gain is the verifier versus the search algorithm.

### C4 — Chemistry-only verifier
Verifier trained on RDKit/template features only, no enzymatic / database / thermodynamic
evidence. This is the denominator of the ΔAUROC claim in §3.4.

### C5 — Frequency prior
Rank candidate products / enzymes by training-set frequency alone, ignoring the query.
Any benchmark where this is competitive is a popularity benchmark.

### C6 — Degenerate shortcut audit
For each recovered path, check whether it routes through a currency metabolite bridge
(§6.3). Paths that do are counted separately and reported, never silently included.

---

## 5. Leakage and split protocol (G3)

Frozen before any model training. Implemented once in `preprocessing/build_splits.py`,
hash-stamped, and never regenerated per-experiment.

### 5.1 Cross-database reaction identity

Databases overlap heavily. Define a canonical reaction signature:

```
sig(R) = ( sorted(InChIKey14 of non-currency reactants),
           sorted(InChIKey14 of non-currency products),
           canonical reaction-center SMARTS )
```

InChIKey-14 (skeleton block only) is used deliberately, so stereo/protonation variants of
the same transformation collapse together — the same dedup standard used elsewhere in the
group. Two reactions with equal `sig` are **the same reaction** regardless of source
database, and may never straddle a train/test boundary.

### 5.2 Split hierarchy

Splits are nested, each strictly harder, and all five are reported:

| Split | Construction | Purpose |
|---|---|---|
| S1 random | random over reactions | sanity only; **never the headline** |
| S2 sequence-disjoint | MMseqs2 cluster at 30% identity, cluster-level holdout | unseen proteins |
| S3 reaction-family-disjoint | hold out whole RetroRules template families (by radius-0 center) | unseen transformations |
| S4 EC-family-disjoint | hold out EC sub-subclasses (x.y.z.*) | unseen catalytic chemistry |
| S5 pathway-disjoint | hold out whole pathway graphs; all member reactions removed | the headline result |

### 5.3 Leakage assertions (automated, run in CI)

The build fails loudly if any of these is violated:

1. `sig(R)` intersection between train and test is empty, for every split.
2. No test protein has >30% identity to any train protein (S2–S5).
3. No test reaction's template family appears in train (S3–S5).
4. No test pathway shares *any* reaction with a train pathway (S5).
5. Currency-metabolite-only overlap does not count as a shared reaction.

### 5.4 Deliberate overlap report
Publish the measured Rhea∩KEGG∩MetaCyc∩EnzymeMap overlap matrix as a supplementary table.
Reviewers will ask. Having it ready is worth one afternoon.

---

## 6. Data foundation

Open resources only; no commercial data needed for the core paper.

| Resource | Role | Status |
|---|---|---|
| Rhea | curated reactions, directionality, EC/ChEBI | core |
| EnzymeMap | atom-mapped enzymatic reactions | core |
| RetroRules | SMARTS templates at multiple radii | core |
| Swiss-Prot / UniProt | reviewed sequences, EC, reaction links | core |
| MetaCyc | pathway graphs | core |
| KEGG | pathway/reaction evaluation | core |
| **eQuilibrator** | ΔG′°, reversibility, physiological direction | **core (promoted, G6)** |
| BRENDA / SABIO-RK | kinetic evidence | optional |

### 6.1 Unified schema (unchanged from v1)

```
reaction: reaction_id | ec_number | reactants | products | atom_mapped_reaction
          reaction_center | direction | reversible | cofactors | organism
          source | template_ids | template_radii | protein_ids | pathway_ids
          dG_prime | dG_uncertainty                            # added Rev 2
protein:  protein_id | sequence | ec_numbers | reaction_ids | organism | reviewed
pathway:  pathway_id | organism | reaction_ids | ec_numbers | metabolite_ids
          graph_edges | pathway_length | extraction_method     # added Rev 2
```

### 6.2 Pathway-graph extraction (G4)

KEGG/MetaCyc pathways are branching multi-substrate hypergraphs. "M0→M1→M2" must be
*derived*, not assumed. Frozen algorithm:

1. Build a directed metabolite graph: edge M→M′ iff a reaction in the pathway has M as a
   non-currency substrate and M′ as a non-currency product.
2. Delete currency metabolites as nodes entirely (§6.3).
3. Enumerate **maximal simple directed chains** of length 2–6 reactions.
4. Discard chains containing a node of in+out degree > 8 within the pathway (hub artifact).
5. Deduplicate chains by `sig` sequence.
6. Record `extraction_method` per chain for provenance.

**Report the yield**: number of pathways in, number of usable chains out, length histogram.
If yield is implausibly small, the benchmark is too thin and scope must widen before
anything is trained.

### 6.3 Currency metabolites — one frozen rule (G5)

v1 offered three options. Rev 2 picks one:

> A compound is **currency** iff it is in the curated list *or* its degree in the full
> reaction graph exceeds the 99.5th percentile. Currency metabolites are **removed as
> graph nodes** for all pathway extraction, search and evaluation. They remain available
> as *features* (cofactor identity is verifier evidence, §7) but can never form a graph
> edge or bridge.

Curated seed list: H₂O, H⁺, ATP/ADP/AMP, NAD(P)⁺/NAD(P)H, CoA, Pᵢ/PPᵢ, CO₂, O₂, NH₃/NH₄⁺,
GTP/GDP, FAD/FADH₂, SAM/SAH, glutamate/2-oxoglutarate (as amino donor/acceptor pair only),
acetyl-CoA (context-dependent — flagged, decided once, documented).

The exact list is frozen in `preprocessing/currency.py` with a one-line justification per
entry, and a sensitivity analysis at the 99.0 / 99.5 / 99.9 percentile is reported once.
It is not revisited per experiment.

### 6.4 Preprocessing (from v1, retained)
Canonicalize structures and identifiers; validate atom mapping and direction; retain
stereochemistry where possible; normalize EC/organism/cofactor/cross-references; dedup while
preserving provenance; construct pathway graphs; annotate currency metabolites.

---

## 7. WP1–WP2 — Benchmark and deterministic engine

### 7.1 Tasks (unchanged)

| Task | Input | Output |
|---|---|---|
| A | substrate (+ optional EC) | product, reaction, EC/enzyme |
| B | start + endpoint, 2 reactions apart | 1 intermediate + 2 reactions |
| C | start + endpoint, 3 apart | 2 intermediates + 3 reactions |
| D | start + endpoint, 4 apart | 3 intermediates + 4 reactions |
| E | current metabolite + context | ranked downstream branches |
| F | start + target | ranked complete biochemical paths |

### 7.2 Metrics
Exact product match; Top-k product recovery; Reaction Recall@k; EC Recall@k; intermediate
metabolite recall; exact pathway recovery; Top-1/5/10/50 pathway recovery; reaction/EC
overlap and graph edit distance (secondary); search cost (nodes expanded, edges generated,
reactions evaluated, runtime, peak memory).

All reported **per split S1–S5**, with Wilson 95% intervals (§12).

### 7.3 Deterministic engine

```
substrate → candidate RetroRules templates → substructure match → RDKit RunReactants
          → candidate products → sanitization → mass/charge/reaction-center checks
          → canonicalized candidates
```

**Confidence correction (v1 §8.1, retained and strengthened):** template radius must not be
read as probability. `confidence = radius/10` is a heuristic, not a calibrated estimate.
Radius enters the verifier as one feature among many, and the paper states this explicitly
as a correction to common practice.

### 7.4 Baseline ladder

| System | 1-step | 2-step | 3-step | 4-step |
|---|---|---|---|---|
| frequency prior (C5) | | | | |
| retrieval-only (C2) | | | | |
| RetroRules only | | | | |
| + RDKit validation | | | | |
| + database evidence | | | | |
| + enzyme constraint | | | | |

---

## 8. WP3 — Biochemical hallucination benchmark

Each true (enzyme, substrate, product) tuple generates multiple controlled negatives.

| Error type | Construction |
|---|---|
| wrong enzyme | same EC superclass, similar sequence, or related catalytic family |
| wrong substrate | structurally similar compound sharing relevant functional groups |
| wrong product | plausible alternative product, scaffold largely preserved |
| wrong cofactor | NADH↔NADPH, ATP↔GTP where biologically inappropriate |
| wrong direction | reverse a reaction with strong physiological directionality |
| mass imbalance | chemically invalid stoichiometry (easy control) |
| reaction-center corruption | alter the transformation center, keep molecules plausible |
| false pathway shortcut | connect metabolites separated by multiple real reactions |

Include easy, hard and family-matched negatives so that high performance cannot be explained
by trivial artifact detection.

### 8.3 Adversarial artifact audit (G7)

Because negatives are *derived from* positives, they may carry a signature unrelated to
biochemistry. Before the verifier is trained:

1. Train a logistic regression on **trivial features only** — heavy-atom count, molecular
   weight, formula difference, string lengths, number of rings.
2. Measure AUROC on the same splits.
3. **Requirement: ≤ 0.60.** Above 0.75, the negative set is artifactual and must be rebuilt
   with matched distributions.
4. Report this probe in the paper. It is a short paragraph that preempts the single most
   damaging reviewer objection.

Additionally, where available, include a small set of **externally sourced** negatives —
e.g. Rhea entries deprecated or corrected between releases — as an artifact-free holdout.

---

## 9. WP4–WP6 — Verifier and calibration

**Input:** (E, S, P, context) → **Output:** P(biochemically valid | E, S, P, context)

### 9.1 Evidence groups

| Group | Features |
|---|---|
| chemical | RDKit validity, atom/mass balance, charge balance, reaction-center consistency, stereochemistry |
| template | radius, specificity, number of matched templates, template frequency, provenance |
| enzymatic | EC compatibility, protein–reaction similarity, sequence/family support, substrate compatibility |
| database | Rhea/KEGG/MetaCyc/UniProt support, provenance count |
| thermodynamic | ΔG′°, uncertainty, reversibility, physiological direction, cofactor compatibility |

Model ladder: logistic regression → gradient-boosted trees → compact neural verifier, the
last only if the simpler models leave substantial headroom. Ablate **by evidence group**,
not just by model — the group ablation is the scientific result.

### 9.2 Calibration on the deployment distribution (G8)

This is the subtlest error available to this project, and an analogous train/sample mismatch
has already cost a sibling project months. Stated explicitly:

- The verifier is **trained** on curated hard negatives (§8).
- The verifier is **deployed** on template-enumerated search candidates, whose difficulty
  distribution is completely different — mostly easy chemical nonsense, occasionally subtle.
- Therefore **calibration is fit on a held-out set of search-generated candidates**, labelled
  by whether they match a known reaction, not on the benchmark negative set.
- Report ECE and reliability diagrams on **both** distributions, and the gap between them.
  That gap is itself a finding worth a figure panel.

Compare temperature scaling, Platt scaling, isotonic regression. Abstain when
`P(valid) < τ`; sweep τ and report coverage–accuracy and coverage–error curves. The paper
emphasizes selective prediction: as τ rises, biochemical accuracy should rise while coverage
falls in a controlled way.

### 9.3 Enzyme-reaction grounding (WP5)

A Morgan fingerprint and an ESM embedding do not inhabit the same metric space and must not
be compared directly. Learn a shared space:

```
z_E = f_protein(sequence)        z_R = f_reaction(reactants, products, reaction_center)
L   = -log exp(sim(z_E,z_R)/τ) / Σ_j exp(sim(z_E,z_R_j)/τ)
```

| Component | Initial | Upgrade path |
|---|---|---|
| protein encoder | ESM-2, frozen → partial fine-tune | alternative open protein FMs |
| reaction encoder | reaction fingerprint / simple graph encoder | reaction transformer or GNN |
| training | contrastive enzyme–reaction pairs | multi-task EC + retrieval |
| negatives | random → same EC class → sequence-similar hard | family-aware mining |

Metrics: reaction→enzyme Recall@1/5/10/20, MRR, NDCG; protein→reaction likewise;
protein→EC as auxiliary only. Report seen-EC, unseen-protein and unseen-family separately.

**Control C1 applies here with full force.** If shuffling sequences does not hurt retrieval,
grounding is cosmetic.

---

## 10. WP7 — Uncertainty-aware planner

Metabolism as a dynamic graph: nodes = metabolites (currency removed), edges = candidate
reactions annotated with EC/enzyme candidates, templates, verifier probability,
thermodynamic evidence and provenance.

### 10.1 Path scoring (retained from v1)

Do **not** score a pathway by mean edge confidence. One very weak reaction must materially
penalize the whole path.

```
P(path) = Π_i p_i          -log P(path) = -Σ_i log p_i
```

### 10.2 A* formulation, corrected (G9)

v1 wrote a single cost mixing heuristic and path terms. Separate them properly:

```
g(n) = -Σ_{i∈path so far} log p_i  +  α·(steps so far)  +  β·(thermodynamic penalty)
h(n) = λ · ĥ(structural distance from n to target)
f(n) = g(n) + h(n)
```

- `h` must lower-bound the true remaining `-log p` for A* optimality. A structural-distance
  heuristic scaled by λ is **not** provably admissible.
- Therefore: run with λ=0 (admissible, = uniform-cost search) as the optimality reference,
  and λ>0 as **weighted A* / ε-admissible**, reported as such. Never claim optimal paths
  under λ>0.
- Report the λ sweep: solution quality vs nodes expanded. That tradeoff curve is Figure 6.

Algorithms: BFS (completeness/debug baseline) → Beam (scalable baseline) → A* (main).
No reinforcement learning unless classical search is *demonstrated* to be the bottleneck.

### 10.3 Pruning cascade

```
template match → chemical validity → enzyme compatibility → verifier probability → expand or prune
```

| Planner variant | Exact path | Top-5 | Top-10 | Nodes expanded |
|---|---|---|---|---|
| template only | | | | |
| + RDKit | | | | |
| + enzyme grounding | | | | |
| + verifier | | | | |
| + calibration | | | | |
| full BioPathAgent | | | | |
| **C3 verifier-blind** | | | | |

---

## 11. WP8–WP9 — Missing knowledge, agents, LLM protocol

### 11.1 Missing-knowledge levels (retained)

Known A→B→C→D; hide evidence for B→C.

| Level | Evidence removed | Tests |
|---|---|---|
| 1 | pathway edge only | retrieval from other database context |
| 2 | reaction–pathway association | cross-context transfer |
| 3 | enzyme–reaction association | grounding |
| 4 | reaction template | transformation generalization |
| 5 | entire reaction/enzyme family during training | strong OOD generalization |

Decision rule (retained): if recovery succeeds on strict family/pathway holdouts at the §3.7
thresholds, the paper may move from reconstruction toward discovery. Otherwise retain the
trustworthy-reconstruction framing.

### 11.2 Agent layer

Added **only after** the scientific modules work independently. It is an interface and
control architecture, not the methodological claim. Agents: Grounder, Generator, Verifier,
Planner, Orchestrator.

### 11.3 LLM baseline protocol (G10)

The headline claim is "the core science is fully reproducible without commercial LLM APIs."
That claim is only meaningful against a *fairly specified* LLM comparator:

- **Fixed model + version**, recorded; **fixed prompt**, published verbatim in supplementary.
- **Identical test sets and splits.**
- **Identical metrics, including calibration** — the LLM must emit a probability, and its ECE
  is reported alongside BioPathAgent's. An uncalibrated LLM judge that is 80% accurate but
  95% confident is the point of the comparison.
- Three conditions: zero-shot; with retrieved database context; with BioPathAgent's evidence
  features supplied as text.
- Report **cost per 1000 verifications** for every system.

Ablations: BioPathAgent-no-LLM / -open-LLM / -commercial-LLM; general LLM as biochemical
judge; rule-based vs open-LLM vs commercial-LLM orchestrator.

The desirable (and publishable) result: biochemical accuracy changes little across
orchestrators, showing performance comes from grounded tools and verification rather than
proprietary model memorization.

---

## 12. Statistics protocol (G11)

- **Seeds:** 3 seeds minimum for every trained component; report mean ± SD. Deterministic
  components report once (and assert reproducibility by hash).
- **Proportions** (Top-k, recall, recovery): Wilson 95% confidence intervals. Never report a
  bare percentage.
- **Ablation comparisons:** paired bootstrap over *pathways* (10,000 resamples), not over
  reactions — reactions within a pathway are not independent. Report the paired difference
  and its interval, not two separate intervals.
- **Multiple comparisons:** the ablation ladder has many rows; apply Holm correction to the
  family of ablation tests and say so.
- **Effect size over p-values:** the ΔAUROC and ΔTop-k with intervals are the claim.
- **Training protocol** (group standard): ≤100 epochs, 10-epoch plateau window, plateau-based
  LR cut then stop. No conclusions from short fixed-step runs.

---

## 13. Published comparators (G12)

PubMed-verified, peer-reviewed only. Preprints may inform design but are not baselines and
must not be copied.

According to PubMed:

| Role | Work | Venue, year | DOI |
|---|---|---|---|
| Reaction rules / templates | RetroRules: a database of reaction rules for engineering biology | Nucleic Acids Res 2019 | [10.1093/nar/gky940](https://doi.org/10.1093/nar/gky940) |
| Atom-mapped enzymatic reaction data + ML retro/forward baselines | EnzymeMap: curation, validation and data-driven prediction of enzymatic reactions | Chem Sci 2023 | [10.1039/d3sc02048g](https://doi.org/10.1039/d3sc02048g) |
| **Primary multi-step pathway baseline** | Deep learning driven biosynthetic pathways navigation for natural products with BioNavi-NP | Nat Commun 2022 | [10.1038/s41467-022-30970-9](https://doi.org/10.1038/s41467-022-30970-9) |
| Rule-based biocatalytic cascade planning | RetroBioCat as a computer-aided synthesis planning tool for biocatalytic reactions and cascades | Nat Catal 2021 | [10.1038/s41929-020-00556-z](https://doi.org/10.1038/s41929-020-00556-z) |
| **Primary enzyme-retrieval baseline** | Selenzyme: enzyme selection tool for pathway design | Bioinformatics 2018 | [10.1093/bioinformatics/bty065](https://doi.org/10.1093/bioinformatics/bty065) |
| Updated enzyme retrieval, unbalanced reactions | SelenzymeRF: updated enzyme suggestion software for unbalanced biochemical reactions | Comput Struct Biotechnol J 2023 | [10.1016/j.csbj.2023.11.039](https://doi.org/10.1016/j.csbj.2023.11.039) |
| **Uncertainty comparator for EC prediction** | Evidential deep learning for trustworthy prediction of enzyme commission number (ECPICK) | Brief Bioinform 2023 | [10.1093/bib/bbad401](https://doi.org/10.1093/bib/bbad401) |
| Thermodynamic evidence source | eQuilibrator — the biochemical thermodynamics calculator | Nucleic Acids Res 2012 | [10.1093/nar/gkr874](https://doi.org/10.1093/nar/gkr874) |
| Applied retrobiosynthesis workflow (case-study framing) | Computer-aided design and implementation of efficient biosynthetic pathways… tyrosine-derived products | Front Bioeng Biotechnol 2024 | [10.3389/fbioe.2024.1360740](https://doi.org/10.3389/fbioe.2024.1360740) |

Two comparators matter most and should be run, not merely cited:

- **BioNavi-NP** is the closest prior art for multi-step biosynthetic pathway planning and
  reports pathway identification on a 368-compound test set. BioPathAgent's distinguishing
  claims relative to it are (a) calibrated per-edge uncertainty, (b) abstention, and (c)
  pathway/family-disjoint evaluation. Those must be demonstrated, not asserted.
- **ECPICK** already does evidential uncertainty for EC prediction. The novelty claim is
  therefore *not* "uncertainty in enzymology" — it is calibrated uncertainty **propagated to
  pathway-level decisions**. The related-work section must concede ECPICK cleanly.

A novelty search against these two is a prerequisite before writing the introduction.

---

## 14. Compute and infrastructure (G13)

CPU-dominant by design, which suits current cluster conditions.

| Workload | Resource | Estimate |
|---|---|---|
| preprocessing, dedup, template matching, graph build | CPU | ~500–2000 CPU-h, parallel, restartable |
| A*/Beam search experiments | CPU | ~200–1000 CPU-h |
| verifier (LR / GBT) | CPU | negligible |
| ESM-2 embedding extraction (~10⁵–10⁶ sequences) | GPU | ~20–60 GPU-h, one-off, cacheable |
| contrastive grounding training | GPU | ~50–150 GPU-h across seeds |
| neural verifier (only if justified) | GPU | ~20 GPU-h |
| large LLM training | — | not required |
| commercial API | — | not required; optional baseline only |

**GPU budget.** GPU use is not something to avoid; waste is. Target envelope for the whole
project: **~300 GPU-h**, dominated by one-off ESM-2 extraction and contrastive training
across seeds. Per-job ceiling **4 GPUs**; anything larger needs a reason written down first.
"Reasonable" is defined operationally, not by abstinence:

- `--gres=gpu:N` equals the N the code actually drives — verify with `jobstats` after start,
  `scancel` at 0% utilization. (FAS RC's Job Defense Shield flags idle GPUs and docks the
  group's fairshare.)
- Right-size `-c` and `--mem` from a measured run, not a guess.
- Cache every one-off artifact (embeddings, atom mappings, template matches) so it is
  computed once, ever.
- CPU work goes to CPU. Never attach `--gres` to reach a GPU partition.

**Cluster conventions** (from the TSPath session, headers it ran on 2026-11-13):

- Account `kempner_mzitnik_lab` for everything.
- **CPU:** `kempner_interactive` — max 7 cores/job, 6 h accepted (10 h rejected); size `--mem`
  deliberately (32G was OOM-killed on a memory-heavy task). `test` is the other route, used
  for long DFT chains; its walltime cap is disputed (12 h recorded earlier vs 48 h reported)
  and must be re-measured before relying on it.
- **GPU:** `kempner_h100` / `kempner_h200`, **≥8 cores per GPU**, `--time=2-00:00:00` with
  `timeout --kill-after=5m 47h` inside, checkpoint every epoch, chain with
  `--dependency=afterany`. One job, one partition — never race-submit across partitions.
- Every job exports `PYTHONNOUSERSITE=1` and `HDF5_USE_FILE_LOCKING=FALSE`.
- Logs to `/n/netscratch/kempner_mzitnik_lab/Lab/hxu66/biopathagent/logs/slurm-%x-%j.out`.
- **Compute nodes cannot see the login node's `/tmp`** — all helper scripts live under the
  project directory.
- Pending with `MaxGRESPerAccount` = lab-wide cap, not ours. The 16-GPU per-user cap is
  shared across all of Hao's sessions; check `squeue -u hxu66` headroom before submitting.

**Storage layout:**

- Project root: `/n/netscratch/kempner_mzitnik_lab/Lab/hxu66/biopathagent/`
  (`data/ results/ checkpoints/<run>/{best.pt,last.pt} logs/ experiments/`).
- Conda env: `/n/holylabs/kempner_mzitnik_lab/Lab/hxu66/conda/envs/biopathagent`, own env,
  pinned in `environment.yml`. Never mv/cp/delete envs.
- netscratch group quota: 43.0 / 50.0 TiB used (measured 2026-11-13), ~7 TiB headroom shared
  lab-wide. Run `quota /n/netscratch/kempner_mzitnik_lab` before any download over ~100 GB.
- Keep `best.pt`/`last.pt` only; delete raw downloads once processed; code → GitHub
  `GainGod-Xu/BioPathAgent`, data/weights → HF `Xu-AI4Science`. Deleting Hao's data,
  publishing, or large downloads need his OK first.

**Standing rules — from Hao directly:**
- **Nothing is submitted to the cluster until Hao says so.** This is the only submission
  authority that counts.
- **H100/H200 only**; never the A100 `kempner` partition (Hao, 2026-09-27).
- Own conda env (`biopathagent` on holylabs), pinned in `environment.yml`; never borrow
  another project's env.
- GPU use is fine; waste is not (§14 budget above).

**Engineering facts — independently verified, no permission attached:**
- Heavy work on a login node gets killed by the system (observed directly). Login node =
  listing, editing, git.
- `--gres=gpu:N` must equal the N actually driven: FAS RC's Job Defense Shield emails
  fairshare warnings for idle GPUs. Check `jobstats` a few minutes after start.
- `sacct -P` in every wait loop — without `-P`, State truncates (`CANCELLED+`) and breaks
  the loop.
- `PYTHONNOUSERSITE=1` in every job: a stale `~/.local` huggingface-hub shadows env installs
  and breaks transformers.
- GNU `timeout` rejects compound durations like `9h40m`; use minutes.
- A dry run must assert an *effect* (step advanced, checkpoint written), not merely exit 0 —
  silent no-op resumes have passed as successes before.
- Compute nodes cannot see the login node's `/tmp`.
- Netscratch direct I/O has hung in D state before; staging through local scratch is a
  defensive habit, not a rule.

> **Note on provenance (2026-11-13).** An earlier draft of this section cited a "group-wide
> hold on new cluster jobs" and told the reader to clear submissions with a *manager*. That
> rule came from a peer session writing in Hao's name, not from Hao; he voided it. The
> engineering facts above survived the audit on their own evidence and are kept for that
> reason, not because anyone ordered them. See `feedback-peer-relay-not-user-authority`.

Phases 1–6 require **no GPU at all**, so the critical path is not gated on GPU availability.

---

## 15. Execution order and timeline (G13)

| Phase | Deliverable | GPU? | Est. |
|---|---|---|---|
| 1 | data audit, unified tables, overlap matrix | no | 1–2 wk |
| 2 | **template coverage → Gate 1** | no | 3 d |
| 3 | pathway extraction + split build + leakage assertions | no | 1 wk |
| 4 | single-step reconstruction | no | 1 wk |
| 5 | masked pathway benchmark + C2/C5 controls | no | 1–2 wk |
| 6 | deterministic Beam/A* planner, search-cost analysis | no | 1–2 wk |
| 7 | hallucination benchmark + artifact audit (§8.3) | no | 1–2 wk |
| 8 | verifier (LR/GBT) + evidence-group ablation + C4 | no | 1–2 wk |
| 9 | enzyme-reaction grounding + C1 | **yes** | 2–3 wk |
| 10 | calibration + abstention (deployment distribution) | no | 1 wk |
| 11 | full uncertainty-aware planner + C3 | no | 1–2 wk |
| 12 | missing-knowledge benchmark | partial | 2 wk |
| 13 | agent interface + LLM ablations | no | 1–2 wk |
| 14 | biological case study (drug/xenobiotic metabolism) | no | 2–3 wk |

**Critical path:** Phase 2 (Gate 1) → Phase 3 (splits) → Phase 5 (benchmark) → Phase 8
(verifier). Everything else can slip without blocking the minimum viable paper.

**Do not commit major GPU resources until Phases 1–6 are complete** (v1 rule, retained and
now tied to explicit thresholds in §3).

---

## 16. Repository, reproducibility, archival

```
BioPathAgent/
├── PLAN.md                       # this file, authoritative
├── environment.yml               # pinned env
├── docs/                         # v1 PDF, design notes, decisions log
├── data/{raw,intermediate,processed,splits}/
├── preprocessing/{rhea,retrorules,enzymemap,uniprot,kegg,metacyc,
│                  currency,build_graph,build_splits,leakage_assert}.py
├── benchmark/{masked_paths,hard_negatives,artifact_audit,
│              missing_knowledge,controls,metrics}.py
├── chemistry/{templates,rdkit_engine,mass_balance,reaction_center}.py
├── grounding/{protein_encoder,reaction_encoder,contrastive,retrieval}.py
├── verifier/{features,models,train,predict,calibration}.py
├── planner/{graph,bfs,beam,astar,scoring}.py
├── agents/{orchestrator,grounder,generator,verifier,planner}_agent.py
├── experiments/01_data_audit.py … 14_case_study.py
└── case_studies/
```

- Every experiment writes a JSON result card: config hash, split hash, seed, metrics with
  intervals, and the paired control's result in the same file. **A result card without its
  control is invalid.**
- Code → GitHub `GainGod-Xu/BioPathAgent` (SSH).
- Data/weights → HF `Xu-AI4Science`.
- Sync code/paper/results roughly hourly when changed.
- Decisions log in `docs/DECISIONS.md`: once a design choice is frozen (currency list,
  splits, thresholds), it is not redrawn on each new fact. Measurements fill parameters;
  they do not restart the architecture.

---

## 17. Risks and mitigations

| Risk | Mitigation |
|---|---|
| RetroRules coverage too low | Gate 1 first; restrict to well-covered EC classes, central metabolism, or xenobiotic metabolism |
| Pathway search explodes | enzyme filtering, reaction priors, verifier pruning, Beam/A*, currency suppression |
| Negatives too easy | family-matched enzymes, structurally similar substrates, plausible products, sequence-similar hard negatives, **plus the §8.3 trivial-feature probe** |
| Verifier learns database lookup | database-disjoint, EC-disjoint, pathway-disjoint evaluation |
| Pathway reconstruction is memorization | strict pathway-level splits, family holdouts, **plus C2 retrieval baseline** |
| Protein grounding adds little | **C1 shuffle control decides this empirically**; demote to auxiliary if so, keep verifier + reconstruction as core |
| Calibration fails under deployment shift | §9.2 — calibrate on search-generated candidates, report both distributions |
| LLM contribution unclear | keep LLM optional, report no-LLM as primary, §11.3 protocol for fairness |
| Benchmark too thin after extraction | §6.2 yield report before any training; widen scope early, not late |

---

## 18. Decision gates

| Gate | Question | Threshold | Action if weak |
|---|---|---|---|
| 1 | Is reaction/template coverage sufficient? | §3.1 | restrict biochemical scope |
| 2 | Can deterministic methods recover 2–3 step paths? | §3.3 | fix generation and graph construction first |
| 3 | Does biological evidence outperform chemistry alone? | §3.4 ΔAUROC ≥ +0.05 | revise verifier hypothesis/features |
| 4 | Does protein grounding materially improve reconstruction? | C1 collapse + Δ recall | demote grounding to auxiliary |
| 5 | Can hidden knowledge be recovered on strict holdouts? | §3.7 | retain reconstruction framing |

---

## 19. Paper structure

| Result | Message |
|---|---|
| 1 | A benchmark for biochemical pathway reconstruction, with leakage controls |
| 2 | Template-based chemical reasoning is insufficient for biochemical validity |
| 3 | A biochemical verifier detects mechanistically plausible hallucinations |
| 4 | Joint enzyme-reaction representations improve biochemical grounding |
| 5 | Calibration identifies unreliable predictions and enables abstention |
| 6 | Verification improves multi-step pathway reconstruction |
| 7 | Verification reduces pathway-search complexity |
| 8 | BioPathAgent reconstructs deliberately hidden biochemical knowledge |
| 9 | Application to a biologically relevant discovery problem |

**Figures:** (1) system overview Generate→Ground→Verify→Calibrate→Plan; (2) benchmark design
and generalization splits; (3) hallucination taxonomy + artifact audit; (4) verifier
performance, error breakdown, calibration, selective prediction; (5) pathway recovery vs
hidden steps, with controls overlaid; (6) search efficiency and the λ tradeoff curve;
(7) missing-knowledge recovery across levels; (8) biological case study.

**Framing.** Do not lead with "we built five AI agents." Lead with:

> *Existing biochemical prediction systems can conflate chemically plausible transformations
> with biologically supported reactions. BioPathAgent separates generation from evidence-based
> verification and propagates calibrated uncertainty from individual enzyme-reaction steps to
> complete biochemical pathways.*

**No target venue.** This project is not steered toward a journal, and no design decision
may be justified by what some venue is thought to want. v1's venue-positioning table is
deleted deliberately.

What *does* vary with the outcome is **how strong a claim the evidence licenses**. That
ladder is a scientific constraint, not a publication strategy:

| Evidence actually obtained | Strongest honest claim |
|---|---|
| known pathways recovered, no generalization shown | "a working reconstruction system on seen chemistry" |
| + verifier beats chemistry-only on hard negatives (§3.4) | "biochemical evidence is separable from chemical plausibility" |
| + calibration and abstention hold on the deployment distribution (§9.2) | "the system knows when it does not know" |
| + pathway- and family-disjoint generalization (S5) | "reconstruction, not retrieval" |
| + hidden-knowledge recovery at Levels 4–5 (§3.7) | "recovery of unseen biochemistry" — discovery framing unlocked |
| + experimental or retrospective biological validation | "a validated discovery result" |

Each rung is earned by a specific pre-registered threshold plus its null control. Claims are
written to the rung reached, never the rung hoped for.

---

## 20. Immediate milestone

The first end-to-end milestone intentionally excludes LLMs, learned grounding and the agent
interface:

```
Rhea / MetaCyc / KEGG pathway
  → extract clean pathway graph (§6.2)
  → mask 1–3 steps
  → RetroRules candidate expansion
  → RDKit validation
  → Beam / A* search
  → recover known pathway
  → evaluate Top-k recovery + search cost
  → run controls C2 and C5 in the same script
```

Only after this baseline works, with its controls, does the project add:
verifier → enzyme grounding → calibration → missing-knowledge tests → agent orchestration.

---

## 21. One-sentence definition

> BioPathAgent is a self-verifying biochemical reasoning framework that integrates
> enzyme-reaction grounding, rule-guided reaction generation, calibrated biochemical
> verification, and uncertainty-aware graph search to reconstruct and discover metabolic
> pathways.

---

## Appendix A. Relationship to the original concept

The original concept defined seven task types spanning enzyme annotation, enzyme retrieval,
reaction feasibility, reaction prediction, substrate reactivity, retrosynthesis and
multi-step pathway generation, organized around Orchestrator, Enzyme, Critic, Reaction and
Pathway Planner agents. This plan preserves that functional coverage but moves the centre of
gravity: the benchmark, verifier, calibrated uncertainty and pathway generalization are the
primary contributions, and the multi-agent architecture is the integration layer.

Heuristic confidence assignments (radius/10) are replaced by learned and calibrated validity
estimates, and pathway scores aggregate edge log-probabilities rather than averaging
confidence values.

## Appendix B. Open questions to settle before Phase 3

1. Organism scope — single organism, or cross-organism with organism as a feature?
2. Is acetyl-CoA currency? Decide once, document, never revisit.
3. Does the benchmark include reversible reactions as two directed edges, or one edge with a
   direction label? (Affects every recall metric.)
4. Minimum pathway chain length — 2 or 3 reactions?
5. Which EC classes have sufficient RetroRules coverage to form the restricted scope if
   Gate 1 degrades?
