# Decision record

Design changes made after a measurement contradicted something the plan had
frozen. One entry per change. A frozen choice is not unchangeable - it is
changeable only with a number attached, and this file is where the number goes.

Entries are append-only. If a decision is later reversed, add a new entry that
says so; do not edit the old one.

---

## D1 - the currency-metabolite rule (2026-11-13)

**What the plan said** (PLAN.md 6.3, v1 and Revision 2):

> A participant is currency if it is in the curated list OR its reaction-degree is
> above the 99.5th percentile. Removed as graph nodes, retained as verifier
> features.

**What the measurement showed** (experiment 03, first run):

The degree half of the rule is correct at the top and wrong below it. Ranks 1-21
by degree are exactly the cofactors one would curate by hand: H+, H2O, O2, CoA,
ATP, NADP+, NADPH, NAD+, NADH, diphosphate, phosphate, CO2, SAM, ADP, FMN, SAH,
AMP, NH4+, H2O2. Further down, the same rule also selects:

| formula | metabolite | degree | breadth |
|---|---|---|---|
| `C5H4O5-2` | 2-oxoglutarate | 474 | 58 |
| `C23H34N7O17P3S-4` | acetyl-CoA | 415 | 55 |
| `C4H4O4-2` | fumarate | 374 | 41 |
| `C5H8NO4-` | L-glutamate | 252 | 57 |
| `C6H12O6` | glucose | 239 | 23 |
| `C3H3O3-` | pyruvate | 221 | 56 |
| `C2H3O2-` | acetate | 191 | 38 |
| `C3H7NO2` | L-alanine | 156 | 31 |
| `C18H33O2-` | oleate | 118 | 2 |

Removing those as graph nodes deletes glycolysis and the TCA cycle. The rule as
written would not have failed loudly; it would have produced pathway results that
looked fine and meant nothing. 53 metabolites selected, 9 of 9 labelled
intermediates wrongly among them.

**The hypothesis tested instead** (experiment 04):

Degree conflates "used throughout metabolism" with "used heavily inside a few
pathways". So try pathway breadth - the number of distinct ordinary KEGG maps a
metabolite participates in - on the grounds that a cofactor is broad and a central
intermediate is deep but narrow.

Result: **breadth does not separate the two classes either.** Reference cofactors
run down to breadth 21 (FAD); reference intermediates run up to 58
(2-oxoglutarate). The hypothesis as stated is false.

**What was adopted, and why it is not the same as tuning a parameter:**

```
currency(m) = breadth(m) >= 59  OR  m in config.CURATED_COFACTORS
```

Three things make this defensible rather than a fitted threshold:

1. **The two errors are not symmetric, so the threshold is set for zero false
   positives, not best accuracy.** A false positive deletes a real metabolite from
   the graph and corrupts every downstream pathway result. A false negative leaves
   a cofactor in the graph as an extra hub node, which costs search time and
   nothing else. Accuracy would trade these against each other as if they were
   comparable. They are not.

2. **59 sits in a gap, not on an edge.** Nothing in the data has breadth between
   58 and 75, so every threshold in `[59, 75]` yields the identical 13-member
   automatic set with zero errors (experiment 03 section E). The threshold is
   detecting a real bimodality in the distribution rather than being fitted to the
   references. Below 59 the error count rises immediately: at 55 it is 4 false
   positives, at 40 it is 5.

3. **The curated list carries the residual, which is the job a curated list is
   for.** 18 cofactors sit below the threshold (AMP, SAM, SAH, FMN, FMNH2, FAD,
   FADH2, NH4+, H2O2, UDP, GTP, GDP, Fe/Fe2S2 species, Na+, Cl-). No automatic
   criterion reaches them without also reaching pyruvate. Listing them is honest;
   pretending a threshold found them would not be.

On the 27 hand-labelled references: breadth alone gives 0 false positives and 5
false negatives; breadth plus the curated list gives 0 and 0.

**Known limitation, stated rather than buried:** the 27 references were labelled by
hand by the author of the rule, and the threshold is not cross-validated. The
claim is only "zero errors on the labelled set, in the direction that matters",
plus the gap argument in (2), which does not depend on the labels. A larger
independent labelling would strengthen this. The full output and the sensitivity
sweep are in `results/03_currency_audit.json` so the rule can be re-judged without
re-deriving it.

**Implementation changes:**

- `config.CURRENCY_NAMES` (25 names) -> `config.CURATED_COFACTORS` (31 entries
  keyed by InChIKey-14, each carrying its molecular formula and justification).
  Names could not be joined to anything: Rhea participants are SMILES and we do
  not download ChEBI names, so the curated half of the original rule was
  unenforceable as written.
- The formula beside each key is **asserted at load time**
  (`chem.verify_curated_formulae`). A curated list keyed by opaque identifiers with
  the molecule named only in a comment is a list that rots silently.
- `config.CURRENCY_DEGREE_PERCENTILE` is retained and experiment 03 still prints
  the superseded rule's output, so this revision stays auditable against what it
  replaced.

**Side effect - PLAN.md appendix B question 1 is now closed.** The plan
deliberately left acetyl-CoA, L-glutamate and 2-oxoglutarate undecided. All three
fall below the threshold and all three carry the transferred carbon skeleton
rather than a transferable group, so all three are **nodes, not currency**. Also
settled the same way: malonyl-CoA (carbon donor) and formate (C1 metabolite).
Recorded in `config.SETTLED_NOT_CURRENCY`, which experiment 03 checks against the
frozen set on every run - if a future change pulls one of them in, the run says so.

---
