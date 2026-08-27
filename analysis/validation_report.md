# CAPYBARA 2 development validation report

## Scope and interpretation

This report records the pre-server development validation. The 10,528 tree-tip
replay uses the same labelled population collection from which markers were
discovered. It is an apparent-performance and regression benchmark, not an
independent estimate of generalization. Independent assembly and metagenomic
testing is required next.

## OBSERVED: database and implementation validation

- The production database contains 188 active markers for 53 nodes.
- Thirty-three represented nodes have at least three selected markers; twelve
  represented nodes currently have one selected marker.
- The terminal rule database contains 45 exclusive rules, two ordered
  override/fallback rules, and five deliberately unresolved nodes.
- Database validation reports zero errors.
- The five warnings are `1.4.9`, `2.2.5`, `2.2.6`, `2.4.2`, and `2.4.7`, for
  which the current automatic rule builder did not find a resolved rule under
  the frozen criteria.
- Ten unit/stress tests pass, including exact wrong-ALT rejection,
  multiallelic expected-ALT AD selection, incompatible parent/child evidence,
  incomplete child panels, sibling conflict, missing loci, non-ESL input,
  mixed samples, and the 2.5.8-over-2.5.6 decision rule.

## OBSERVED: tree-tip replay

| Classifier | Terminal accuracy | Macro F1 | Resolved | Unresolved | False terminal assignment |
|---|---:|---:|---:|---:|---:|
| Original CAPYBARA behavior | 88.70% | 70.91% | 98.18% | 1.82% | 10.35% |
| Exact-ALT bug-fixed single-SNP | 88.81% | 70.98% | 98.18% | 1.82% | 10.24% |
| Hierarchical multi-SNP/rule classifier | 84.94% | 90.38% | 84.07% | 15.93% | 0.00% |

The exact values and per-sample predictions are in `benchmark_summary.tsv`,
`benchmark_old.tsv`, `benchmark_bugfixed.tsv`, and `benchmark_new.tsv`.
Confusion matrices are stored under `confusion_matrices/`.

## INFERRED

- Correcting REF/ALT handling alone removes only a small fraction of erroneous
  terminal calls in this replay. Barcode design and non-hierarchical terminal
  override behavior are the dominant failure modes.
- The frozen hierarchical rules trade classification rate for a large reduction
  in false fine-scale assignments, matching the stated scientific objective.
- Macro F1 rises despite lower raw terminal accuracy because the new classifier
  avoids systematically assigning large clades to incorrect legacy labels.
- The zero false-terminal result is encouraging but is optimistic because the
  same tree informed discovery. It must not be presented as independent
  performance.

## HYPOTHESIS to test on server data

- Independent assemblies should retain high specificity but may show additional
  unresolved calls from assembly gaps, mapping ambiguity, or population
  diversity absent from the discovery tree.
- Metagenomic samples should yield more `MIXED`, `AMBIGUOUS`, and
  `INSUFFICIENT_COVERAGE` outcomes than isolates; the initial DP/AF defaults may
  require recalibration.
- Some of the five unresolved clades may be recoverable with explicitly ordered
  positive-only fallback rules without increasing cross-clade false positives.

## Next action

The development version was deployed to `/home/naclist/Sc/Capybara_debug` and
all remote output was restricted to its `test/` directory.

## OBSERVED: independent server validation

- Remote database validation completed with zero errors and the same five known
  unresolved-node warnings.
- All eleven current local tests pass; the first ten also passed on the deployed
  server snapshot before the input-count regression test was added.
- Seventeen unique external NCBI assemblies were tested. Sixteen returned PASS
  and one returned UNRESOLVED. No sample produced a Mash/barcode GC conflict.
- The sixteen resolved assemblies comprised eleven 2.5.6, two 2.5.2, one
  2.4.10, and two 2.4.13 calls. These samples are outside the supplied labelled
  discovery metadata, so the result establishes technical operation and path
  consistency but not independent clade accuracy.
- The repository smoke assembly for 2.5.6 returned `GC2 > 2.5 > 2.5.6`.
- The repository non-ESL assembly had an ST4 nearest Mash reference and returned
  NON_ESL. The GC1 smoke assembly passed the Mash GC1 gate but remained
  unresolved because its current node lacks a frozen rule.
- Three small metagenomes were rerun after true FASTQ record counting was
  implemented. All returned NON_ESL with zero matched barcode markers:
  ERR1332582 had 62,073 reads/1 primary mapped, ERR1332586 had 6,406,763/485,
  and SRR15204137 had 958,451/976. Earlier pre-fix `INSUFFICIENT_COVERAGE`
  outputs remain archived for regression provenance but are not final results.
- The full positive metagenome SRR16308649 was UNRESOLVED under the new mapping
  and barcode rules. The old result had assigned legacy label 1.3.8 using
  `3649 C>A` at reported DP 12 and AF 1.0. In the reconstructed hierarchy this
  allele belongs to 1.4.8. The new `minimap2 -x sr --secondary=no` pipeline found
  the locus uncallable even when MAPQ and base-quality thresholds were set to
  zero, indicating a mapping-strategy difference rather than a demonstrated
  base-quality effect.

The collected remote artifacts and combined table are stored under
`analysis/server_validation/` and `analysis/server_validation_summary.tsv`.

## INFERRED from server validation

- The ordered 2.5.8/2.5.6 rule works end-to-end on assemblies: samples carrying
  the fallback anchor without the override marker consistently return 2.5.6.
- The Mash gate supplies the missing NON_ESL behavior and independently checks
  GC-level path consistency without defining the terminal clade.
- Conservative unresolved output prevented reproduction of at least one legacy
  fine-scale assignment whose label and mapping evidence are inconsistent with
  the new tree-aware scheme.

## Remaining validation limitations

- Independent truth labels were unavailable for the 17 newly sampled NCBI
  assemblies; their calls must not be counted as accuracy estimates.
- Only one previously positive metagenome was run at full depth. Metagenomic AF
  and multi-strain mixture thresholds remain provisional.
- Geographic, temporal, and clone-grouped holdouts have not yet been completed.
