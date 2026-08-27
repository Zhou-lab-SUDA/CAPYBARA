# CAPYBARA 2 (development)

CAPYBARA is a lightweight hierarchical canonical-SNP classifier for epidemic
super-lineages of *Acinetobacter baumannii*. Version 2 replaces independent
single-position score lookup with exact-allele, hierarchy-aware barcode rules.

This development database was reconstructed from `r4.labelled.nwk`, strain
metadata, and the corresponding branch-mutation history. It must be benchmarked
on independent genomes before release as a stable clinical or surveillance tool.

## Key behavior

- The expected `REF>ALT` allele is checked exactly.
- Multiallelic AD is indexed using the expected ALT, not fixed `AD[1]`.
- Assembly mapping uses minimap2 `asm5`; short reads use `sr`.
- A child clade cannot override strong evidence for an incompatible parent.
- Sibling support produces `AMBIGUOUS` or `MIXED`, not an arbitrary winner.
- A clade can use one canonical SNP when that is the available phylogenetic
  information; confidence remains limited.
- Ordered `OVERRIDE_THEN_FALLBACK` rules represent shared anchors. For example,
  `2860886 T>C` identifies 2.5.8 before `1733723 T>G` falls back to 2.5.6.
- Missing calls never count as positive marker matches.
- Marker definitions have one source of truth in external TSV files.

## Requirements

- Python 3.9 or newer
- minimap2
- samtools
- bcftools
- Mash (assembly ESL/GC gate)

The executables must be available on `PATH`.

## Database validation

```bash
python capy.py --validate-db
```

Validation checks hierarchy integrity, REF/ALT validity, reference coordinates,
homoplasy limits, duplicate rows, rule targets, missing rule markers, and
unresolved nodes. Errors prevent classification; scientifically explicit
limitations are warnings.

## Assembly mode

```bash
python capy.py --assembly genome.fna.gz -o results/sample -t 4
```

Backward-compatible single-assembly input remains available:

```bash
python capy.py -i genome.fna.gz -o results/sample
```

Assembly defaults are DP >= 1, AF >= 0.90, MAPQ >= 20, and base quality >= 20.
Alignment depth represents callable alignment support, not sequencing coverage.
Before barcode traversal, Mash compares the assembly to the labelled reference
sketch database. The development defaults require distance <= 0.02 and a
nearest-class margin >= 0.0005. A nearest non-ST1/ST2 reference returns
`NON_ESL`; Mash does not define the terminal clade.

## Isolate paired reads

```bash
python capy.py -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -o results/sample -t 4
```

Read defaults are DP >= 10, expected-ALT AF >= 0.90, MAPQ >= 20, and base
quality >= 20.

## Metagenomic mode

```bash
python capy.py --metagenomic -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -o results/sample -t 4
```

Single-file metagenomic input is also accepted with `--metagenomic -i`. Initial
development defaults are DP >= 2 and expected-ALT AF >= 0.20. These thresholds
are provisional and must be calibrated with the validation data. Simultaneously
supported clades are reported as `MIXED` rather than collapsed to one label.

## Outputs

For prefix `results/sample`, CAPYBARA writes:

- `results/sample.summary.tsv`: compact hierarchical result and status.
- `results/sample.markers.tsv`: every evaluated marker, expected and observed
  alleles, DP, expected-ALT depth, AF, quality filter, and marker status.
- `results/sample.evidence.json`: thresholds and node-level evidence.

Possible statuses include `PASS`, `UNRESOLVED`, `AMBIGUOUS`, `MIXED`,
`INSUFFICIENT_COVERAGE`, `CONFLICTING_MARKERS`, and `NON_ESL`.

## Authoritative database

The runtime database is:

```text
capydb/esl/
├── esl_ref.fna
├── hierarchy.tsv
├── markers.tsv
└── clade_rules.tsv
```

`configure.py` contains only paths, executable discovery, and default settings.
It contains no biological marker dictionaries.

## Rebuilding the barcode

The reproducible analytical scripts are kept separate from `capy.py`:

```text
scripts/parse_metadata.py
scripts/parse_mutations.py
scripts/build_hierarchy.py
scripts/discover_markers.py
scripts/build_clade_rules.py
scripts/select_markers.py
```

Candidate selection uses homoplasy <= 10. Exclusive candidates are checked
against the complete reconstructed mutation distribution, including events from
other clades that would themselves be rejected for high homoplasy. Shared
anchors can only be used through explicit ordered override/fallback rules.

## Tests

```bash
python -m unittest discover -s tests -v
```

The tests include wrong ALT, multiallelic ALT depth, weak child versus
incompatible parent, incomplete child panels, sibling conflict, missing loci,
non-ESL input, mixed samples, and the 2.5.8-over-2.5.6 rule.

See `docs/technical_audit.md` and `docs/marker_selection_policy.md` for the
current audit evidence and marker policy.
