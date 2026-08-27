# CAPYBARA technical audit before redesign

Date: 2026-08-27  
Audited upstream snapshot: `main` at `4ff09dada3d42700fff6ba92609d20955e5b8578`  
Scope: current code and database, README, examples, and selected historical revisions. No production classifier code was modified during this audit.

## Executive conclusion

**OBSERVED:** The current program is not a hierarchical phylogenetic classifier. It independently scores hard-coded lineage and terminal-marker dictionaries, discards lineage evidence whenever any terminal marker passes, and uses `DP * AF` as the score. The lineage and terminal panels do not check the observed REF/ALT allele at all.

**OBSERVED:** A wrong allele at the correct coordinate can support a label. Multiallelic AF is calculated from the first ALT and `AD[1]`, regardless of the expected ALT. A single two-read terminal call can suppress a 50-unit incompatible parent call.

**OBSERVED:** The external `variant.SNP` contains the exact conflicting barcode `123384 C>T` for both `1.3.2` and `1.4.2`. The active hard-coded dictionary silently retains only `1.3.2`, so the text database and runtime database disagree.

**OBSERVED:** The lineage records for `2.2` (`592890 A>G`) and `2.3` (`266845 G>T`) disagree with the current ESL reference, whose bases are respectively `G` and `T`. Their declared ALT is therefore the reference allele. Because the caller emits variant sites only, these parent-lineage records cannot work as declared in the present coordinate/reference system.

**IMPLICATION:** Fixing only REF/ALT comparison is necessary but insufficient. The database, hierarchy, evidence model, input-specific calling, broad ESL/GC gate, abstention policy, and output semantics all require redesign and validation against the supplied strain data.

## Repository and provenance

The current GitHub repository could not be cloned through the terminal Git transport because `github.com:443` timed out. The current `main` archive was retrieved successfully through GitHub codeload and verified locally with SHA-256:

`0886F2ED103CFDA0C69F133BC48CB1E60BEDAC0859106B404435F381BF673573`

The GitHub API reported 79 commits, with the current head dated 2026-03-25. Historical source snapshots used for this audit were fetched by immutable commit SHA:

- `c9ba7b19aaaf12005a0e0e7bbf30cc69de4aaae7` (2025-09 metagenomic revision)
- `119c20b18e8faf6a3d7d1f779b1fb8f62e6ce023` (2024-07 Mash-based revision)
- `b9c742918bbfeaaeb38a7c43d8e0fc9ccf459c27` (2024-05 Mash-based revision)

Current repository inventory:

- `capy.py`, `configure.py`, `README.md`, `LICENSE`, and workflow image
- 3 assembly examples: `1.3.9.fa`, `2.5.6.fa`, and `Non-esl.fa`
- ESL reference: one chromosome, 3,964,912 bases
- HC1030 reference: chromosome plus plasmid; marker evaluation targets the first FASTA contig
- 22 HC1030 markers, 8 lineage rows, and 51 terminal rows at 50 unique terminal positions
- 5,823 Mash sketch files, but only 5,818 entries in `cc.mash.list`; five sketches are not indexed by that list

No user-supplied strain metadata or per-strain mutation dataset was present in `D:\Codex_tools` at audit time. Files belonging to unrelated projects were not assumed to be CAPYBARA inputs.

## Actual current execution path

1. `configure.py` resolves executable and database paths, reads `cc.mash.list`, and then defines all HC1030, lineage, and terminal biology as Python dictionaries.
2. `capy.py` imports those dictionaries. It does not load `hc1030.snp`, `lineage.SNP`, or `variant.SNP`.
3. The query is mapped to the ESL reference with the same minimap2 command for assemblies and reads; no `asm5` or `sr` preset is selected.
4. Only panel coordinates are sent to `bcftools mpileup` and `bcftools call -mv`.
5. For lineage and terminal panels, the program decides support from depth and AF alone.
6. Terminal scores are computed first. If any terminal score survives relative filtering, lineage scores are not considered.
7. In isolate mode, a second mapping to the HC1030 reference counts exact HC marker matches. Any nonzero votes for exactly one of GC1 or GC2 makes the HC gate "trusted"; the winning GC is not checked against the lineage/terminal label.
8. A compact summary, a marker table, and JSON are written. Lineage/terminal marker status is always `.` and observed alleles are not reported.

## Requested problem audit

### Problem A — expected REF/ALT is not checked

**OBSERVED, confirmed.** In `evaluate_panel()`, `ref_obs` and `alt_obs` are parsed but never used in `pass_flag`. The returned tuple does not retain either observed allele. `score_variants()` therefore cannot distinguish `T>G` from `T>A` at position 1,733,723.

A direct synthetic invocation of the current scoring function gave a score of `19.0` to a `2.5.6` record with `DP=20`, `AF=0.95`, and `ok=True`; the tuple has no field in which the wrong observed ALT could be represented.

### Problem B — ALT-specific AF

**OBSERVED, confirmed.** `parse_vcf_with_af()` selects `alt.split(",")[0]`, reads only `AD[1]`, and uses `AD[0] + AD[1]` as the denominator. It neither locates the expected ALT index nor includes other alleles in the depth denominator.

For `ALT=A,G; AD=1,3,20`, the current code reports AF `3/(1+3)=0.75` for the first ALT. It cannot calculate the expected-G AF `20/(1+3+20)=0.8333`.

### Problem C — terminal clade overrides lineage evidence

**OBSERVED, confirmed.** `filtered_lineages` is computed only when `filtered_variants` is empty. A synthetic strong `2.3` score of 50 and weak `2.5.6` score of 2 retained `2.5.6` and discarded all lineage evidence.

The lineage printed for a terminal label is reconstructed by joining the first two dot-separated components. It is not an independently supported call.

### Problem D — single SNP per terminal clade

**OBSERVED, confirmed.** Every active terminal node has exactly one coordinate. Every active lineage node also has exactly one coordinate. Only the broad HC panels contain multiple markers.

### Problem E — duplicated/non-unique barcode

**OBSERVED, confirmed.** `variant.SNP` assigns `123384 C>T` to both `1.3.2` and `1.4.2`. This is emitted in `analysis/marker_conflicts.tsv` by the reproducible audit script.

### Problem F — scheme drift

**OBSERVED, confirmed.** The current runtime uses hard-coded dictionaries in `configure.py`; the external SNP files are not parsed. The external terminal table has 51 rows, while the active dictionary has 50 positions and silently omits the `1.4.2` assignment at the duplicated position.

**OBSERVED, historical regression.** The 2024-07 configuration loaded the external SNP files into dictionaries. The 2026 configuration changed back to hard-coded biological dictionaries while retaining external files and their paths. The README still describes behavior that is absent from current code.

### Problem G — incomplete hierarchy

**OBSERVED, confirmed.** Fifteen `2.4.x` terminal nodes exist, but no `2.4` parent marker exists in the lineage panel. The current program can still emit a `2.4.x` terminal call and fabricate lineage `2.4` from its string prefix.

No explicit hierarchy database exists. Parentage is inferred only for output formatting and, in the 2025 version, by string-prefix filtering.

### Problem H — classification must allow abstention

**OBSERVED, confirmed.** Current output has no status vocabulary for `PASS`, `AMBIGUOUS`, `UNRESOLVED`, `INSUFFICIENT_COVERAGE`, `CONFLICTING_MARKERS`, or `NON_ESL`. It can emit multiple summary rows for nodes within 20% of the maximum score, but does not explain this as ambiguity.

The HC gate clears all fine-scale calls when both GC groups have at least one vote or neither group has a vote. It does not distinguish insufficient coverage, non-ESL biology, contamination, mixed samples, or marker conflict.

## Additional observed defects

### Reference/database inconsistency

**OBSERVED:** Two of eight lineage REF alleles disagree with `esl_ref.fna`:

| Node | Position | Declared change | ESL reference base | Consequence |
|---|---:|---|---|---|
| 2.2 | 592890 | A>G | G | expected ALT is reference allele |
| 2.3 | 266845 | G>T | T | expected ALT is reference allele |

All 22 HC records and all 51 terminal rows have REF alleles matching their respective current reference at the declared coordinate.

### GC/path inconsistency is not tested

**OBSERVED:** An exclusive GC1 vote makes `hc_trusted=True`, but a GC2 lineage or terminal result remains eligible. The code does not enforce `GC1 -> 1.x` or `GC2 -> 2.x`.

### Isolate depth threshold is ignored

**OBSERVED:** In `evaluate_panel()`, isolate mode uses `dp >= 1` even when the user supplies `--min-snp-depth`. The default is also one. This is inconsistent with the intended isolate-read default of approximately ten and makes one-read errors eligible.

### Assembly and read evidence are conflated

**OBSERVED:** One `-i/--input` argument accepts everything. The same mapper command, DP logic, AF cutoff, and variant caller are applied to assemblies and reads. Paired `-1/-2` input is unsupported. Assembly alignment multiplicity can be interpreted as read depth.

### Mapping/calling quality is not controlled

**OBSERVED:** No explicit minimap2 preset, MAPQ filter, base-quality filter, BAQ policy, duplicate policy, strand requirement, or ambiguous-mapping exclusion is present. The `threads` argument of `targeted_call_with_af()` is unused.

### DP dominates phylogenetic evidence

**OBSERVED:** `DP * AF` gives unlimited weight to a deeply covered single locus. It does not count independent markers, callable markers, missingness, contradictions, sensitivity, specificity, or phylogenetic consistency.

### Marker output is not an observation table

**OBSERVED:** The marker report prints expected REF and ALT under generic `ref` and `alt` columns, discards observed alleles and ALT depth, and prints `.` for every lineage/terminal status. This can make a nonmatching locus appear like the expected mutation was observed.

### Broad ESL logic differs from documentation

**OBSERVED:** The README says the current workflow screens against 5,824 pre-sketched genomes. Current `capy.py` never calls Mash. It uses only the HC1030 panel in isolate mode, and skips that gate entirely in metagenomic mode.

**OBSERVED:** The 2024 Mash implementation selected the nearest sketch and accepted its clonal-complex label without an absolute distance threshold or nearest-neighbour margin. The initial comparison value was merely a programming sentinel, not a calibrated ESL threshold.

### Metagenomic mode is not mixture-aware

**OBSERVED:** Metagenomic mode skips HC gating yet still applies terminal-first fine-scale scoring. The default depth is two and AF is 0.8. It has no multi-marker AF coherence rule and no mixed-lineage status.

### Operational robustness

**OBSERVED:** External commands are assembled as unquoted strings and executed with `shell=True`, so paths containing spaces can fail and untrusted path text can be interpreted by the shell. Executable discovery may return `None`, which is then interpolated into a command rather than rejected at startup.

**OBSERVED:** No database schema version, startup validation, marker provenance, threshold provenance, or standalone validation command exists.

## Historical behavior

**OBSERVED:** The 2024-07 implementation performed Mash nearest-neighbour CC assignment, then scanned any called variant whose coordinate appeared in the terminal dictionary. It checked position but not expected REF/ALT and had no hierarchy.

**OBSERVED:** The 2025-09 implementation introduced a partial parent-first rule: it chose one lineage, then considered terminal labels with that string prefix. However, both lineage and terminal decisions still used only DP and AF, not the expected allele. Each node still had one SNP.

**OBSERVED:** The 2026-03 implementation regressed from that partial parent gate: terminal evidence is evaluated first, and any passing terminal score prevents lineage scoring.

## What cannot yet be concluded

**HYPOTHESIS:** Some current SNPs are likely homoplastic, recurrent, recombinant, or cross-lineage beyond the known `123384 C>T` conflict. This cannot be established from the shipped database alone.

**HYPOTHESIS:** Multi-SNP panels will reduce false terminal assignments with only a modest rise in unresolved calls. The magnitude cannot be reported until metadata and callable per-strain alleles are supplied and leakage-aware validation is run.

**NOT YET MEASURABLE:** sensitivity, specificity, missingness, sister-node occurrence, number of robust markers per node, optimal panel size, old-versus-new performance, false terminal-clade rate, and performance cost of abstention.

## Reproducible audit artifacts

- `scripts/audit_existing_db.py`: strict audit of external tables, active Python dictionaries, reference alleles, parent presence, and duplicate alleles
- `analysis/current_marker_audit.tsv`: one row per shipped external marker
- `analysis/current_scheme_summary.tsv`: panel-level counts and discrepancies
- `analysis/marker_conflicts.tsv`: exact duplicate allele/node conflicts

The audit script intentionally exits nonzero when conflicts are present, so it can become a CI/startup validation building block.

## Implications and next action

1. Do not patch the existing terminal-first score in place.
2. Obtain the strain metadata and mutation/callability data and inspect their actual formats.
3. Resolve reference identity and coordinate conventions before marker discovery, especially the `2.2` and `2.3` reversed/reference-allele records.
4. Build a metadata-backed hierarchy and explicitly review exceptions rather than relying only on dot splitting.
5. Quantify all current markers on held-out data before selecting new panels.
6. Benchmark the original, an exact-allele bug-fixed single-SNP baseline, and hierarchical multi-SNP candidates before changing production classification.
7. Only after threshold selection, replace the duplicated dictionaries with one structured, validated database and implement mode-specific evidence calling.
