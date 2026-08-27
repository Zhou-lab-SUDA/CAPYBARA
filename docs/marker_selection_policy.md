# CAPYBARA marker selection policy

This document records the hard eligibility rules applied before panel ranking.
Ranking must never restore a marker rejected by these rules.

## Hard eligibility rules

A node-specific canonical SNP is selectable only when all of the following hold:

1. Its SNPPar homoplasy value is no greater than 10.
2. For an exclusive marker, its reconstructed expected ALT occurs in no sample
   outside the target node.
3. Its sensitivity meets the configured discovery threshold.
4. Its specificity meets the configured discovery threshold.
5. Its REF and ALT are valid, distinct nucleotides and the event is assigned to
   the target MRCA branch.

The cross-node conflict check is performed against the complete mutation history
and reconstructed tip states. Events belonging to another lineage are included
in this check even when those events are themselves ineligible because their
homoplasy exceeds 10. Filtering high-homoplasy candidates before conflict
detection is prohibited.

Occurrences in descendants are not conflicts for a parent-node marker because
descendants are members of the parent node. Any occurrence in a sister or other
unrelated node is a conflict for an exclusive marker.

Ancestral reconstruction supplies inferred tip alleles, not original per-sample
callability. Therefore discovery output reports `missing_rate=NA` with
`callability_basis=ANCESTRAL_RECONSTRUCTION`; missingness must be measured from
the original variant/callability inputs when those data are available.

## Consequence for 2.5.6

The historical `1733723 T>G` marker has homoplasy 3 and is reconstructed in
3239/3241 labelled 2.5.6 samples. It is also reconstructed in all 78 labelled
2.5.8 samples. It is therefore rejected as an exclusive 2.5.6 marker with
`CROSS_NODE_CONFLICT`, despite its high sensitivity.

Under the hard rules above, the current tree/mutation dataset contains no
eligible *exclusive* marker for 2.5.6. The historical SNP is instead retained as
a conditional anchor. Because 2.5.8 has an eligible exclusive marker at
`2860886 T>C`, the database rule is:

```text
2860886 T>C matched
=> 2.5.8

otherwise, 1733723 T>G matched
=> 2.5.6
```

This is an ordered override/fallback rule. If the 2.5.8 marker is not returned,
including when its locus is uncallable, the 2.5.6 anchor remains the fallback.
The output must retain the evidence details so downstream users can distinguish
an explicit reference observation from a missing guard call.

This conditional-exclusion construction is generated from the database for all
clades. It must not be hard-coded as a special Python branch for 2.5.6.

## Panel selection

Among eligible markers, selection should prefer lower homoplasy, higher
sensitivity, adequate genomic separation, and independent genomic regions.
Three to five markers per node remain a preference, not a requirement. A node
with one eligible SNP may use a validated single-marker panel with explicitly
limited confidence. A node with no eligible SNP remains unresolved at that
level.
