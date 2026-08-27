"""Interpretable hierarchical classification over evaluated canonical markers."""

from __future__ import annotations

from dataclasses import dataclass, asdict

from .database import CladeRule, Database, Marker


MATCH = "MATCH"
CALLABLE_STATUSES = frozenset({"MATCH", "REFERENCE", "WRONG_ALT", "LOW_AF", "CONFLICT"})
CONFLICT_STATUSES = frozenset({"WRONG_ALT", "CONFLICT"})


@dataclass(frozen=True)
class MarkerCall:
    position: int
    observed_ref: str
    observed_alt: str
    depth: int
    expected_alt_depth: int
    allele_fraction: float
    mapq: float | None
    status: str


@dataclass
class NodeEvidence:
    node: str
    marker_n: int
    callable_n: int
    matched_n: int
    conflicting_n: int
    support: float
    passed: bool
    rule_type: str
    matched_markers: list[str]
    conflicting_markers: list[str]


@dataclass
class Classification:
    sample: str
    mode: str
    esl: str
    gc: str
    lineage: str
    sublineage: str
    status: str
    confidence: str
    lineage_basis: str
    gc_support: str
    lineage_support: str
    sublineage_support: str
    callable_markers: int
    matched_markers: int
    conflicting_markers: int
    reads_total: int
    reads_mapped: int
    alternatives: str
    gate_method: str = ""
    gate_result: str = ""
    gate_distance: str = ""
    gate_margin: str = ""

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def marker_id(marker: Marker) -> str:
    return f"{marker.position}:{marker.ref}>{marker.alt}"


def marker_key(marker: Marker) -> tuple[int, str, str]:
    return marker.position, marker.ref, marker.alt


def minimum_evidence(marker_n: int) -> tuple[int, float]:
    if marker_n <= 1:
        return 1, 1.0
    if marker_n == 2:
        return 2, 1.0
    return 2, 0.60


def evidence_for_node(
    node: str,
    markers: list[Marker],
    calls: dict[tuple[int, str, str], MarkerCall],
    rule_type: str,
) -> NodeEvidence:
    callable_n = matched_n = conflicting_n = 0
    matched_markers: list[str] = []
    conflicting_markers: list[str] = []
    for marker in markers:
        call = calls.get(marker_key(marker))
        if call is None:
            continue
        if call.status in CALLABLE_STATUSES:
            callable_n += 1
        if call.status == MATCH:
            matched_n += 1
            matched_markers.append(marker_id(marker))
        elif call.status in CONFLICT_STATUSES:
            conflicting_n += 1
            conflicting_markers.append(marker_id(marker))
    support = matched_n / callable_n if callable_n else 0.0
    required_matches, required_support = minimum_evidence(len(markers))
    passed = (
        callable_n >= required_matches
        and matched_n >= required_matches
        and support >= required_support
        and conflicting_n == 0
    )
    return NodeEvidence(
        node=node, marker_n=len(markers), callable_n=callable_n,
        matched_n=matched_n, conflicting_n=conflicting_n, support=support,
        passed=passed, rule_type=rule_type, matched_markers=matched_markers,
        conflicting_markers=conflicting_markers,
    )


def _marker_lookup(database: Database) -> dict[tuple[str, str], Marker]:
    return {(marker.node, marker_id(marker)): marker for marker in database.markers}


def evaluate_nodes(
    database: Database,
    calls: dict[tuple[int, str, str], MarkerCall],
) -> dict[str, NodeEvidence]:
    by_node = database.markers_by_node
    evidence: dict[str, NodeEvidence] = {}
    for node, markers in by_node.items():
        canonical = [marker for marker in markers if marker.marker_set == "CANONICAL"]
        if canonical:
            rule_type = database.rules[node].rule_type if node in database.rules else "CANONICAL_PANEL"
            evidence[node] = evidence_for_node(node, canonical, calls, rule_type)

    lookup = _marker_lookup(database)
    for node, rule in database.rules.items():
        if rule.rule_type != "OVERRIDE_THEN_FALLBACK":
            continue
        anchor = lookup[(node, rule.anchor_marker)]
        anchor_call = calls.get(marker_key(anchor))
        override_matches: list[str] = []
        for override_node, marker_text in rule.override_markers:
            override = lookup[(override_node, marker_text)]
            call = calls.get(marker_key(override))
            if call is not None and call.status == MATCH:
                override_matches.append(f"{override_node}={marker_text}")
        anchor_match = anchor_call is not None and anchor_call.status == MATCH
        evidence[node] = NodeEvidence(
            node=node,
            marker_n=1,
            callable_n=int(anchor_call is not None and anchor_call.status in CALLABLE_STATUSES),
            matched_n=int(anchor_match),
            conflicting_n=len(override_matches),
            support=1.0 if anchor_match else 0.0,
            passed=bool(anchor_match and not override_matches),
            rule_type=rule.rule_type,
            matched_markers=[rule.anchor_marker] if anchor_match else [],
            conflicting_markers=override_matches,
        )
    return evidence


def _support_text(item: NodeEvidence | None) -> str:
    if item is None:
        return "0/0"
    return f"{item.matched_n}/{item.callable_n}"


def _confidence(item: NodeEvidence, rule: CladeRule | None) -> str:
    if item.matched_n >= 3 and item.support >= 0.80 and item.conflicting_n == 0:
        confidence = "HIGH"
    else:
        confidence = "MEDIUM"
    if rule and rule.confidence_ceiling == "MEDIUM":
        confidence = "MEDIUM"
    return confidence


def classify(
    database: Database,
    calls: dict[tuple[int, str, str], MarkerCall],
    sample: str,
    mode: str,
    reads_total: int = 0,
    reads_mapped: int = 0,
) -> tuple[Classification, dict[str, NodeEvidence]]:
    evidence = evaluate_nodes(database, calls)
    passed_lineages = sorted(
        node for node, item in evidence.items()
        if item.passed and database.hierarchy[node].level == "lineage"
    )
    passed_clades = sorted(
        node for node, item in evidence.items()
        if item.passed and database.hierarchy[node].level == "sublineage"
    )
    callable_total = sum(call.status in CALLABLE_STATUSES for call in calls.values())
    matched_total = sum(call.status == MATCH for call in calls.values())
    conflict_total = sum(call.status in CONFLICT_STATUSES for call in calls.values())

    if len(passed_clades) == 1:
        clade = passed_clades[0]
        lineage = database.hierarchy[clade].parent
        gc = database.hierarchy[lineage].parent
        item = evidence[clade]
        lineage_item = evidence.get(lineage)
        incompatible_lineages = [node for node in passed_lineages if node != lineage]
        if incompatible_lineages:
            result = Classification(
                sample=sample, mode=mode, esl="ESL", gc="NA", lineage="NA",
                sublineage="NA", status="CONFLICTING_MARKERS", confidence="LOW",
                lineage_basis="INCOMPATIBLE_PARENT_AND_CHILD",
                gc_support="0/0", lineage_support="0/0",
                sublineage_support=_support_text(item),
                callable_markers=callable_total, matched_markers=matched_total,
                conflicting_markers=conflict_total, reads_total=reads_total,
                reads_mapped=reads_mapped,
                alternatives=";".join(incompatible_lineages + [clade]),
            )
            return result, evidence
        result = Classification(
            sample=sample, mode=mode, esl="ESL", gc=gc, lineage=lineage,
            sublineage=clade, status="PASS", confidence=_confidence(item, database.rules.get(clade)),
            lineage_basis="INFERRED_FROM_CHILD_CLADE",
            gc_support=_support_text(evidence.get(gc)),
            lineage_support=_support_text(lineage_item),
            sublineage_support=_support_text(item),
            callable_markers=callable_total, matched_markers=matched_total,
            conflicting_markers=conflict_total, reads_total=reads_total,
            reads_mapped=reads_mapped, alternatives="",
        )
        return result, evidence

    if len(passed_clades) > 1:
        parents = {database.hierarchy[node].parent for node in passed_clades}
        lineage = next(iter(parents)) if len(parents) == 1 else "NA"
        gc_values = {
            database.hierarchy[parent].parent for parent in parents
            if parent in database.hierarchy
        }
        gc = next(iter(gc_values)) if len(gc_values) == 1 else "NA"
        status = "MIXED" if mode == "metagenomic" else "AMBIGUOUS"
        result = Classification(
            sample=sample, mode=mode, esl="ESL", gc=gc, lineage=lineage,
            sublineage="NA", status=status, confidence="LOW",
            lineage_basis="CONFLICTING_CHILD_CLADES",
            gc_support="0/0", lineage_support="0/0", sublineage_support="0/0",
            callable_markers=callable_total, matched_markers=matched_total,
            conflicting_markers=conflict_total, reads_total=reads_total,
            reads_mapped=reads_mapped, alternatives=";".join(passed_clades),
        )
        return result, evidence

    if len(passed_lineages) == 1:
        lineage = passed_lineages[0]
        gc = database.hierarchy[lineage].parent
        item = evidence[lineage]
        result = Classification(
            sample=sample, mode=mode, esl="ESL", gc=gc, lineage=lineage,
            sublineage="NA", status="UNRESOLVED", confidence=_confidence(item, None),
            lineage_basis="DIRECT_MARKERS", gc_support=_support_text(evidence.get(gc)),
            lineage_support=_support_text(item), sublineage_support="0/0",
            callable_markers=callable_total, matched_markers=matched_total,
            conflicting_markers=conflict_total, reads_total=reads_total,
            reads_mapped=reads_mapped, alternatives="",
        )
        return result, evidence

    mapped_fraction = reads_mapped / reads_total if reads_total else 0.0
    non_esl = callable_total == 0 and mapped_fraction < 0.01
    insufficient = callable_total == 0 and not non_esl
    result = Classification(
        sample=sample, mode=mode, esl="NON_ESL" if non_esl else "ESL",
        gc="NA", lineage="NA", sublineage="NA",
        status="NON_ESL" if non_esl else ("INSUFFICIENT_COVERAGE" if insufficient else "UNRESOLVED"),
        confidence="LOW", lineage_basis="NONE", gc_support="0/0",
        lineage_support="0/0", sublineage_support="0/0",
        callable_markers=callable_total, matched_markers=matched_total,
        conflicting_markers=conflict_total, reads_total=reads_total,
        reads_mapped=reads_mapped,
        alternatives=";".join(passed_lineages) if passed_lineages else "",
    )
    return result, evidence
