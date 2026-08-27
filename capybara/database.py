"""Load and validate the external CAPYBARA barcode database."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


DNA = frozenset("ACGT")


@dataclass(frozen=True)
class HierarchyNode:
    node: str
    parent: str
    level: str
    classifiable: bool


@dataclass(frozen=True)
class Marker:
    node: str
    parent: str
    level: str
    position: int
    ref: str
    alt: str
    homoplasy: int
    sensitivity: float
    specificity: float
    weight: float
    marker_set: str
    status: str


@dataclass(frozen=True)
class CladeRule:
    node: str
    parent: str
    rule_type: str
    anchor_marker: str
    conflicting_clades: tuple[str, ...]
    override_markers: tuple[tuple[str, str], ...]
    confidence_ceiling: str
    status: str


@dataclass
class Database:
    root: Path
    hierarchy: dict[str, HierarchyNode]
    markers: list[Marker]
    rules: dict[str, CladeRule]
    reference: Path

    @property
    def markers_by_node(self) -> dict[str, list[Marker]]:
        result: dict[str, list[Marker]] = {}
        for marker in self.markers:
            result.setdefault(marker.node, []).append(marker)
        return result

    @property
    def markers_by_position(self) -> dict[int, list[Marker]]:
        result: dict[int, list[Marker]] = {}
        for marker in self.markers:
            result.setdefault(marker.position, []).append(marker)
        return result


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise ValueError(f"Required database file is missing: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _parse_marker_text(value: str) -> tuple[int, str, str]:
    try:
        position_text, change = value.split(":", 1)
        ref, alt = change.split(">", 1)
        return int(position_text), ref, alt
    except Exception as exc:
        raise ValueError(f"Invalid marker expression: {value!r}") from exc


def _read_reference(path: Path) -> str:
    sequence = "".join(
        line.strip() for line in path.open(encoding="utf-8")
        if line and not line.startswith(">")
    ).upper()
    if not sequence:
        raise ValueError(f"Reference FASTA is empty: {path}")
    return sequence


def is_ancestor(ancestor: str, descendant: str, hierarchy: dict[str, HierarchyNode]) -> bool:
    current = descendant
    seen: set[str] = set()
    while current in hierarchy and current not in seen:
        if current == ancestor:
            return True
        seen.add(current)
        current = hierarchy[current].parent
    return False


def load_database(db_dir: Path, validate: bool = True) -> Database:
    hierarchy_rows = _read_tsv(db_dir / "hierarchy.tsv")
    marker_rows = _read_tsv(db_dir / "markers.tsv")
    rule_rows = _read_tsv(db_dir / "clade_rules.tsv")
    hierarchy = {
        row["node"]: HierarchyNode(
            node=row["node"],
            parent=row["parent"],
            level=row["level"],
            classifiable=row["classifiable"].lower() == "true",
        )
        for row in hierarchy_rows
    }
    markers = [
        Marker(
            node=row["node"], parent=row["parent"], level=row["level"],
            position=int(row["position"]), ref=row["ref"].upper(), alt=row["alt"].upper(),
            homoplasy=int(row["homoplasy"]), sensitivity=float(row["sensitivity"]),
            specificity=float(row["specificity"]), weight=float(row["weight"]),
            marker_set=row["marker_set"], status=row["status"],
        )
        for row in marker_rows
    ]
    rules: dict[str, CladeRule] = {}
    for row in rule_rows:
        overrides: list[tuple[str, str]] = []
        for item in row["override_markers"].split(";") if row["override_markers"] else []:
            label, marker = item.split("=", 1)
            overrides.append((label, marker))
        rules[row["node"]] = CladeRule(
            node=row["node"], parent=row["parent"], rule_type=row["rule_type"],
            anchor_marker=row["anchor_marker"],
            conflicting_clades=tuple(filter(None, row["conflicting_clades"].split(";"))),
            override_markers=tuple(overrides),
            confidence_ceiling=row["confidence_ceiling"], status=row["status"],
        )
    database = Database(
        root=db_dir, hierarchy=hierarchy, markers=markers, rules=rules,
        reference=db_dir / "esl_ref.fna",
    )
    if validate:
        errors, warnings = validate_database(database)
        if errors:
            raise ValueError("Invalid CAPYBARA database:\n- " + "\n- ".join(errors))
        database.validation_warnings = warnings  # type: ignore[attr-defined]
    return database


def validate_database(database: Database) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    hierarchy = database.hierarchy
    if "ESL" not in hierarchy:
        errors.append("hierarchy lacks ESL root")
    for node, item in hierarchy.items():
        if item.parent != "ROOT" and item.parent not in hierarchy:
            errors.append(f"unknown parent {item.parent!r} for {node!r}")
        seen: set[str] = set()
        current = node
        while current in hierarchy:
            if current in seen:
                errors.append(f"hierarchy cycle involving {node!r}")
                break
            seen.add(current)
            current = hierarchy[current].parent

    if not database.reference.exists():
        errors.append(f"reference FASTA missing: {database.reference}")
        reference = ""
    else:
        reference = _read_reference(database.reference)

    exact_rows: set[tuple[str, int, str, str, str]] = set()
    allele_nodes: dict[tuple[int, str], set[str]] = {}
    for marker in database.markers:
        if marker.node not in hierarchy:
            errors.append(f"marker refers to unknown node {marker.node!r}")
            continue
        if marker.parent != hierarchy[marker.node].parent or marker.level != hierarchy[marker.node].level:
            errors.append(f"marker hierarchy fields disagree for {marker.node!r}")
        if marker.ref not in DNA or marker.alt not in DNA or marker.ref == marker.alt:
            errors.append(f"invalid REF/ALT at {marker.node}:{marker.position} {marker.ref}>{marker.alt}")
        if marker.homoplasy > 10:
            errors.append(f"active marker exceeds homoplasy 10: {marker.node}:{marker.position}")
        if marker.position < 1 or (reference and marker.position > len(reference)):
            errors.append(f"marker position outside reference: {marker.node}:{marker.position}")
        elif reference and reference[marker.position - 1] != marker.ref:
            errors.append(
                f"reference mismatch at {marker.node}:{marker.position}; "
                f"database={marker.ref}, FASTA={reference[marker.position - 1]}"
            )
        key = (marker.node, marker.position, marker.ref, marker.alt, marker.marker_set)
        if key in exact_rows:
            errors.append(f"duplicate marker row: {key}")
        exact_rows.add(key)
        allele_nodes.setdefault((marker.position, marker.alt), set()).add(marker.node)

    for allele, nodes in allele_nodes.items():
        unrelated = [
            (left, right) for left in nodes for right in nodes
            if left < right
            and not is_ancestor(left, right, hierarchy)
            and not is_ancestor(right, left, hierarchy)
        ]
        if unrelated:
            warnings.append(f"allele {allele[0]}->{allele[1]} shared by unrelated nodes: {unrelated}")

    marker_keys = {(m.node, m.position, m.ref, m.alt) for m in database.markers}
    for node, rule in database.rules.items():
        if node not in hierarchy:
            errors.append(f"rule refers to unknown node {node!r}")
            continue
        if rule.parent != hierarchy[node].parent:
            errors.append(f"rule parent disagrees for {node!r}")
        if rule.rule_type == "OVERRIDE_THEN_FALLBACK":
            position, ref, alt = _parse_marker_text(rule.anchor_marker)
            if (node, position, ref, alt) not in marker_keys:
                errors.append(f"fallback anchor absent from markers.tsv: {node}={rule.anchor_marker}")
            for override_node, marker_text in rule.override_markers:
                position, ref, alt = _parse_marker_text(marker_text)
                if (override_node, position, ref, alt) not in marker_keys:
                    errors.append(f"override marker absent from markers.tsv: {override_node}={marker_text}")
        elif rule.rule_type == "NO_RESOLVABLE_RULE":
            warnings.append(f"node lacks resolvable marker rule: {node}")
        elif rule.rule_type != "EXCLUSIVE_MARKER":
            errors.append(f"unknown rule type for {node}: {rule.rule_type}")

    for node, item in hierarchy.items():
        if item.classifiable and item.level == "sublineage" and node not in database.rules:
            errors.append(f"classifiable sublineage lacks rule: {node}")
    return errors, warnings
