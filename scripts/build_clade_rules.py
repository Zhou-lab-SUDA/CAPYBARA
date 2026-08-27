#!/usr/bin/env python3
"""Build database-driven exclusive and conditional rules for terminal clades."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_counts(value: str) -> dict[str, int]:
    result: dict[str, int] = {}
    for item in value.split(";") if value else []:
        name, count = item.rsplit(":", 1)
        result[name] = int(count)
    return result


def marker_text(row: dict[str, str]) -> str:
    return f"{row['position']}:{row['ref']}>{row['alt']}"


def rank_key(row: dict[str, str]) -> tuple[float, float, int, int]:
    return (
        -float(row["sensitivity"]),
        -float(row["specificity"]),
        int(row["homoplasy"]),
        int(row["position"]),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--hierarchy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-homoplasy", type=int, default=10)
    parser.add_argument("--min-anchor-sensitivity", type=float, default=0.98)
    args = parser.parse_args()

    candidates = read_tsv(args.candidates)
    hierarchy = read_tsv(args.hierarchy)
    clades = {
        row["node"]: row
        for row in hierarchy
        if row["level"] == "sublineage" and row["classifiable"] == "true"
    }
    by_node: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        if row["node"] in clades:
            by_node[row["node"]].append(row)

    exclusive: dict[str, list[dict[str, str]]] = {}
    for node in clades:
        rows = [row for row in by_node[node] if row["candidate_status"] == "PASS"]
        exclusive[node] = sorted(rows, key=rank_key)

    output: list[dict[str, object]] = []
    for node, hierarchy_row in sorted(clades.items()):
        if exclusive[node]:
            anchor = exclusive[node][0]
            output.append(
                {
                    "node": node,
                    "parent": hierarchy_row["parent"],
                    "rule_type": "EXCLUSIVE_MARKER",
                    "anchor_marker": marker_text(anchor),
                    "anchor_homoplasy": anchor["homoplasy"],
                    "anchor_sensitivity": anchor["sensitivity"],
                    "anchor_specificity": anchor["specificity"],
                    "conflicting_clades": "",
                    "override_markers": "",
                    "required_guard_state": "",
                    "uncallable_guard_result": "",
                    "confidence_ceiling": "MEDIUM" if len(exclusive[node]) == 1 else "HIGH",
                    "status": "READY",
                    "notes": f"eligible_exclusive_markers={len(exclusive[node])}",
                }
            )
            continue

        conditional: list[tuple[dict[str, str], dict[str, int]]] = []
        for row in by_node[node]:
            conflicts = parse_counts(row.get("outside_sublineages", ""))
            if (
                int(row["homoplasy"]) <= args.max_homoplasy
                and float(row["sensitivity"]) >= args.min_anchor_sensitivity
                and conflicts
                and all(conflict in clades and exclusive[conflict] for conflict in conflicts)
            ):
                conditional.append((row, conflicts))

        if conditional:
            anchor, conflicts = sorted(conditional, key=lambda item: rank_key(item[0]))[0]
            guards = [exclusive[conflict][0] for conflict in sorted(conflicts)]
            output.append(
                {
                    "node": node,
                    "parent": hierarchy_row["parent"],
                    "rule_type": "OVERRIDE_THEN_FALLBACK",
                    "anchor_marker": marker_text(anchor),
                    "anchor_homoplasy": anchor["homoplasy"],
                    "anchor_sensitivity": anchor["sensitivity"],
                    "anchor_specificity": anchor["specificity"],
                    "conflicting_clades": ";".join(sorted(conflicts)),
                    "override_markers": ";".join(
                        f"{conflict}={marker_text(guard)}"
                        for conflict, guard in zip(sorted(conflicts), guards)
                    ),
                    "required_guard_state": "NO_OVERRIDE_MARKER_RETURNED",
                    "uncallable_guard_result": "ASSIGN_FALLBACK",
                    "confidence_ceiling": "MEDIUM",
                    "status": "READY_CONDITIONAL",
                    "notes": (
                        "evaluate override markers first; when none is returned, "
                        "assign the anchor clade as fallback"
                    ),
                }
            )
        else:
            output.append(
                {
                    "node": node,
                    "parent": hierarchy_row["parent"],
                    "rule_type": "NO_RESOLVABLE_RULE",
                    "anchor_marker": "",
                    "anchor_homoplasy": "",
                    "anchor_sensitivity": "",
                    "anchor_specificity": "",
                    "conflicting_clades": "",
                    "override_markers": "",
                    "required_guard_state": "",
                    "uncallable_guard_result": "UNRESOLVED_AT_PARENT",
                    "confidence_ceiling": "",
                    "status": "UNRESOLVED",
                    "notes": "no exclusive marker and no conditionally resolvable anchor",
                }
            )

    fields = [
        "node", "parent", "rule_type", "anchor_marker", "anchor_homoplasy",
        "anchor_sensitivity", "anchor_specificity", "conflicting_clades",
        "override_markers", "required_guard_state", "uncallable_guard_result",
        "confidence_ceiling", "status", "notes",
    ]
    write_tsv(args.output, output, fields)
    counts: defaultdict[str, int] = defaultdict(int)
    for row in output:
        counts[str(row["rule_type"])] += 1
    print(" ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
