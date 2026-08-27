#!/usr/bin/env python3
"""Select compact, spatially separated CAPYBARA marker panels."""

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


def parse_marker(value: str) -> tuple[int, str, str]:
    position, change = value.split(":", 1)
    ref, alt = change.split(">", 1)
    return int(position), ref, alt


def rank_key(row: dict[str, str]) -> tuple[int, float, float, int]:
    return (
        int(row["homoplasy"]),
        -float(row["specificity"]),
        -float(row["sensitivity"]),
        int(row["position"]),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--rules", type=Path, required=True)
    parser.add_argument("--hierarchy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--max-per-node", type=int, default=5)
    parser.add_argument("--min-distance", type=int, default=50_000)
    args = parser.parse_args()

    reference = "".join(
        line.strip() for line in args.reference.open(encoding="utf-8")
        if not line.startswith(">")
    ).upper()

    candidates = read_tsv(args.candidates)
    rules = read_tsv(args.rules)
    hierarchy = {row["node"]: row for row in read_tsv(args.hierarchy)}
    by_node: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    lookup: dict[tuple[str, int, str, str], dict[str, str]] = {}
    for row in candidates:
        by_node[row["node"]].append(row)
        lookup[(row["node"], int(row["position"]), row["ref"], row["alt"])] = row

    selected: list[dict[str, object]] = []
    selected_keys: set[tuple[str, int, str, str, str]] = set()

    def add(row: dict[str, str], marker_set: str, status: str) -> None:
        position = int(row["position"])
        reference_base = reference[position - 1]
        if reference_base != row["ref"] or reference_base == row["alt"]:
            return
        key = (row["node"], position, row["ref"], row["alt"], marker_set)
        if key in selected_keys:
            return
        selected_keys.add(key)
        hrow = hierarchy[row["node"]]
        selected.append(
            {
                "node": row["node"],
                "parent": hrow["parent"],
                "level": hrow["level"],
                "position": row["position"],
                "ref": row["ref"],
                "alt": row["alt"],
                "homoplasy": row["homoplasy"],
                "sensitivity": row["sensitivity"],
                "specificity": row["specificity"],
                "weight": "1.0",
                "marker_set": marker_set,
                "status": status,
            }
        )

    # Exclusive panels for every discoverable node, not only terminal clades.
    for node, rows in sorted(by_node.items()):
        eligible = sorted(
            (row for row in rows if row["candidate_status"] == "PASS"),
            key=rank_key,
        )
        chosen: list[dict[str, str]] = []
        for row in eligible:
            position = int(row["position"])
            if all(abs(position - int(other["position"])) >= args.min_distance for other in chosen):
                chosen.append(row)
            if len(chosen) == args.max_per_node:
                break
        for row in chosen:
            add(row, "CANONICAL", "ACTIVE")

    # Preserve conditional fallback anchors and their override markers even when
    # the anchor is non-exclusive by itself.
    for rule in rules:
        if rule["rule_type"] != "OVERRIDE_THEN_FALLBACK":
            continue
        position, ref, alt = parse_marker(rule["anchor_marker"])
        anchor = lookup[(rule["node"], position, ref, alt)]
        add(anchor, "FALLBACK_ANCHOR", "ACTIVE_CONDITIONAL")
        for item in rule["override_markers"].split(";"):
            override_node, marker = item.split("=", 1)
            o_position, o_ref, o_alt = parse_marker(marker)
            override = lookup[(override_node, o_position, o_ref, o_alt)]
            add(override, "OVERRIDE", "ACTIVE")

    selected.sort(key=lambda row: (row["level"], row["node"], int(row["position"])))
    fields = [
        "node", "parent", "level", "position", "ref", "alt", "homoplasy",
        "sensitivity", "specificity", "weight", "marker_set", "status",
    ]
    write_tsv(args.output, selected, fields)
    node_counts: defaultdict[str, int] = defaultdict(int)
    for row in selected:
        node_counts[str(row["node"])] += 1
    print(f"markers={len(selected)} nodes={len(node_counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
