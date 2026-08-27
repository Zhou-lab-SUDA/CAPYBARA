#!/usr/bin/env python3
"""Build and validate the CAPYBARA hierarchy from normalized metadata."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def numeric_label(label: str) -> tuple[int, ...]:
    return tuple(int(part) for part in label.split("."))


def load_metadata(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_scheme_nodes(repo: Path) -> tuple[set[str], set[str]]:
    lineage_nodes: set[str] = set()
    sublineage_nodes: set[str] = set()
    for relative, target in [
        ("capydb/esl/lineage.SNP", lineage_nodes),
        ("capydb/esl/variant.SNP", sublineage_nodes),
    ]:
        with (repo / relative).open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    target.add(line.rstrip("\r\n").split("\t")[1])
    return lineage_nodes, sublineage_nodes


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--hierarchy", type=Path, required=True)
    parser.add_argument("--node-statistics", type=Path, required=True)
    parser.add_argument("--issues", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()

    rows = load_metadata(args.metadata)
    metadata_lineages = sorted({row["lineage"] for row in rows if row["lineage"]}, key=numeric_label)
    reported_sublineages = sorted(
        {row["reported_sublineage"] for row in rows if row["reported_sublineage"]},
        key=numeric_label,
    )
    scheme_lineages, scheme_sublineages = load_scheme_nodes(args.repo.resolve())

    hierarchy: list[dict[str, object]] = [
        {"node": "ESL", "parent": "ROOT", "level": "ESL", "classifiable": "true", "source": "definition"},
        {"node": "GC1", "parent": "ESL", "level": "GC", "classifiable": "true", "source": "metadata"},
        {"node": "GC2", "parent": "ESL", "level": "GC", "classifiable": "true", "source": "metadata"},
    ]
    for lineage in metadata_lineages:
        hierarchy.append(
            {
                "node": lineage,
                "parent": "GC1" if lineage.startswith("1.") else "GC2",
                "level": "lineage",
                "classifiable": "true",
                "source": "metadata",
            }
        )
    for sublineage in reported_sublineages:
        unresolved_bucket = sublineage.endswith(".0")
        hierarchy.append(
            {
                "node": sublineage,
                "parent": sublineage.rsplit(".", 1)[0],
                "level": "unresolved_bucket" if unresolved_bucket else "sublineage",
                "classifiable": "false" if unresolved_bucket else "true",
                "source": "metadata",
            }
        )
    write_tsv(
        args.hierarchy,
        hierarchy,
        ["node", "parent", "level", "classifiable", "source"],
    )

    members: defaultdict[str, set[str]] = defaultdict(set)
    countries: defaultdict[str, set[str]] = defaultdict(set)
    years: defaultdict[str, set[str]] = defaultdict(set)
    pas_sts: defaultdict[str, set[str]] = defaultdict(set)
    for row in rows:
        nodes: list[str] = []
        if row["esl_status"] == "ESL":
            nodes.append("ESL")
        if row["gc"]:
            nodes.append(row["gc"])
        if row["lineage"]:
            nodes.append(row["lineage"])
        if row["reported_sublineage"]:
            nodes.append(row["reported_sublineage"])
        for node in nodes:
            members[node].add(row["sample_id"])
            if row["country"]:
                countries[node].add(row["country"])
            if row["year"]:
                years[node].add(row["year"])
            if row["pas_st"]:
                pas_sts[node].add(row["pas_st"])
    node_statistics = [
        {
            "node": item["node"],
            "parent": item["parent"],
            "level": item["level"],
            "classifiable": item["classifiable"],
            "member_n": len(members[item["node"]]),
            "country_n": len(countries[item["node"]]),
            "year_value_n": len(years[item["node"]]),
            "pas_st_n": len(pas_sts[item["node"]]),
        }
        for item in hierarchy
        if item["node"] != "ROOT"
    ]
    write_tsv(
        args.node_statistics,
        node_statistics,
        ["node", "parent", "level", "classifiable", "member_n", "country_n", "year_value_n", "pas_st_n"],
    )

    issues: list[dict[str, object]] = []
    metadata_classifiable_sublineages = {
        item["node"] for item in hierarchy if item["level"] == "sublineage"
    }
    for node in sorted(set(metadata_lineages) - scheme_lineages, key=numeric_label):
        issues.append({"severity": "HIGH", "code": "METADATA_LINEAGE_MISSING_FROM_SCHEME", "node": node, "detail": "No current lineage marker"})
    for node in sorted(scheme_lineages - set(metadata_lineages), key=numeric_label):
        issues.append({"severity": "HIGH", "code": "SCHEME_LINEAGE_ABSENT_FROM_METADATA", "node": node, "detail": "Current scheme node not observed"})
    for node in sorted(metadata_classifiable_sublineages - scheme_sublineages, key=numeric_label):
        issues.append({"severity": "HIGH", "code": "METADATA_SUBLINEAGE_MISSING_FROM_SCHEME", "node": node, "detail": "No current terminal marker"})
    for node in sorted(scheme_sublineages - metadata_classifiable_sublineages, key=numeric_label):
        issues.append({"severity": "HIGH", "code": "SCHEME_SUBLINEAGE_ABSENT_FROM_METADATA", "node": node, "detail": "Current scheme node not observed in supplied ground truth"})
    for item in hierarchy:
        if item["level"] == "unresolved_bucket":
            issues.append({"severity": "INFO", "code": "UNRESOLVED_BUCKET", "node": item["node"], "detail": "Preserved from metadata but excluded as a classifiable child"})
    write_tsv(args.issues, issues, ["severity", "code", "node", "detail"])
    counts = Counter(item["code"] for item in issues)
    print("\n".join(f"{code}\t{count}" for code, count in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

