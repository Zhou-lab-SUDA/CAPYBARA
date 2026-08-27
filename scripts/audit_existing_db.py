#!/usr/bin/env python3
"""Audit the legacy CAPYBARA barcode sources without importing configure.py."""

from __future__ import annotations

import argparse
import ast
import csv
from collections import defaultdict
from pathlib import Path


PANELS = {
    "hc1030": ("capydb/hc1030.snp", "hc1030_snp", "capydb/hc1030_ref.fna"),
    "lineage": ("capydb/esl/lineage.SNP", "lineage_snp", "capydb/esl/esl_ref.fna"),
    "variant": ("capydb/esl/variant.SNP", "variant_snp", "capydb/esl/esl_ref.fna"),
}


def read_fasta(path: Path) -> tuple[str, str]:
    name = ""
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    break
                name = line[1:].split()[0]
            elif name:
                sequence.append(line.upper())
    if not name or not sequence:
        raise ValueError(f"Invalid FASTA: {path}")
    return name, "".join(sequence)


def read_panel(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields or all(not value.strip() for value in fields):
                continue
            if len(fields) != 4:
                raise ValueError(f"{path}:{line_number}: expected 4 columns, got {len(fields)}")
            position, node, ref, alt = fields
            rows.append(
                {
                    "position": int(position),
                    "node": node.strip(),
                    "ref": ref.strip().upper(),
                    "alt": alt.strip().upper(),
                    "line_number": line_number,
                }
            )
    return rows


def literal_dict_assignment(path: Path, variable: str) -> dict[int, list[str]]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for statement in tree.body:
        if not isinstance(statement, (ast.Assign, ast.AnnAssign)):
            continue
        targets = statement.targets if isinstance(statement, ast.Assign) else [statement.target]
        if any(isinstance(target, ast.Name) and target.id == variable for target in targets):
            value = ast.literal_eval(statement.value)
            return {int(key): list(item) for key, item in value.items()}
    raise ValueError(f"No literal assignment for {variable} in {path}")


def parent_for(node: str, panel: str) -> str:
    if panel == "hc1030":
        return "ESL"
    if panel == "lineage":
        return "GC1" if node.startswith("1.") else "GC2" if node.startswith("2.") else "UNKNOWN"
    return node.rsplit(".", 1)[0] if "." in node else "UNKNOWN"


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output_dir = (args.output_dir or repo / "analysis").resolve()

    configure_path = repo / "configure.py"
    lineage_nodes: set[str] = set()
    audited: list[dict[str, object]] = []
    conflicts: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    loaded_panels: dict[str, list[dict[str, object]]] = {}
    configured_panels: dict[str, dict[int, list[str]]] = {}
    for panel, (relative_panel, variable, relative_reference) in PANELS.items():
        rows = read_panel(repo / relative_panel)
        configured = literal_dict_assignment(configure_path, variable)
        _, reference = read_fasta(repo / relative_reference)
        loaded_panels[panel] = rows
        configured_panels[panel] = configured
        if panel == "lineage":
            lineage_nodes = {str(row["node"]) for row in rows}

        by_allele: dict[tuple[int, str, str], list[str]] = defaultdict(list)
        for row in rows:
            by_allele[(int(row["position"]), str(row["ref"]), str(row["alt"]))].append(str(row["node"]))

        for row in rows:
            position = int(row["position"])
            expected = [str(row["node"]), str(row["ref"]), str(row["alt"])]
            configured_value = configured.get(position)
            duplicate_nodes = sorted(set(by_allele[(position, str(row["ref"]), str(row["alt"]))]))
            parent = parent_for(str(row["node"]), panel)
            parent_present = panel != "variant" or parent in lineage_nodes
            reference_base = reference[position - 1] if 1 <= position <= len(reference) else "OUT_OF_RANGE"
            audited.append(
                {
                    "panel": panel,
                    "node": row["node"],
                    "position": position,
                    "ref": row["ref"],
                    "alt": row["alt"],
                    "line_number": row["line_number"],
                    "reference_base": reference_base,
                    "reference_ref_match": reference_base == row["ref"],
                    "configured_value": "|".join(configured_value or []),
                    "external_matches_config": configured_value == expected,
                    "parent": parent,
                    "parent_marker_present": parent_present,
                    "duplicate_allele_nodes": "|".join(duplicate_nodes),
                    "duplicate_allele_conflict": len(duplicate_nodes) > 1,
                }
            )

        for (position, ref, alt), nodes in sorted(by_allele.items()):
            unique_nodes = sorted(set(nodes))
            if len(unique_nodes) > 1:
                conflicts.append(
                    {
                        "panel": panel,
                        "position": position,
                        "ref": ref,
                        "alt": alt,
                        "nodes": "|".join(unique_nodes),
                        "conflict_type": "same_allele_multiple_nodes",
                    }
                )

        external_by_position = {int(row["position"]): row for row in rows}
        summaries.append(
            {
                "panel": panel,
                "external_rows": len(rows),
                "external_unique_positions": len(external_by_position),
                "configured_positions": len(configured),
                "external_rows_not_active": sum(
                    configured.get(int(row["position"]))
                    != [str(row["node"]), str(row["ref"]), str(row["alt"])]
                    for row in rows
                ),
                "reference_ref_mismatches": sum(
                    not bool(row["reference_ref_match"])
                    for row in audited
                    if row["panel"] == panel
                ),
                "duplicate_allele_conflicts": sum(1 for row in conflicts if row["panel"] == panel),
            }
        )

    write_tsv(
        output_dir / "current_marker_audit.tsv",
        audited,
        [
            "panel", "node", "position", "ref", "alt", "line_number",
            "reference_base", "reference_ref_match", "configured_value",
            "external_matches_config", "parent", "parent_marker_present",
            "duplicate_allele_nodes", "duplicate_allele_conflict",
        ],
    )
    write_tsv(
        output_dir / "marker_conflicts.tsv",
        conflicts,
        ["panel", "position", "ref", "alt", "nodes", "conflict_type"],
    )
    write_tsv(
        output_dir / "current_scheme_summary.tsv",
        summaries,
        [
            "panel", "external_rows", "external_unique_positions", "configured_positions",
            "external_rows_not_active", "reference_ref_mismatches", "duplicate_allele_conflicts",
        ],
    )

    for summary in summaries:
        print("\t".join(str(summary[key]) for key in summary))
    return 1 if conflicts else 0


if __name__ == "__main__":
    raise SystemExit(main())
