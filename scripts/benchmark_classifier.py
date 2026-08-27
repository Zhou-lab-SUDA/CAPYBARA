#!/usr/bin/env python3
"""Replay old, exact-ALT single-SNP, and new CAPYBARA on tree tip states."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from capybara.classifier import MarkerCall, classify, marker_key  # noqa: E402
from capybara.database import load_database  # noqa: E402
from newick_utils import Node, postorder, preorder, read_newick  # noqa: E402


def open_tsv(path: Path):
    handle = gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(encoding="utf-8", newline="")
    return handle, csv.DictReader(handle, delimiter="\t")


def canonical_tip(name: str) -> str:
    for suffix in (".result.fastq", ".fna"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_old_panel(path: Path) -> list[tuple[int, str, str, str]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            position, label, ref, alt = line.rstrip().split("\t")
            rows.append((int(position), label, ref, alt))
    return rows


def reconstruct_alleles(
    tree_path: Path, mutation_path: Path, positions: set[int]
) -> tuple[list[str], dict[int, dict[str, int]]]:
    root = read_newick(tree_path)
    nodes = list(preorder(root))
    named = {node.name: node for node in nodes if node.name}
    depth = {root: 0}
    for node in nodes:
        for child in node.children:
            depth[child] = depth[node] + 1
    tips = [node for node in nodes if node.is_leaf]
    tip_ids = [canonical_tip(node.name) for node in tips]
    tip_index = {sample: index for index, sample in enumerate(tip_ids)}
    bits: dict[Node, int] = {}
    for node in postorder(root):
        if node.is_leaf:
            bits[node] = 1 << tip_index[canonical_tip(node.name)]
        else:
            bits[node] = 0
            for child in node.children:
                bits[node] |= bits[child]
    all_bits = bits[root]
    events_by_position: defaultdict[int, list[dict[str, str]]] = defaultdict(list)
    handle, reader = open_tsv(mutation_path)
    for row in reader:
        position = int(row["position"])
        if position in positions:
            events_by_position[position].append(row)
    handle.close()

    states_by_position: dict[int, dict[str, int]] = {}
    for position in positions:
        events = events_by_position.get(position, [])
        if not events:
            continue
        event_nodes = {named[row["node"]] for row in events}
        top = []
        for row in events:
            ancestor = named[row["node"]].parent
            while ancestor is not None and ancestor not in event_nodes:
                ancestor = ancestor.parent
            if ancestor is None:
                top.append(row)
        root_states = {row["ref"] for row in top}
        if len(root_states) != 1:
            continue
        states = {base: 0 for base in "ACGT"}
        states[next(iter(root_states))] = all_bits
        valid = True
        for row in sorted(events, key=lambda item: depth[named[item["node"]]]):
            event_bits = bits[named[row["node"]]]
            current = [base for base, state_bits in states.items() if state_bits & event_bits == event_bits]
            if current != [row["ref"]]:
                valid = False
                break
            for base in states:
                states[base] &= ~event_bits
            states[row["alt"]] |= event_bits
        if valid:
            states_by_position[position] = states
    return tip_ids, states_by_position


def allele_for(index: int, position: int, states: dict[int, dict[str, int]]) -> str:
    bit = 1 << index
    for base, sample_bits in states[position].items():
        if sample_bits & bit:
            return base
    return "N"


def old_predict(
    index: int,
    states: dict[int, dict[str, int]],
    reference: str,
    lineage_panel: list[tuple[int, str, str, str]],
    variant_panel: list[tuple[int, str, str, str]],
    exact_alt: bool,
) -> tuple[str, str, str]:
    def supported(panel):
        labels = []
        for position, label, _expected_ref, expected_alt in panel:
            if position not in states:
                continue
            allele = allele_for(index, position, states)
            present = allele == expected_alt if exact_alt else allele != reference[position - 1]
            if present:
                labels.append(label)
        return labels

    variants = supported(variant_panel)
    if variants:
        sublineage = sorted(variants)[0]
        lineage = ".".join(sublineage.split(".")[:2])
        return ("GC1" if lineage.startswith("1.") else "GC2"), lineage, sublineage
    lineages = supported(lineage_panel)
    if lineages:
        lineage = sorted(lineages)[0]
        return ("GC1" if lineage.startswith("1.") else "GC2"), lineage, "NA"
    return "NA", "NA", "NA"


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--mutations", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--old-lineage", type=Path, required=True)
    parser.add_argument("--old-variant", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    database = load_database(args.db)
    old_lineage = read_old_panel(args.old_lineage)
    old_variant = read_old_panel(args.old_variant)
    positions = {marker.position for marker in database.markers}
    positions.update(row[0] for row in old_lineage + old_variant)
    tip_ids, states = reconstruct_alleles(args.tree, args.mutations, positions)
    missing_positions = positions - states.keys()
    if missing_positions:
        raise ValueError(f"could not reconstruct marker positions: {sorted(missing_positions)[:20]}")
    metadata_handle, metadata_reader = open_tsv(args.metadata)
    metadata = {row["sample_id"]: row for row in metadata_reader}
    metadata_handle.close()
    reference = "".join(
        line.strip() for line in database.reference.open(encoding="utf-8")
        if not line.startswith(">")
    ).upper()

    schemes = {"old": [], "bugfixed": [], "new": []}
    fields = [
        "sample", "truth_gc", "truth_lineage", "truth_sublineage",
        "pred_gc", "pred_lineage", "pred_sublineage", "status",
    ]
    for index, sample in enumerate(tip_ids):
        truth = metadata[sample]
        truth_sublineage = truth["sublineage"] or "NA"
        for scheme, exact in (("old", False), ("bugfixed", True)):
            gc, lineage, sublineage = old_predict(
                index, states, reference, old_lineage, old_variant, exact,
            )
            schemes[scheme].append(
                {
                    "sample": sample, "truth_gc": truth["gc"],
                    "truth_lineage": truth["lineage"], "truth_sublineage": truth_sublineage,
                    "pred_gc": gc, "pred_lineage": lineage, "pred_sublineage": sublineage,
                    "status": "PASS" if sublineage != "NA" else "UNRESOLVED",
                }
            )

        calls = {}
        for marker in database.markers:
            allele = allele_for(index, marker.position, states)
            if allele == marker.alt:
                status = "MATCH"
            elif allele == marker.ref:
                status = "REFERENCE"
            else:
                status = "WRONG_ALT"
            calls[marker_key(marker)] = MarkerCall(
                marker.position, marker.ref, allele if allele != marker.ref else "",
                1, int(status == "MATCH"), float(status == "MATCH"), None, status,
            )
        result, _ = classify(database, calls, sample, "assembly", 1, 1)
        schemes["new"].append(
            {
                "sample": sample, "truth_gc": truth["gc"],
                "truth_lineage": truth["lineage"], "truth_sublineage": truth_sublineage,
                "pred_gc": result.gc, "pred_lineage": result.lineage,
                "pred_sublineage": result.sublineage, "status": result.status,
            }
        )

    summary_rows = []
    for scheme, rows in schemes.items():
        write_tsv(args.output_dir / f"benchmark_{scheme}.tsv", rows, fields)
        total = len(rows)
        for level in ("gc", "lineage", "sublineage"):
            truth_key, pred_key = f"truth_{level}", f"pred_{level}"
            correct = sum(row[truth_key] == row[pred_key] for row in rows)
            resolved = [row for row in rows if row[pred_key] != "NA"]
            labels = sorted({row[truth_key] for row in rows if row[truth_key] != "NA"})
            precision_values = []
            recall_values = []
            f1_values = []
            for label in labels:
                tp = sum(row[truth_key] == label and row[pred_key] == label for row in rows)
                fp = sum(row[truth_key] != label and row[pred_key] == label for row in rows)
                fn = sum(row[truth_key] == label and row[pred_key] != label for row in rows)
                precision = tp / (tp + fp) if tp + fp else 0.0
                recall = tp / (tp + fn) if tp + fn else 0.0
                f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
                precision_values.append(precision)
                recall_values.append(recall)
                f1_values.append(f1)
            wrong_resolved = sum(
                row[pred_key] != "NA" and row[pred_key] != row[truth_key] for row in rows
            )
            false_negative = sum(
                row[truth_key] != "NA" and row[pred_key] == "NA" for row in rows
            )
            false_terminal = (
                wrong_resolved / total
                if level == "sublineage" else ""
            )
            summary_rows.append(
                {
                    "scheme": scheme, "level": level, "n": total,
                    "accuracy": f"{correct / total:.8f}",
                    "macro_precision": f"{sum(precision_values) / len(precision_values):.8f}",
                    "macro_recall": f"{sum(recall_values) / len(recall_values):.8f}",
                    "macro_f1": f"{sum(f1_values) / len(f1_values):.8f}",
                    "resolved_rate": f"{len(resolved) / total:.8f}",
                    "unresolved_rate": f"{1 - len(resolved) / total:.8f}",
                    "ambiguous_rate": f"{sum(row['status'] in {'AMBIGUOUS', 'MIXED'} for row in rows) / total:.8f}",
                    "false_positive_rate": f"{wrong_resolved / total:.8f}",
                    "false_negative_rate": f"{false_negative / total:.8f}",
                    "false_terminal_assignment_rate": (
                        f"{false_terminal:.8f}" if false_terminal != "" else ""
                    ),
                }
            )
            matrix = Counter((row[truth_key], row[pred_key]) for row in rows)
            write_tsv(
                args.output_dir / "confusion_matrices" / f"{scheme}.{level}.tsv",
                [
                    {"truth": truth_label, "prediction": prediction, "n": count}
                    for (truth_label, prediction), count in sorted(matrix.items())
                ],
                ["truth", "prediction", "n"],
            )
    write_tsv(
        args.output_dir / "benchmark_summary.tsv", summary_rows,
        [
            "scheme", "level", "n", "accuracy", "macro_precision", "macro_recall",
            "macro_f1", "resolved_rate", "unresolved_rate", "ambiguous_rate",
            "false_positive_rate", "false_negative_rate",
            "false_terminal_assignment_rate",
        ],
    )
    print(f"samples={len(tip_ids)} output={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
