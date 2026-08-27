#!/usr/bin/env python3
"""Parse ancestral branch-mutation records and preserve their tree semantics."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


MUTATION_RE = re.compile(r"^([ACGT])->([ACGT])$")


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def canonical_sample_id(node: str) -> tuple[str, str]:
    if re.fullmatch(r"N_\d+", node):
        return "internal", ""
    if node.endswith(".result.fastq"):
        return "leaf", node[: -len(".result.fastq")]
    if node.endswith(".fna"):
        return "leaf", node[: -len(".fna")]
    return "unknown", ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--missing-regions", type=Path, required=True)
    parser.add_argument("--sequence-lengths", type=Path, required=True)
    parser.add_argument("--node-stats", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    args = parser.parse_args()

    sequence_lengths: list[dict[str, object]] = []
    missing_regions: list[dict[str, object]] = []
    node_counts: Counter[str] = Counter()
    node_classes: dict[str, str] = {}
    node_samples: dict[str, str] = {}
    sequence_counts: Counter[str] = Counter()
    homoplasy_counts: Counter[int] = Counter()
    mutation_counts: Counter[str] = Counter()
    duplicate_keys: Counter[tuple[str, str, int]] = Counter()
    invalid_rows = 0
    data_rows = 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.input, "r") as source, open_text(args.output, "w") as target:
        writer = csv.DictWriter(
            target,
            fieldnames=[
                "node", "node_type", "sample_id", "sequence", "position",
                "homoplasy", "ref", "alt", "mutation",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        data_started = False
        for line_number, line in enumerate(source, 1):
            stripped = line.rstrip("\r\n")
            if stripped.startswith("## Sequence_length:"):
                fields = stripped.split()
                if len(fields) != 4:
                    raise ValueError(f"Malformed Sequence_length at line {line_number}")
                sequence_lengths.append(
                    {"sequence": fields[2], "length": int(fields[3]), "source_line": line_number}
                )
                continue
            if stripped.startswith("## Missing_region:"):
                fields = stripped.split()
                if len(fields) != 5:
                    raise ValueError(f"Malformed Missing_region at line {line_number}")
                missing_regions.append(
                    {
                        "sequence": fields[2], "start": int(fields[3]), "end": int(fields[4]),
                        "source_line": line_number,
                    }
                )
                continue
            if stripped == "#Node\t#Seq\t#Site\t#Homoplasy\t#Mutation":
                data_started = True
                continue
            if not stripped:
                continue
            if not data_started:
                raise ValueError(f"Unexpected content before header at line {line_number}")
            fields = stripped.split("\t")
            if len(fields) != 5:
                invalid_rows += 1
                continue
            node, sequence, position_text, homoplasy_text, mutation = fields
            match = MUTATION_RE.fullmatch(mutation)
            try:
                position = int(position_text)
                homoplasy = int(homoplasy_text)
            except ValueError:
                invalid_rows += 1
                continue
            if not match or position < 1 or homoplasy < 1:
                invalid_rows += 1
                continue
            ref, alt = match.groups()
            node_type, sample_id = canonical_sample_id(node)
            writer.writerow(
                {
                    "node": node,
                    "node_type": node_type,
                    "sample_id": sample_id,
                    "sequence": sequence,
                    "position": position,
                    "homoplasy": homoplasy,
                    "ref": ref,
                    "alt": alt,
                    "mutation": mutation,
                }
            )
            data_rows += 1
            node_counts[node] += 1
            node_classes[node] = node_type
            node_samples[node] = sample_id
            sequence_counts[sequence] += 1
            homoplasy_counts[homoplasy] += 1
            mutation_counts[mutation] += 1
            duplicate_keys[(node, sequence, position)] += 1

    sequence_length_by_name = {str(row["sequence"]): int(row["length"]) for row in sequence_lengths}
    out_of_range = 0
    # The normalized output is streamed above; verify ranges in a second streaming pass.
    with open_text(args.output, "r") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length = sequence_length_by_name.get(row["sequence"])
            if length is None or int(row["position"]) > length:
                out_of_range += 1

    def write_plain_tsv(path: Path, rows, fields: list[str]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    write_plain_tsv(
        args.sequence_lengths,
        sequence_lengths,
        ["sequence", "length", "source_line"],
    )
    write_plain_tsv(
        args.missing_regions,
        missing_regions,
        ["sequence", "start", "end", "source_line"],
    )
    node_rows = [
        {
            "node": node,
            "node_type": node_classes[node],
            "sample_id": node_samples[node],
            "mutation_events": count,
        }
        for node, count in sorted(node_counts.items())
    ]
    write_plain_tsv(
        args.node_stats,
        node_rows,
        ["node", "node_type", "sample_id", "mutation_events"],
    )
    stats = {
        "input": str(args.input.resolve()),
        "semantics": "ancestral_branch_mutations_requires_matching_named_tree",
        "sequence_lengths": sequence_length_by_name,
        "missing_region_count": len(missing_regions),
        "mutation_rows": data_rows,
        "invalid_rows": invalid_rows,
        "distinct_nodes_with_mutations": len(node_counts),
        "node_type_counts": dict(sorted(Counter(node_classes.values()).items())),
        "sequence_counts": dict(sorted(sequence_counts.items())),
        "homoplasy_one_rows": homoplasy_counts[1],
        "homoplasy_one_rate": homoplasy_counts[1] / data_rows if data_rows else 0,
        "mutation_counts": dict(sorted(mutation_counts.items())),
        "duplicate_node_sequence_position_rows": sum(
            count - 1 for count in duplicate_keys.values() if count > 1
        ),
        "out_of_range_rows": out_of_range,
    }
    args.stats.parent.mkdir(parents=True, exist_ok=True)
    args.stats.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(stats, sort_keys=True))
    return 1 if invalid_rows or out_of_range else 0


if __name__ == "__main__":
    raise SystemExit(main())

