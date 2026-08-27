#!/usr/bin/env python3
"""Discover tree-defined canonical SNPs on MRCA branches for CAPYBARA nodes."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

from newick_utils import Node, postorder, preorder, read_newick


def open_tsv(path: Path):
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open(encoding="utf-8", newline="")
    return handle, csv.DictReader(handle, delimiter="\t")


def canonical_tip(name: str) -> str:
    for suffix in (".result.fastq", ".fna"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float) -> str:
    return f"{value:.8f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--mutations", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--hierarchy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mrca-output", type=Path, required=True)
    parser.add_argument("--node-summary", type=Path, required=True)
    parser.add_argument("--min-sensitivity", type=float, default=0.98)
    parser.add_argument("--min-specificity", type=float, default=0.999)
    parser.add_argument(
        "--max-homoplasy",
        type=int,
        default=10,
        help="Hard upper bound for selectable marker events (default: 10)",
    )
    parser.add_argument(
        "--allow-outside-alt",
        action="store_true",
        help=(
            "Allow expected ALT outside the target node. By default any reconstructed "
            "outside occurrence is a hard cross-node conflict, including occurrences "
            "caused by events that are themselves excluded for high homoplasy."
        ),
    )
    args = parser.parse_args()

    root = read_newick(args.tree)
    nodes = list(preorder(root))
    named_nodes = {node.name: node for node in nodes if node.name}
    depth: dict[Node, int] = {root: 0}
    for node in nodes:
        for child in node.children:
            depth[child] = depth[node] + 1

    tips = [node for node in nodes if node.is_leaf]
    tip_ids = [canonical_tip(node.name) for node in tips]
    if len(tip_ids) != len(set(tip_ids)):
        raise ValueError("Canonicalized tree tip identifiers are not unique")
    tip_index = {sample_id: index for index, sample_id in enumerate(tip_ids)}
    node_bits: dict[Node, int] = {}
    for node in postorder(root):
        if node.is_leaf:
            node_bits[node] = 1 << tip_index[canonical_tip(node.name)]
        else:
            bits = 0
            for child in node.children:
                bits |= node_bits[child]
            node_bits[node] = bits
    all_bits = node_bits[root]

    metadata_handle, metadata_reader = open_tsv(args.metadata)
    metadata = list(metadata_reader)
    metadata_handle.close()
    metadata_by_id = {row["sample_id"]: row for row in metadata}
    if set(metadata_by_id).intersection(tip_index) != set(tip_index):
        missing = sorted(set(tip_index) - set(metadata_by_id))
        raise ValueError(f"Tree tips absent from metadata: {missing[:10]}")

    hierarchy_handle, hierarchy_reader = open_tsv(args.hierarchy)
    hierarchy = list(hierarchy_reader)
    hierarchy_handle.close()
    hierarchy_by_node = {row["node"]: row for row in hierarchy}
    discoverable = [
        row for row in hierarchy
        if row["classifiable"] == "true" and row["level"] in {"GC", "lineage", "sublineage"}
    ]

    group_bits: dict[str, int] = {}
    for item in discoverable:
        label = item["node"]
        bits = 0
        for sample_id, index in tip_index.items():
            row = metadata_by_id[sample_id]
            belongs = (
                (item["level"] == "GC" and row["gc"] == label)
                or (item["level"] == "lineage" and row["lineage"] == label)
                or (item["level"] == "sublineage" and row["sublineage"] == label)
            )
            if belongs:
                bits |= 1 << index
        if not bits:
            raise ValueError(f"Hierarchy node has no tree tips: {label}")
        group_bits[label] = bits

    def find_mrca(target_bits: int) -> Node:
        node = root
        while True:
            covering = [child for child in node.children if node_bits[child] & target_bits == target_bits]
            if not covering:
                return node
            if len(covering) != 1:
                raise AssertionError("Disjoint children cannot both contain all target tips")
            node = covering[0]

    mrca_rows: list[dict[str, object]] = []
    mrca_by_label: dict[str, Node] = {}
    for item in discoverable:
        label = item["node"]
        target = group_bits[label]
        mrca = find_mrca(target)
        descendant_n = node_bits[mrca].bit_count()
        target_n = target.bit_count()
        purity = target_n / descendant_n
        mrca_by_label[label] = mrca
        mrca_rows.append(
            {
                "node": label,
                "parent": item["parent"],
                "level": item["level"],
                "target_n": target_n,
                "mrca_node": mrca.name or "<root>",
                "mrca_descendant_n": descendant_n,
                "mrca_purity": fmt(purity),
                "monophyletic": str(target == node_bits[mrca]).lower(),
                "mrca_depth": depth[mrca],
            }
        )
    write_tsv(
        args.mrca_output,
        mrca_rows,
        [
            "node", "parent", "level", "target_n", "mrca_node", "mrca_descendant_n",
            "mrca_purity", "monophyletic", "mrca_depth",
        ],
    )

    labels_by_mrca: defaultdict[str, list[str]] = defaultdict(list)
    for label, node in mrca_by_label.items():
        if node.name:
            labels_by_mrca[node.name].append(label)
    mrca_names = set(labels_by_mrca)

    candidate_events: list[dict[str, str]] = []
    mutation_handle, mutation_reader = open_tsv(args.mutations)
    for row in mutation_reader:
        if row["node"] in mrca_names:
            candidate_events.append(row)
    mutation_handle.close()
    candidate_positions = {int(row["position"]) for row in candidate_events}

    events_by_position: defaultdict[int, list[dict[str, str]]] = defaultdict(list)
    mutation_handle, mutation_reader = open_tsv(args.mutations)
    for row in mutation_reader:
        position = int(row["position"])
        if position in candidate_positions:
            events_by_position[position].append(row)
    mutation_handle.close()

    allele_bits_by_position: dict[int, dict[str, int]] = {}
    position_notes: defaultdict[int, list[str]] = defaultdict(list)
    for position, events in events_by_position.items():
        event_nodes = {named_nodes[row["node"]] for row in events}
        top_events: list[dict[str, str]] = []
        for row in events:
            ancestor = named_nodes[row["node"]].parent
            has_ancestor_event = False
            while ancestor is not None:
                if ancestor in event_nodes:
                    has_ancestor_event = True
                    break
                ancestor = ancestor.parent
            if not has_ancestor_event:
                top_events.append(row)
        root_states = {row["ref"] for row in top_events}
        if len(root_states) != 1:
            position_notes[position].append("INCONSISTENT_ROOT_STATE")
            continue
        root_state = next(iter(root_states))
        states = {base: 0 for base in "ACGT"}
        states[root_state] = all_bits
        valid = True
        for row in sorted(events, key=lambda item: depth[named_nodes[item["node"]]]):
            event_node = named_nodes[row["node"]]
            event_bits = node_bits[event_node]
            current_bases = [base for base, bits in states.items() if bits & event_bits == event_bits]
            if current_bases != [row["ref"]]:
                position_notes[position].append(
                    f"STATE_CONFLICT_AT_{row['node']}_{row['ref']}->{row['alt']}"
                )
                valid = False
                break
            for base in states:
                states[base] &= ~event_bits
            states[row["alt"]] |= event_bits
        if valid:
            allele_bits_by_position[position] = states

    hierarchy_children: defaultdict[str, list[str]] = defaultdict(list)
    for row in hierarchy:
        hierarchy_children[row["parent"]].append(row["node"])

    candidates: list[dict[str, object]] = []
    for event in candidate_events:
        position = int(event["position"])
        states = allele_bits_by_position.get(position)
        if states is None:
            continue
        for label in labels_by_mrca[event["node"]]:
            item = hierarchy_by_node[label]
            target = group_bits[label]
            outside = all_bits ^ target
            alt_bits = states[event["alt"]]
            target_n = target.bit_count()
            non_target_n = outside.bit_count()
            target_alt_n = (alt_bits & target).bit_count()
            outside_alt_n = (alt_bits & outside).bit_count()
            sensitivity = target_alt_n / target_n
            specificity = 1 - (outside_alt_n / non_target_n) if non_target_n else 1.0
            sibling_labels = [
                sibling for sibling in hierarchy_children[item["parent"]]
                if sibling != label and sibling in group_bits
            ]
            sister_bits = 0
            for sibling in sibling_labels:
                sister_bits |= group_bits[sibling]
            sister_alt_n = (alt_bits & sister_bits).bit_count()
            elsewhere: list[str] = []
            for other_label, other_bits in group_bits.items():
                if other_label == label:
                    continue
                overlap = (alt_bits & other_bits).bit_count()
                if overlap and not (other_bits & target == other_bits):
                    elsewhere.append(f"{other_label}:{overlap}/{other_bits.bit_count()}")
            outside_sublineages: dict[str, int] = {}
            outside_alt_bits = alt_bits & outside
            for other_label, other_group_bits in group_bits.items():
                other_item = hierarchy_by_node[other_label]
                if other_item["level"] != "sublineage" or other_label == label:
                    continue
                overlap = (outside_alt_bits & other_group_bits).bit_count()
                if overlap:
                    outside_sublineages[other_label] = overlap
            notes = list(position_notes[position])
            if int(event["homoplasy"]) > 1:
                notes.append("HOMOPLASY_GT_1")
            if int(event["homoplasy"]) > args.max_homoplasy:
                notes.append("HOMOPLASY_GT_MAX")
            if sister_alt_n:
                notes.append("OBSERVED_IN_SISTER")
            if outside_alt_n:
                notes.append("CROSS_NODE_CONFLICT")

            rejection_reasons: list[str] = []
            if int(event["homoplasy"]) > args.max_homoplasy:
                rejection_reasons.append("HOMOPLASY_GT_MAX")
            if sensitivity < args.min_sensitivity:
                rejection_reasons.append("LOW_SENSITIVITY")
            if specificity < args.min_specificity:
                rejection_reasons.append("LOW_SPECIFICITY")
            if outside_alt_n and not args.allow_outside_alt:
                rejection_reasons.append("CROSS_NODE_CONFLICT")
            status = "PASS" if not rejection_reasons else "REJECT"
            candidates.append(
                {
                    "node": label,
                    "parent": item["parent"],
                    "level": item["level"],
                    "mrca_node": event["node"],
                    "position": position,
                    "ref": event["ref"],
                    "alt": event["alt"],
                    "homoplasy": int(event["homoplasy"]),
                    "sensitivity": fmt(sensitivity),
                    "specificity": fmt(specificity),
                    "target_alt_frequency": fmt(sensitivity),
                    "non_target_alt_frequency": fmt(1 - specificity),
                    # A branch-mutation reconstruction infers an allele for every tip;
                    # it does not measure whether the original isolate locus was callable.
                    "missing_rate": "NA",
                    "callability_basis": "ANCESTRAL_RECONSTRUCTION",
                    "target_n": target_n,
                    "non_target_n": non_target_n,
                    "supporting_genomes": target_alt_n,
                    "contradictory_genomes": target_n - target_alt_n,
                    "outside_alt_n": outside_alt_n,
                    "sister_alt_n": sister_alt_n,
                    "elsewhere_nodes": ";".join(elsewhere),
                    "outside_sublineages": ";".join(
                        f"{name}:{count}" for name, count in sorted(outside_sublineages.items())
                    ),
                    "candidate_status": status,
                    "rejection_reasons": ";".join(rejection_reasons),
                    "marker_rank": 0,
                    "notes": ";".join(notes),
                }
            )

    by_label: defaultdict[str, list[dict[str, object]]] = defaultdict(list)
    for row in candidates:
        by_label[str(row["node"])].append(row)
    for rows in by_label.values():
        rows.sort(
            key=lambda row: (
                row["candidate_status"] != "PASS",
                -float(row["specificity"]),
                -float(row["sensitivity"]),
                int(row["homoplasy"]),
                int(row["position"]),
            )
        )
        for rank, row in enumerate(rows, 1):
            row["marker_rank"] = rank
    candidates.sort(key=lambda row: (row["level"], row["node"], int(row["marker_rank"])))
    candidate_fields = [
        "node", "parent", "level", "mrca_node", "position", "ref", "alt", "homoplasy",
        "sensitivity", "specificity", "target_alt_frequency", "non_target_alt_frequency",
        "missing_rate", "callability_basis", "target_n", "non_target_n", "supporting_genomes",
        "contradictory_genomes", "outside_alt_n", "sister_alt_n", "elsewhere_nodes",
        "outside_sublineages",
        "candidate_status", "rejection_reasons", "marker_rank", "notes",
    ]
    write_tsv(args.output, candidates, candidate_fields)

    summary_rows: list[dict[str, object]] = []
    mrca_lookup = {row["node"]: row for row in mrca_rows}
    for item in discoverable:
        label = item["node"]
        node_candidates = by_label.get(label, [])
        passing = [row for row in node_candidates if row["candidate_status"] == "PASS"]
        homoplasy_one = [row for row in passing if int(row["homoplasy"]) == 1]
        summary_rows.append(
            {
                "node": label,
                "parent": item["parent"],
                "level": item["level"],
                "target_n": group_bits[label].bit_count(),
                "mrca_node": mrca_lookup[label]["mrca_node"],
                "monophyletic": mrca_lookup[label]["monophyletic"],
                "mrca_purity": mrca_lookup[label]["mrca_purity"],
                "branch_mutation_n": len(node_candidates),
                "passing_candidate_n": len(passing),
                "passing_homoplasy_one_n": len(homoplasy_one),
                "max_homoplasy": args.max_homoplasy,
                "requires_zero_outside_alt": str(not args.allow_outside_alt).lower(),
                "has_at_least_3_robust": str(len(passing) >= 3).lower(),
            }
        )
    write_tsv(
        args.node_summary,
        summary_rows,
        [
            "node", "parent", "level", "target_n", "mrca_node", "monophyletic",
            "mrca_purity", "branch_mutation_n", "passing_candidate_n",
            "passing_homoplasy_one_n", "max_homoplasy", "requires_zero_outside_alt",
            "has_at_least_3_robust",
        ],
    )
    print(f"nodes={len(discoverable)} candidates={len(candidates)} passing={sum(r['candidate_status'] == 'PASS' for r in candidates)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
