"""Command-line interface for CAPYBARA 2 hierarchical SNP typing."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict
from pathlib import Path

from . import __version__
from .calling import call_markers_from_inputs
from .classifier import classify, marker_key
from .database import load_database, validate_database
from .mash_gate import MashGateResult, evaluate_mash_gate


SUMMARY_FIELDS = [
    "sample", "mode", "esl", "gc", "lineage", "sublineage", "status",
    "confidence", "lineage_basis", "gc_support", "lineage_support",
    "sublineage_support", "callable_markers", "matched_markers",
    "conflicting_markers", "reads_total", "reads_mapped", "alternatives",
    "gate_method", "gate_result", "gate_distance", "gate_margin",
]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="CAPYBARA 2 hierarchical multi-SNP classifier for A. baumannii ESLs"
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--db", type=Path, default=Path(__file__).resolve().parent.parent / "capydb" / "esl")
    parser.add_argument("--validate-db", action="store_true")
    parser.add_argument("--assembly", type=Path, help="assembled genome FASTA/FASTA.GZ")
    parser.add_argument("-i", "--input", type=Path, help="backward-compatible single input")
    parser.add_argument("-1", "--read1", type=Path, help="read 1 FASTQ/FASTQ.GZ")
    parser.add_argument("-2", "--read2", type=Path, help="read 2 FASTQ/FASTQ.GZ")
    parser.add_argument("--metagenomic", action="store_true", help="use metagenomic evidence thresholds")
    parser.add_argument("-o", "--output", type=Path, default=Path("Capy"))
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--min-depth", type=int)
    parser.add_argument("--min-af", type=float)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--skip-mash-gate", action="store_true")
    parser.add_argument("--mash-max-distance", type=float, default=0.02)
    parser.add_argument("--mash-min-margin", type=float, default=0.0005)
    return parser


def resolve_inputs(args: argparse.Namespace) -> tuple[str, list[Path]]:
    supplied = int(args.assembly is not None) + int(args.input is not None) + int(args.read1 is not None)
    if supplied != 1:
        raise ValueError("provide exactly one of --assembly, -i/--input, or -1/--read1")
    if args.read2 is not None and args.read1 is None:
        raise ValueError("-2/--read2 requires -1/--read1")
    if args.assembly is not None:
        if args.metagenomic:
            raise ValueError("--assembly cannot be combined with --metagenomic")
        mode, inputs = "assembly", [args.assembly]
    elif args.read1 is not None:
        mode = "metagenomic" if args.metagenomic else "isolate_reads"
        inputs = [args.read1] + ([args.read2] if args.read2 is not None else [])
    else:
        mode = "metagenomic" if args.metagenomic else "assembly"
        inputs = [args.input]
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise ValueError(f"input files do not exist: {', '.join(missing)}")
    return mode, inputs


def default_thresholds(mode: str) -> tuple[int, float]:
    if mode == "assembly":
        return 1, 0.90
    if mode == "isolate_reads":
        return 10, 0.90
    return 2, 0.20


def run_cli(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    database = load_database(args.db, validate=False)
    errors, warnings = validate_database(database)
    if args.validate_db:
        for warning in warnings:
            print(f"WARNING\t{warning}")
        for error in errors:
            print(f"ERROR\t{error}")
        print(f"validation\terrors={len(errors)}\twarnings={len(warnings)}")
        return 1 if errors else 0
    if errors:
        raise ValueError("Invalid CAPYBARA database:\n- " + "\n- ".join(errors))
    for warning in warnings:
        print(f"[DB WARNING] {warning}", file=sys.stderr)

    mode, inputs = resolve_inputs(args)
    default_depth, default_af = default_thresholds(mode)
    min_depth = args.min_depth if args.min_depth is not None else default_depth
    min_af = args.min_af if args.min_af is not None else default_af
    gate = MashGateResult("NOT_RUN", "", "", 1.0, None, "not applicable")
    if mode == "assembly" and not args.skip_mash_gate:
        gate = evaluate_mash_gate(
            inputs[0], args.db.parent / "msh",
            max_distance=args.mash_max_distance, min_margin=args.mash_min_margin,
        )
    calls, mapping = call_markers_from_inputs(
        markers=database.markers, reference=database.reference, inputs=inputs,
        mode=mode, threads=args.threads, min_depth=min_depth,
        min_allele_fraction=min_af, min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
    )
    sample = inputs[0].name
    result, evidence = classify(
        database, calls, sample=sample, mode=mode,
        reads_total=mapping.total, reads_mapped=mapping.mapped,
    )
    result.gate_method = "MASH" if gate.result not in {"NOT_RUN", "UNAVAILABLE"} else gate.result
    result.gate_result = gate.result
    result.gate_distance = f"{gate.distance:.8f}" if gate.result not in {"NOT_RUN", "UNAVAILABLE"} else ""
    result.gate_margin = f"{gate.margin:.8f}" if gate.margin is not None else ""
    if gate.result == "NON_ESL":
        result.esl, result.gc, result.lineage, result.sublineage = "NON_ESL", "NA", "NA", "NA"
        result.status, result.confidence, result.lineage_basis = "NON_ESL", "HIGH", "MASH_GATE"
        result.alternatives = gate.reason
    elif gate.result in {"GC1", "GC2"} and result.gc not in {"NA", gate.result}:
        result.status, result.confidence = "CONFLICTING_MARKERS", "LOW"
        result.alternatives = f"MASH={gate.result};BARCODE={result.gc}"

    summary_path = Path(f"{args.output}.summary.tsv")
    marker_path = Path(f"{args.output}.markers.tsv")
    evidence_path = Path(f"{args.output}.evidence.json")
    write_tsv(summary_path, [result.to_dict()], SUMMARY_FIELDS)

    marker_rows: list[dict[str, object]] = []
    for marker in database.markers:
        call = calls[marker_key(marker)]
        marker_rows.append(
            {
                "sample": sample, "node": marker.node, "position": marker.position,
                "expected_ref": marker.ref, "expected_alt": marker.alt,
                "observed_ref": call.observed_ref, "observed_alt": call.observed_alt or ".",
                "DP": call.depth, "expected_alt_depth": call.expected_alt_depth,
                "AF": f"{call.allele_fraction:.6f}", "MAPQ": f">={args.min_mapq}",
                "status": call.status, "marker_set": marker.marker_set,
            }
        )
    write_tsv(
        marker_path, marker_rows,
        [
            "sample", "node", "position", "expected_ref", "expected_alt",
            "observed_ref", "observed_alt", "DP", "expected_alt_depth", "AF",
            "MAPQ", "status", "marker_set",
        ],
    )
    payload = {
        "classifier_version": __version__, "thresholds": {
            "min_depth": min_depth, "min_af": min_af,
            "min_mapq": args.min_mapq, "min_baseq": args.min_baseq,
        },
        "mash_gate": asdict(gate),
        "classification": result.to_dict(),
        "node_evidence": {node: asdict(item) for node, item in evidence.items()},
    }
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    with evidence_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    print(f"[OK] {result.gc}>{result.lineage}>{result.sublineage} {result.status} {result.confidence}")
    print(f"[OK] wrote {summary_path}, {marker_path}, {evidence_path}")
    return 0


def main() -> None:
    try:
        raise SystemExit(run_cli())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
