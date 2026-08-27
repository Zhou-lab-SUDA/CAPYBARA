from __future__ import annotations

import tempfile
import gzip
import unittest
from pathlib import Path

from capybara.calling import VariantObservation, count_input_records, evaluate_markers
from capybara.classifier import MarkerCall, classify, marker_key
from capybara.database import Marker, load_database


ROOT = Path(__file__).resolve().parent.parent
DB = load_database(ROOT / "capydb" / "esl")


def reference_calls(depth: int = 20):
    return {
        marker_key(marker): MarkerCall(
            position=marker.position, observed_ref=marker.ref, observed_alt="",
            depth=depth, expected_alt_depth=0, allele_fraction=0.0,
            mapq=None, status="REFERENCE",
        )
        for marker in DB.markers
    }


def set_node_matches(calls, node: str, limit: int | None = None):
    markers = [m for m in DB.markers if m.node == node and m.marker_set == "CANONICAL"]
    for marker in markers[:limit]:
        calls[marker_key(marker)] = MarkerCall(
            position=marker.position, observed_ref=marker.ref, observed_alt=marker.alt,
            depth=20, expected_alt_depth=20, allele_fraction=1.0,
            mapq=None, status="MATCH",
        )


def set_marker_match(calls, node: str, position: int):
    marker = next(m for m in DB.markers if m.node == node and m.position == position)
    calls[marker_key(marker)] = MarkerCall(
        position=position, observed_ref=marker.ref, observed_alt=marker.alt,
        depth=20, expected_alt_depth=20, allele_fraction=1.0,
        mapq=None, status="MATCH",
    )


class AlleleEvaluationTests(unittest.TestCase):
    def test_fastq_record_count(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("@r1\nAC\n+\nII\n@r2\nGT\n+\nII\n")
            self.assertEqual(count_input_records([path], "isolate_reads"), 2)

    def test_wrong_alt_at_correct_position(self):
        marker = Marker("2.5.6", "2.5", "sublineage", 1733723, "T", "G", 3, 1, 1, 1, "CANONICAL", "ACTIVE")
        calls = evaluate_markers(
            [marker], {1733723: 20},
            {1733723: VariantObservation("T", ("A",), (1, 19), 20)}, 10, 0.9,
        )
        self.assertEqual(calls[marker_key(marker)].status, "WRONG_ALT")

    def test_multiallelic_expected_alt_depth(self):
        marker = Marker("2.5.6", "2.5", "sublineage", 1733723, "T", "G", 3, 1, 1, 1, "CANONICAL", "ACTIVE")
        calls = evaluate_markers(
            [marker], {1733723: 24},
            {1733723: VariantObservation("T", ("A", "G"), (1, 3, 20), 24)}, 10, 0.8,
        )
        call = calls[marker_key(marker)]
        self.assertEqual(call.expected_alt_depth, 20)
        self.assertAlmostEqual(call.allele_fraction, 20 / 24)
        self.assertEqual(call.status, "MATCH")


class HierarchicalRuleTests(unittest.TestCase):
    def test_256_fallback_without_258_override(self):
        calls = reference_calls()
        set_marker_match(calls, "2.5.6", 1733723)
        result, _ = classify(DB, calls, "x", "assembly", 100, 100)
        self.assertEqual(result.sublineage, "2.5.6")

    def test_258_overrides_256(self):
        calls = reference_calls()
        set_marker_match(calls, "2.5.6", 1733723)
        set_marker_match(calls, "2.5.8", 2860886)
        result, _ = classify(DB, calls, "x", "assembly", 100, 100)
        self.assertEqual(result.sublineage, "2.5.8")

    def test_sibling_conflict(self):
        calls = reference_calls()
        set_node_matches(calls, "2.5.5")
        set_marker_match(calls, "2.5.8", 2860886)
        result, _ = classify(DB, calls, "x", "assembly", 100, 100)
        self.assertEqual(result.status, "AMBIGUOUS")
        self.assertEqual(result.lineage, "2.5")

    def test_strong_incompatible_parent_blocks_weak_child(self):
        calls = reference_calls()
        set_node_matches(calls, "2.1")
        set_marker_match(calls, "2.5.6", 1733723)
        result, _ = classify(DB, calls, "x", "assembly", 100, 100)
        self.assertEqual(result.status, "CONFLICTING_MARKERS")
        self.assertNotEqual(result.sublineage, "2.5.6")

    def test_incomplete_child_stops_at_lineage(self):
        calls = reference_calls()
        set_node_matches(calls, "2.5")
        set_node_matches(calls, "2.5.5", limit=1)
        result, _ = classify(DB, calls, "x", "assembly", 100, 100)
        self.assertEqual(result.lineage, "2.5")
        self.assertEqual(result.sublineage, "NA")
        self.assertEqual(result.status, "UNRESOLVED")

    def test_mixed_metagenome(self):
        calls = reference_calls()
        set_node_matches(calls, "2.5.5")
        set_marker_match(calls, "2.5.8", 2860886)
        result, _ = classify(DB, calls, "x", "metagenomic", 100, 50)
        self.assertEqual(result.status, "MIXED")

    def test_non_esl(self):
        calls = {
            key: MarkerCall(value.position, value.observed_ref, value.observed_alt, 0, 0, 0, None, "UNCALLABLE")
            for key, value in reference_calls().items()
        }
        result, _ = classify(DB, calls, "x", "assembly", 100, 0)
        self.assertEqual(result.status, "NON_ESL")

    def test_missing_loci_are_insufficient_coverage(self):
        calls = {
            key: MarkerCall(value.position, value.observed_ref, value.observed_alt, 0, 0, 0, None, "UNCALLABLE")
            for key, value in reference_calls().items()
        }
        result, _ = classify(DB, calls, "x", "assembly", 100, 80)
        self.assertEqual(result.status, "INSUFFICIENT_COVERAGE")


if __name__ == "__main__":
    unittest.main()
