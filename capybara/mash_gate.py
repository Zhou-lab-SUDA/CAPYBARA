"""Broad assembly gate using Mash distances to labelled reference sketches."""

from __future__ import annotations

import csv
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MashGateResult:
    result: str
    nearest_id: str
    nearest_label: str
    distance: float
    margin: float | None
    reason: str


def normalize_id(value: str) -> str:
    name = Path(value).name
    for suffix in (".fna.gz", ".fasta.gz", ".fa.gz", ".msh"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def load_labels(path: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            genome, label = line.rstrip().split("\t")
            labels[genome] = label
    return labels


def label_class(label: str) -> str:
    if label == "Clonal_complex_of_ST1":
        return "GC1"
    if label == "Clonal_complex_of_ST2":
        return "GC2"
    return "NON_ESL"


def evaluate_mash_gate(
    query: Path,
    mash_dir: Path,
    max_distance: float = 0.02,
    min_margin: float = 0.0005,
) -> MashGateResult:
    mash = shutil.which("mash")
    if mash is None:
        return MashGateResult("UNAVAILABLE", "", "", 1.0, None, "mash executable unavailable")
    sketch = mash_dir / "capy_all.msh"
    label_file = mash_dir / "cc.mash.list"
    if not sketch.exists() or not label_file.exists():
        return MashGateResult("UNAVAILABLE", "", "", 1.0, None, "combined Mash database unavailable")
    labels = load_labels(label_file)
    process = subprocess.run(
        [mash, "dist", str(sketch), str(query)], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if process.returncode != 0:
        raise RuntimeError(f"mash dist failed ({process.returncode}):\n{process.stderr}")
    rows: list[tuple[float, str, str, str]] = []
    for line in process.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        genome = normalize_id(fields[0])
        label = labels.get(genome, "UNKNOWN")
        rows.append((float(fields[2]), genome, label, label_class(label)))
    if not rows:
        return MashGateResult("NON_ESL", "", "", 1.0, None, "no Mash hits")
    rows.sort()
    distance, genome, label, nearest_class = rows[0]
    competing = [row[0] for row in rows[1:] if row[3] != nearest_class]
    margin = (min(competing) - distance) if competing else None
    if distance > max_distance:
        result, reason = "NON_ESL", "nearest distance exceeds absolute threshold"
    elif nearest_class == "NON_ESL":
        result, reason = "NON_ESL", "nearest labelled reference is outside ST1/ST2"
    elif margin is not None and margin < min_margin:
        result, reason = "AMBIGUOUS", "nearest-class margin below threshold"
    else:
        result, reason = nearest_class, "nearest reference and margin support ESL gate"
    return MashGateResult(result, genome, label, distance, margin, reason)
