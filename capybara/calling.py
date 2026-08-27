"""Assembly/read mapping and exact diagnostic-allele evaluation."""

from __future__ import annotations

import os
import gzip
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .classifier import MarkerCall, marker_key
from .database import Marker


@dataclass(frozen=True)
class MappingStatistics:
    total: int
    mapped: int


@dataclass(frozen=True)
class VariantObservation:
    ref: str
    alts: tuple[str, ...]
    allele_depths: tuple[int, ...]
    depth: int


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required executable is not available on PATH: {name}")
    return path


def run_checked(command: list[str], **kwargs) -> subprocess.CompletedProcess:
    process = subprocess.run(command, text=True, capture_output=True, **kwargs)
    if process.returncode != 0:
        raise RuntimeError(
            f"Command failed ({process.returncode}): {' '.join(command)}\n{process.stderr}"
        )
    return process


def map_query(
    reference: Path,
    inputs: list[Path],
    mode: str,
    threads: int,
    output_bam: Path,
) -> None:
    minimap2 = require_executable("minimap2")
    samtools = require_executable("samtools")
    preset = "asm5" if mode == "assembly" else "sr"
    map_command = [
        minimap2, "-x", preset, "--secondary=no", "-a", "-t", str(threads),
        str(reference), *(str(path) for path in inputs),
    ]
    view_command = [samtools, "view", "-b", "-F", "2308", "-"]
    sort_command = [samtools, "sort", "-@", str(threads), "-o", str(output_bam), "-"]
    with subprocess.Popen(map_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as mapper:
        assert mapper.stdout is not None
        with subprocess.Popen(
            view_command, stdin=mapper.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ) as viewer:
            mapper.stdout.close()
            assert viewer.stdout is not None
            sorter = subprocess.run(
                sort_command, stdin=viewer.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            viewer.stdout.close()
            viewer_stderr = viewer.stderr.read().decode(errors="replace") if viewer.stderr else ""
            viewer_code = viewer.wait()
        mapper_stderr = mapper.stderr.read().decode(errors="replace") if mapper.stderr else ""
        mapper_code = mapper.wait()
    if mapper_code != 0:
        raise RuntimeError(f"minimap2 failed ({mapper_code}):\n{mapper_stderr}")
    if viewer_code != 0:
        raise RuntimeError(f"samtools view failed ({viewer_code}):\n{viewer_stderr}")
    if sorter.returncode != 0:
        raise RuntimeError(
            f"samtools sort failed ({sorter.returncode}):\n{sorter.stderr.decode(errors='replace')}"
        )
    run_checked([samtools, "index", str(output_bam)])


def flagstat(bam: Path) -> MappingStatistics:
    samtools = require_executable("samtools")
    output = run_checked([samtools, "flagstat", str(bam)]).stdout
    total = mapped = 0
    for line in output.splitlines():
        if " in total " in line:
            total = int(line.split()[0])
        elif " mapped (" in line and " primary mapped (" not in line:
            mapped = int(line.split()[0])
    return MappingStatistics(total=total, mapped=mapped)


def count_input_records(inputs: list[Path], mode: str) -> int:
    total = 0
    for path in inputs:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            if mode == "assembly":
                total += sum(line.startswith(">") for line in handle)
            else:
                lines = sum(1 for _line in handle)
                if lines % 4:
                    raise ValueError(f"FASTQ line count is not divisible by four: {path}")
                total += lines // 4
    return total


def fasta_chromosome(reference: Path) -> str:
    with reference.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].strip().split()[0]
    raise ValueError(f"No FASTA header found: {reference}")


def write_bed(path: Path, chrom: str, positions: set[int]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for position in sorted(positions):
            handle.write(f"{chrom}\t{position - 1}\t{position}\n")


def depth_at_positions(
    bam: Path, bed: Path, min_mapq: int, min_baseq: int
) -> dict[int, int]:
    samtools = require_executable("samtools")
    output = run_checked(
        [
            samtools, "depth", "-a", "-b", str(bed), "-q", str(min_baseq),
            "-Q", str(min_mapq), str(bam),
        ]
    ).stdout
    result: dict[int, int] = {}
    for line in output.splitlines():
        _chrom, position, depth = line.split("\t")[:3]
        result[int(position)] = int(depth)
    return result


def call_variants(
    bam: Path,
    reference: Path,
    bed: Path,
    output_vcf: Path,
    min_mapq: int,
    min_baseq: int,
) -> None:
    bcftools = require_executable("bcftools")
    mpileup_command = [
        bcftools, "mpileup", "-Ou", "-f", str(reference), "-R", str(bed),
        "-q", str(min_mapq), "-Q", str(min_baseq), "-a", "FORMAT/AD,FORMAT/DP",
        str(bam),
    ]
    call_command = [
        bcftools, "call", "-mv", "--ploidy", "1", "-Ov", "-o", str(output_vcf),
    ]
    with subprocess.Popen(mpileup_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as pileup:
        assert pileup.stdout is not None
        caller = subprocess.run(
            call_command, stdin=pileup.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        pileup.stdout.close()
        pileup_stderr = pileup.stderr.read().decode(errors="replace") if pileup.stderr else ""
        pileup_code = pileup.wait()
    if pileup_code != 0:
        raise RuntimeError(f"bcftools mpileup failed ({pileup_code}):\n{pileup_stderr}")
    if caller.returncode != 0:
        raise RuntimeError(
            f"bcftools call failed ({caller.returncode}):\n{caller.stderr.decode(errors='replace')}"
        )


def parse_vcf(path: Path) -> dict[int, VariantObservation]:
    observations: dict[int, VariantObservation] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            position = int(fields[1])
            ref = fields[3].upper()
            alts = tuple(alt.upper() for alt in fields[4].split(",") if alt not in {".", "<*>", "<NON_REF>"})
            format_names = fields[8].split(":") if len(fields) > 8 else []
            sample_values = fields[9].split(":") if len(fields) > 9 else []
            values = dict(zip(format_names, sample_values))
            allele_depths: tuple[int, ...] = ()
            if values.get("AD") not in {None, "."}:
                try:
                    allele_depths = tuple(int(value) if value != "." else 0 for value in values["AD"].split(","))
                except ValueError:
                    allele_depths = ()
            try:
                depth = int(values.get("DP", "0"))
            except ValueError:
                depth = sum(allele_depths)
            observations[position] = VariantObservation(
                ref=ref, alts=alts, allele_depths=allele_depths, depth=depth,
            )
    return observations


def evaluate_markers(
    markers: list[Marker],
    depths: dict[int, int],
    variants: dict[int, VariantObservation],
    min_depth: int,
    min_allele_fraction: float,
) -> dict[tuple[int, str, str], MarkerCall]:
    calls: dict[tuple[int, str, str], MarkerCall] = {}
    for marker in markers:
        observation = variants.get(marker.position)
        depth = max(depths.get(marker.position, 0), observation.depth if observation else 0)
        observed_ref = observation.ref if observation else marker.ref
        observed_alt = ",".join(observation.alts) if observation else ""
        expected_alt_depth = 0
        allele_fraction = 0.0
        if depth < min_depth:
            status = "LOW_COVERAGE" if depth > 0 else "UNCALLABLE"
        elif observation is None or not observation.alts:
            status = "REFERENCE"
        elif observation.ref != marker.ref:
            status = "CONFLICT"
        elif marker.alt in observation.alts:
            alt_index = observation.alts.index(marker.alt) + 1
            if alt_index < len(observation.allele_depths):
                expected_alt_depth = observation.allele_depths[alt_index]
            denominator = sum(observation.allele_depths) or depth
            allele_fraction = expected_alt_depth / denominator if denominator else 0.0
            status = "MATCH" if allele_fraction >= min_allele_fraction else "LOW_AF"
        else:
            status = "WRONG_ALT"
        calls[marker_key(marker)] = MarkerCall(
            position=marker.position, observed_ref=observed_ref,
            observed_alt=observed_alt, depth=depth,
            expected_alt_depth=expected_alt_depth, allele_fraction=allele_fraction,
            mapq=None, status=status,
        )
    return calls


def call_markers_from_inputs(
    markers: list[Marker],
    reference: Path,
    inputs: list[Path],
    mode: str,
    threads: int,
    min_depth: int,
    min_allele_fraction: float,
    min_mapq: int,
    min_baseq: int,
) -> tuple[dict[tuple[int, str, str], MarkerCall], MappingStatistics]:
    with tempfile.TemporaryDirectory(prefix="capybara-") as temporary:
        temporary_path = Path(temporary)
        bam = temporary_path / "query.sorted.bam"
        bed = temporary_path / "markers.bed"
        vcf = temporary_path / "markers.vcf"
        map_query(reference, inputs, mode, threads, bam)
        mapped_stats = flagstat(bam)
        stats = MappingStatistics(total=count_input_records(inputs, mode), mapped=mapped_stats.mapped)
        write_bed(bed, fasta_chromosome(reference), {marker.position for marker in markers})
        depths = depth_at_positions(bam, bed, min_mapq, min_baseq)
        call_variants(bam, reference, bed, vcf, min_mapq, min_baseq)
        variants = parse_vcf(vcf)
        calls = evaluate_markers(markers, depths, variants, min_depth, min_allele_fraction)
        return calls, stats
