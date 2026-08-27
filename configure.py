"""Runtime paths and executable configuration only.

Biological barcode definitions live exclusively in capydb/esl/hierarchy.tsv,
markers.tsv, and clade_rules.tsv. Do not add marker dictionaries here.
"""

from pathlib import Path
from shutil import which


PROJECT_ROOT = Path(__file__).resolve().parent
DATABASE_DIR = PROJECT_ROOT / "capydb" / "esl"
REFERENCE = DATABASE_DIR / "esl_ref.fna"

EXECUTABLES = {
    "minimap2": which("minimap2"),
    "samtools": which("samtools"),
    "bcftools": which("bcftools"),
    "mash": which("mash"),
}

DEFAULTS = {
    "threads": 4,
    "min_mapq": 20,
    "min_baseq": 20,
    "assembly_min_depth": 1,
    "assembly_min_af": 0.90,
    "isolate_reads_min_depth": 10,
    "isolate_reads_min_af": 0.90,
    "metagenomic_min_depth": 2,
    "metagenomic_min_af": 0.20,
}
