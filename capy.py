"""Author: Naclist
Contact: Zxzzzhu@gmail.com
Update time: 2026/3/25
Version: 1.1.4

This file is part of Capybara, a Core-snp Assignment PYthon tool for Acinetobacter baumannii, which is totally free to use and redevelop based on the GNU General Public License. Capybara is designed for pathogen detection and is not for profit. For more detailed information on the usage of the program, please refer to the original GNU license text: http://www.gnu.org/licenses/

This version is modified to fit the debug from Lean.

"""

#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, os, subprocess, sys, tempfile
from collections import defaultdict
from typing import Dict, Iterable, Tuple
from configure import exe, esl_list, lineage_snp, variant_snp

try:
    from configure import hc1030_snp, hc_ref
    HAS_HC = True
except Exception:
    hc1030_snp = {}
    hc_ref = None
    HAS_HC = False

MINIMAP2 = exe.get("minimap2", "minimap2")
SAMTOOLS = exe.get("samtools", "samtools")
BCFTOOLS = exe.get("bcftools", "bcftools")
ESL_REF = esl_list.get("esl_ref")


############################################################
# Utils
############################################################

def run(cmd: str, check: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True)
    if check and proc.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\nSTDERR:\n{proc.stderr}")
    return proc


def first_fasta_header(fasta_path: str) -> str:
    with open(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].strip().split()[0]
    raise ValueError(f"No FASTA header in {fasta_path}")


def write_bed_for_positions(bed_path: str, chrom: str,
                            positions: Iterable[int]) -> None:
    with open(bed_path, "wt") as out:
        for p in sorted(set(int(x) for x in positions)):
            out.write(f"{chrom}\t{p-1}\t{p}\n")


def depth_for_positions(bam: str, ref: str, bed_path: str) -> Dict[int, int]:
    with tempfile.TemporaryDirectory() as td:
        out_tsv = os.path.join(td, "dp.tsv")
        cmd = (
            f"{BCFTOOLS} mpileup -a DP -f {ref} -R {bed_path} {bam} -Ou | "
            f"{BCFTOOLS} query -f '%CHROM\\t%POS\\t%DP\\n' > {out_tsv}"
        )
        run(cmd)
        dp = {}
        with open(out_tsv, "rt") as fh:
            for line in fh:
                _, pos, dp_str = line.rstrip().split("\t")
                try:
                    dp[int(pos)] = int(dp_str)
                except:
                    dp[int(pos)] = 0
        return dp


def targeted_call_with_af(bam: str, ref: str, bed: str,
                          out_vcf: str, threads: int) -> None:
    cmd = (
        f"{BCFTOOLS} mpileup -Ou -f {ref} -R {bed} "
        f"-a FORMAT/AD,FORMAT/DP {bam} | "
        f"{BCFTOOLS} call -mv -Ov -o {out_vcf}"
    )
    run(cmd)


def parse_vcf_with_af(vcf_path: str) -> Dict[int, Tuple[str, str, int, int, float]]:
    calls = {}
    with open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            pos = int(parts[1])
            ref, alt = parts[3], parts[4]
            alt0 = alt.split(",")[0]
            fmt = parts[8].split(":") if len(parts) > 8 else []
            smp = parts[9].split(":") if len(parts) > 9 else []
            ad_ref = ad_alt = dp = 0

            if fmt and smp:
                idx = {k: i for i, k in enumerate(fmt)}
                if "AD" in idx and idx["AD"] < len(smp):
                    try:
                        ad = smp[idx["AD"]].split(",")
                        ad_ref = int(ad[0])
                        ad_alt = int(ad[1]) if len(ad) > 1 else 0
                    except:
                        ad_ref, ad_alt = 0, 0
                if "DP" in idx and idx["DP"] < len(smp):
                    try:
                        dp = int(smp[idx["DP"]])
                    except:
                        dp = ad_ref + ad_alt

            denom = max(1, ad_ref + ad_alt)
            af = ad_alt / denom
            calls[pos] = (ref, alt0, ad_ref, ad_alt, af)
    return calls


def map_to_bam(query: str, ref: str, threads: int,
               subsample: float, out_bam: str) -> None:
    subs = f"-s {subsample}" if subsample < 1.0 else ""
    cmd = (
        f"{MINIMAP2} -a {ref} {query} -t {threads} | "
        f"{SAMTOOLS} view -bS {subs} - | "
        f"{SAMTOOLS} sort -@ {threads} -o {out_bam} -"
    )
    run(cmd)
    run(f"{SAMTOOLS} index {out_bam}")


def samtools_flagstat(bam: str) -> Dict[str, int]:
    out = run(f"{SAMTOOLS} flagstat {bam}").stdout
    total = mapped = 0
    for line in out.splitlines():
        if " in total " in line:
            total = int(line.split(" ")[0])
        if " mapped (" in line and " primary " not in line:
            mapped = int(line.split(" ")[0])
    return {"total": total, "mapped": mapped}


############################################################
# Panel evaluation
############################################################

def evaluate_panel(
    bam: str,
    ref: str,
    panel: Dict[int, Tuple[str, str, str]],
    min_dp: int,
    min_af: float,
    is_meta: bool
) -> Dict[int, Tuple[str, str, str, int, float, bool]]:

    chrom = first_fasta_header(ref)
    with tempfile.TemporaryDirectory() as td:
        bed = os.path.join(td, "panel.bed")
        write_bed_for_positions(bed, chrom, panel.keys())
        dp_map = depth_for_positions(bam, ref, bed)
        vcf = os.path.join(td, "panel.vcf")
        targeted_call_with_af(bam, ref, bed, vcf, threads=1)
        vcf_calls = parse_vcf_with_af(vcf)

    enriched = {}
    for pos, (label, exp_ref, exp_alt) in panel.items():
        ref_obs, alt_obs, ad_ref, ad_alt, af = vcf_calls.get(
            pos, (".", "", 0, 0, 0.0)
        )
        dp = dp_map.get(pos, ad_ref + ad_alt)
        if is_meta:
            pass_flag = (dp >= min_dp and af >= min_af)
        else:
            pass_flag = (dp >= 1 and af >= min_af)
        enriched[pos] = (label, exp_ref, exp_alt, dp, af, pass_flag)
    return enriched


def evaluate_hc1030(
    bam: str,
    ref: str,
    panel: Dict[int, Tuple[str, str, str]],
    min_dp: int,
    min_af: float,
    is_meta: bool
):

    chrom = first_fasta_header(ref)
    with tempfile.TemporaryDirectory() as td:
        bed = os.path.join(td, "hc1030.bed")
        write_bed_for_positions(bed, chrom, panel.keys())
        dp_map = depth_for_positions(bam, ref, bed)
        vcf = os.path.join(td, "hc1030.vcf")
        targeted_call_with_af(bam, ref, bed, vcf, threads=1)
        vcf_calls = parse_vcf_with_af(vcf)

    gc1_votes = 0
    gc2_votes = 0
    hc_calls = {}

    depth_threshold = min_dp if is_meta else 1

    for pos, (label, exp_ref, exp_alt) in panel.items():
        ref_obs, alt_obs, ad_ref, ad_alt, af = vcf_calls.get(
            pos, (".", "", 0, 0, 0.0)
        )
        dp = dp_map.get(pos, ad_ref + ad_alt)
        status = "lowcov"

        if dp >= depth_threshold and af >= min_af:

            if ref_obs == exp_ref and alt_obs == exp_alt:
                status = "match"
                if label == "GC1":
                    gc1_votes += 1
                elif label == "GC2":
                    gc2_votes += 1
            else:
                status = "mismatch"

        hc_calls[pos] = (label, exp_ref, exp_alt, dp, af, status)

    return hc_calls, gc1_votes, gc2_votes


############################################################
# Scoring
############################################################

def score_variants(
    variant_calls: Dict[int, Tuple[str, str, str, int, float, bool]]
) -> Dict[str, float]:
    scores = defaultdict(float)
    for _pos, (label, _r, _a, dp, af, ok) in variant_calls.items():
        if ok:
            scores[label] += dp * af
    return scores


def filter_scores(scores: Dict[str, float], cutoff: float) -> Dict[str, float]:
    if not scores:
        return {}
    max_score = max(scores.values())
    threshold = max_score * cutoff
    return {k: v for k, v in scores.items() if v >= threshold and v > 0}


############################################################
# Main
############################################################

def main():
    p = argparse.ArgumentParser(
        description="Capybara lineage/sublineage caller with hc1030 gating"
    )
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--output", default="Capy")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("-m", "--metagenomic", action="store_true")
    p.add_argument("--subsample", type=float, default=1.0)
    p.add_argument("--min-snp-depth", type=int, default=None)
    p.add_argument("--min-allele-frac", type=float, default=None)
    p.add_argument("--min-score-frac", type=float, default=0.2,
                   help="minimum fraction of max total score to retain (default=0.2)")
    args = p.parse_args()

    if not os.path.exists(ESL_REF):
        sys.exit(f"ESL reference not found: {ESL_REF}")

    if args.metagenomic:
        if args.min_snp_depth is None:
            args.min_snp_depth = 2
        if args.min_allele_frac is None:
            args.min_allele_frac = 0.8
    else:
        if args.min_snp_depth is None:
            args.min_snp_depth = 1
        if args.min_allele_frac is None:
            args.min_allele_frac = 0.9

    sample = os.path.basename(args.input)
    mode = "metagenomic" if args.metagenomic else "isolate"

    gc1_votes = 0
    gc2_votes = 0
    hc_trusted = True
    hc_calls = {}

    with tempfile.TemporaryDirectory() as td:
        # ESL mapping
        bam_esl = os.path.join(td, "esl.sorted.bam")
        map_to_bam(args.input, ESL_REF, args.threads, args.subsample, bam_esl)
        flag_esl = samtools_flagstat(bam_esl)

        # ESL panels
        lineage_calls = evaluate_panel(
            bam_esl, ESL_REF, lineage_snp,
            args.min_snp_depth, args.min_allele_frac, args.metagenomic
        )
        variant_calls = evaluate_panel(
            bam_esl, ESL_REF, variant_snp,
            args.min_snp_depth, args.min_allele_frac, args.metagenomic
        )

        # Scores
        variant_scores = score_variants(variant_calls)
        lineage_scores = score_variants(lineage_calls)

        filtered_variants = filter_scores(
            variant_scores, args.min_score_frac
        )
        filtered_lineages = {}
        if not filtered_variants:
            filtered_lineages = filter_scores(
                lineage_scores, args.min_score_frac
            )

        if (not args.metagenomic) and HAS_HC and hc_ref:
            bam_hc = os.path.join(td, "hc1030.sorted.bam")
            map_to_bam(args.input, hc_ref, args.threads,
                       args.subsample, bam_hc)
            hc_calls, gc1_votes, gc2_votes = evaluate_hc1030(
                bam_hc, hc_ref, hc1030_snp,
                args.min_snp_depth, args.min_allele_frac, args.metagenomic
            )

            if gc1_votes == 0 and gc2_votes == 0:
                hc_trusted = False
            elif gc1_votes > 0 and gc2_votes > 0:

                hc_trusted = False
            else:
                hc_trusted = True

            if not hc_trusted:
                filtered_variants = {}
                filtered_lineages = {}

        with open(f"{args.output}.summary.tsv", "wt") as outs:
            outs.write(
                "sample\tmode\tlineage\tsublineage\tGC1_votes\tGC2_votes"
                "\treads_total\treads_mapped\tbest_score\tvariant_total\n"
            )

            if filtered_variants:
                for subl, score in sorted(
                    filtered_variants.items(),
                    key=lambda x: x[1],
                    reverse=True
                ):
                    lineage = (
                        ".".join(subl.split(".")[:2])
                        if "." in subl else subl
                    )
                    outs.write(
                        f"{sample}\t{mode}\t{lineage}\t{subl}\t"
                        f"{gc1_votes}\t{gc2_votes}\t"
                        f"{flag_esl['total']}\t{flag_esl['mapped']}\t"
                        f"{round(score,3)}\t{len(variant_snp)}\n"
                    )

            elif filtered_lineages:
                for lin, score in sorted(
                    filtered_lineages.items(),
                    key=lambda x: x[1],
                    reverse=True
                ):
                    outs.write(
                        f"{sample}\t{mode}\t{lin}\tNA\t"
                        f"{gc1_votes}\t{gc2_votes}\t"
                        f"{flag_esl['total']}\t{flag_esl['mapped']}\t"
                        f"{round(score,3)}\t{len(lineage_snp)}\n"
                    )
            else:
                outs.write(
                    f"{sample}\t{mode}\tNA\tNA\t"
                    f"{gc1_votes}\t{gc2_votes}\t"
                    f"{flag_esl['total']}\t{flag_esl['mapped']}\t"
                    f"0\t0\n"
                )

        with open(f"{args.output}.markers.tsv", "wt") as outm:
            outm.write("panel\tlabel\tpos\tref\talt\tdp\taf\tstatus\n")
            for pos, (label, exp_ref, exp_alt, dp, af, _ok) in lineage_calls.items():
                outm.write(
                    f"lineage\t{label}\t{pos}\t"
                    f"{exp_ref}\t{exp_alt}\t{dp}\t{af:.3f}\t.\n"
                )
            for pos, (label, exp_ref, exp_alt, dp, af, _ok) in variant_calls.items():
                outm.write(
                    f"variant\t{label}\t{pos}\t"
                    f"{exp_ref}\t{exp_alt}\t{dp}\t{af:.3f}\t.\n"
                )
            if hc_calls:
                for pos, (label, exp_ref, exp_alt, dp, af, status) in hc_calls.items():
                    outm.write(
                        f"hc1030\t{label}\t{pos}\t"
                        f"{exp_ref}\t{exp_alt}\t{dp}\t{af:.3f}\t{status}\n"
                    )

        # coverage.json
        cov = {
            "sample": sample,
            "mode": mode,
            "filtered_variants": filtered_variants,
            "filtered_lineages": filtered_lineages,
            "variant_scores": variant_scores,
            "lineage_scores": lineage_scores,
            "reads_total": flag_esl.get("total", 0),
            "reads_mapped": flag_esl.get("mapped", 0),
            "gc1_votes": gc1_votes,
            "gc2_votes": gc2_votes,
            "hc_trusted": hc_trusted,
        }

        with open(f"{args.output}.coverage.json", "wt") as outj:
            json.dump(cov, outj, indent=2)

    print(f"[OK] Wrote: {args.output}.summary.tsv, {args.output}.markers.tsv, {args.output}.coverage.json")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[ERROR] {e}\n")
        sys.exit(1)

