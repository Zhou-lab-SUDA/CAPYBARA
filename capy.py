"""Author: Naclist
Contact: Zxzzzhu@gmail.com
Update time: 2024/5/6
version: 1.1 Beta"""

import os
import subprocess
import tempfile, argparse
from configure import genome_cc_dict, genomes_msh_list, exe, esl_list, variant_snp, hc1030_snp, hc_ref

mash = exe['mash']
minimap2_exe = exe['minimap2']
samtools_exe = exe['samtools']
bcftools_exe = exe['bcftools']
esl_ref = esl_list['esl_ref']

parser = argparse.ArgumentParser(
    description='Capybara, a Core-snp Assignment PYthon tool for Acinetobacter baumannii',
    add_help=True,
    usage='Capy.py [OPTIONS]')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', '--query',
                   type=str,
                   help='-i/--input [Required] Input data, both assembled genome or short reads are acceptable.')
group.add_argument('-l', '--list',
                   type=str,
                   help='-l/--list [Optional] A file containing list of query files, one per line.')
parser.add_argument(
    '-t', '--threads',
    type=int,
    default=8,
    required=False,
    help='-t/--threads [Optional] Number of process to use. default: 8'
)
parser.add_argument(
    '-p', '--prefix',
    required=False,
    default='Capy',
    type=str,
    help='-p/--prefix [Optional] Prefix for output file. Default as Capy'
)

args = parser.parse_args()
if not args.threads:
    threads = 8
else:
    threads = args.threads


def find_hc1030(genome):
    with tempfile.TemporaryDirectory() as temp_dir:
        input_msh = os.path.join(temp_dir, "query.msh")
        input_msh_list = os.path.join(temp_dir, "mash.list")
        output_dist = os.path.join(temp_dir, "dist.out")
        with open(input_msh_list, 'w') as iml:
            for i in genomes_msh_list:
                print(i, file=iml)
        # Generate files for mash dist.
        step1 = f"{mash} sketch {genome} -o {input_msh}"
        result = subprocess.run(step1, shell=True, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error executing '{step1}': {result.stderr}")
            return
        step2 = f"{mash} dist {input_msh} -l {input_msh_list} > {output_dist}"
        result = subprocess.run(step2, shell=True, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error executing '{step2}': {result.stderr}")
            return
        closest = ['', '', 1, 0.05]
        # Assume a potential genome with 0.05 genetic distance between query by mash.
        with open(output_dist, 'r') as dists:
            for i in dists:
                part = i.rstrip().split('\t')
                genome = part[1]
                mash_dist = float(part[2])
                p = float(part[3])
                if mash_dist < closest[2] and p < closest[3]:
                    closest = [genome, genome_cc_dict[genome.replace('.fna.gz', '')], mash_dist, p]
        return closest
        # return a list containing query, the most closet genome, mash distance and p-value.


def parse_vcf(vcf_file):
    # parse vcf file for substitution.
    vcf_dict = {}
    with open(vcf_file, 'rt') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                parts = line.rstrip().split('\t')
                pos, ref_base, alt_base, quality = int(parts[1]), parts[3], parts[4], parts[5]
                vcf_dict[pos] = [ref_base, alt_base, quality]
    return vcf_dict


def esl_snp(genome):
    with tempfile.TemporaryDirectory() as temp_dir:
        mp_sam = os.path.join(temp_dir, "mp.sam")
        mp_bam = os.path.join(temp_dir, "mp.bam")
        sorted_bam = os.path.join(temp_dir, "sorted.bam")
        vcf_file = os.path.join(temp_dir, "vcf.file")

        step1 = f"{minimap2_exe} -a {esl_ref} {genome} -o {mp_sam} -t {threads}"
        step2 = f"{samtools_exe} view -bS -o {mp_bam} {mp_sam}"
        step3 = f"{samtools_exe} sort {mp_bam} -o {sorted_bam}"
        step4 = f"{samtools_exe} index {sorted_bam}"
        step5 = f"{bcftools_exe} mpileup -Ou -f {esl_ref} {sorted_bam} | {bcftools_exe} call -mv -Ov -o {vcf_file}"
        # Steps for SNP calling
        for step in [step1, step2, step3, step4, step5]:
            result = subprocess.run(step, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"Error executing command: {step}")
                print(result.stderr)
                return

        vcf_dict = parse_vcf(vcf_file)
        mapped_esl = []
        for i in vcf_dict.keys():
            site = int(i)
            if site in variant_snp.keys():
                mapped_esl.append((variant_snp[site][0], vcf_dict[site][2]))
        return mapped_esl


def find_hc1030_mapping(genome):
    with tempfile.TemporaryDirectory() as temp_dir:
        mp_sam = os.path.join(temp_dir, "mp.sam")
        mp_bam = os.path.join(temp_dir, "mp.bam")
        sorted_bam = os.path.join(temp_dir, "sorted.bam")
        vcf_file = os.path.join(temp_dir, "vcf.file")

        step1 = f"{minimap2_exe} -a {hc_ref} {genome} -o {mp_sam} -t {threads}"
        step2 = f"{samtools_exe} view -bS -o {mp_bam} {mp_sam}"
        step3 = f"{samtools_exe} sort {mp_bam} -o {sorted_bam}"
        step4 = f"{samtools_exe} index {sorted_bam}"
        step5 = f"{bcftools_exe} mpileup -Ou -f {hc_ref} {sorted_bam} | {bcftools_exe} call -mv -Ov -o {vcf_file}"

        for step in [step1, step2, step3, step4, step5]:
            result = subprocess.run(step, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"Error executing command: {step}")
                print(result.stderr)
                return

        vcf_dict = parse_vcf(vcf_file)
        mapped_hcs = []
        for i in vcf_dict.keys():
            site = int(i)
            if site in hc1030_snp.keys():
                if hc1030_snp[site][1] == vcf_dict[site][0] and hc1030_snp[site][2] == vcf_dict[site][1]:
                    mapped_hcs.append((hc1030_snp[site][0], vcf_dict[site][2]))
        return mapped_hcs


def main():
    output_rows = []
    header = ['query', 'ESL', 'Variant', 'Coverage()']
    esl_marker = 0

    if args.query:
        files_to_process = [args.query]
    elif args.list:
        with open(args.list, 'r') as list_file:
            files_to_process = list_file.read().splitlines()

    # Parse file(s)

    for query in files_to_process:
        mapped_hcs = find_hc1030_mapping(query)
        GC1_count = 0
        GC2_count = 0
        for i in mapped_hcs:
            if i[0] == 'GC1':
                GC1_count += 1
            elif i[0] == 'GC2':
                GC2_count += 1

        if GC1_count > 12 or GC2_count > 4:
            esl_detected = True
            esl_marker += 1
        else:
            esl_detected = False
        if esl_detected:
            for j in esl_snp(query):
                [lineage, variant] = j
                output_rows.append([query, str(esl_detected), lineage, variant])
        else:
            [lineage, variant] = ['-', '-']
            output_rows.append([query, str(esl_detected), lineage, variant])

    # Print to file.
    with open(args.prefix + '.tsv', 'w') as capout:
        print('\t'.join(header), file=capout)
        for row in output_rows:
            print('\t'.join(row), file=capout)

    print(f'--Report generated: {args.prefix}.tsv--')
    if esl_marker > 0:
        print('--ESL detected--')
    else:
        print('--No ESL detected--')

if __name__ == '__main__':
    main()
