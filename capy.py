import os
import subprocess
import sys
import tempfile, argparse
from configure import genome_cc_dict, genomes_msh_list, exe, esl_list, lineage_snp, variant_snp, hc1030_snp, hc_ref

mash = exe['mash']
minimap2_exe = exe['minimap2']
samtools_exe = exe['samtools']
bcftools_exe = exe['bcftools']
esl_ref = esl_list['esl_ref']

parser = argparse.ArgumentParser(
    description='Capybara, chipichipi.',
    add_help=True,
    usage='CHIPICHIPI CHAPACHAPA DUBIDUBI DABADABA')
parser.add_argument(
    '-t', '--threads',
    type=str,
    required=False,
    help='''-t/--threads Tree file in newick format.'''
)
parser.add_argument(
    '-o', '--output',
    required=False,
    default='Capy',
    type=str,
    help='-p/--prefix Prefix for output file. Default as Capy.'
)
parser.add_argument(
    '-i', '--input',
    required=True,
    default='Capy',
    type=str,
    help='-i/--input Input data, both assembled genome or short reads are acceptable.'
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

        with open(output_dist, 'r') as dists:
            for i in dists:
                # print(i.rstrip())
                part = i.rstrip().split('\t')
                genome = part[1]
                mash_dist = float(part[2])
                p = float(part[3])
                if mash_dist < closest[2] and p < closest[3]:
                    closest = [genome, genome_cc_dict[genome.replace('.fna.gz', '')], mash_dist, p]
        return closest


def parse_vcf(vcf_file):
    vcf_dict = {}
    with open(vcf_file, 'rt') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                parts = line.rstrip().split('\t')
                # print(parts)
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
            # print(i)
            site = int(i)
            if site in hc1030_snp.keys():
                if hc1030_snp[site][1] == vcf_dict[site][0] and hc1030_snp[site][2] == vcf_dict[site][1]:
                    mapped_hcs.append((hc1030_snp[site][0], vcf_dict[site][2]))
        return mapped_hcs


def parse_seqs():
    query = args.input
    mapped_hcs = find_hc1030_mapping(query)
    GC1_count = 0
    GC2_count = 0
    for i in mapped_hcs:
        if i[0] == 'GC1':
            GC1_count += 1
        elif i[0] == 'GC2':
            GC2_count += 1
    # print('GC1_count', GC1_count/15, 'GC2_count', GC2_count/7)
    if GC1_count > 12 or GC2_count > 4:
        print(query + '\tESL1 detected') if GC1_count > GC2_count else print(query + '\tESL2 detected')
        with open(args.output, 'w') as capout:
            for j in esl_snp(query):
                print('\t'.join([j[0], j[1]]))
                print('\t'.join([query, j[0], j[1]]), file=capout)
        print('Abundance of ESL was saved in ' + args.output)
    else:
        print(query + '\tESL not detected')



parse_seqs()
