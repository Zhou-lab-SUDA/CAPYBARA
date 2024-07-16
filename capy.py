"""Author: Naclist
Contact: Zxzzzhu@gmail.com
Update time: 2024/7/14
Version: 1.1.2 Beta

This file is part of Capybara, a Core-snp Assignment PYthon tool for Acinetobacter baumannii, which is totally free to use and redevelop based on the GNU General Public License. Capybara is designed for pathogen detection and is not for profit. For more detailed information on the usage of the program, please refer to the original GNU license text: http://www.gnu.org/licenses/


"""

import os
import subprocess
import tempfile
import argparse
from configure import genome_cc_dict, genomes_msh_list, exe, esl_list, genotype_snp, hc1030_snp, hc_ref

# External executable paths
mash = exe['mash']
minimap2_exe = exe['minimap2']
samtools_exe = exe['samtools']
bcftools_exe = exe['bcftools']
esl_ref = esl_list['esl_ref']

# Argument parser setup
parser = argparse.ArgumentParser(
    description='Capybara, a Core-snp Assignment PYthon tool for Acinetobacter baumannii',
    add_help=True,
    usage='Capy.py [OPTIONS]')

# Mutually exclusive group for input arguments
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', '--query', type=str, help='-i/--query [Required] Input data, both assembled genome or short reads are acceptable.')
group.add_argument('-l', '--list', type=str, help='-l/--list [Optional] A file containing list of query files, one per line.')

# Additional optional arguments
parser.add_argument('-t', '--threads', type=int, default=8, help='-t/--threads [Optional] Number of processes to use. Default: 8')
parser.add_argument('-p', '--prefix', required=False, default='Capy', type=str, help='-p/--prefix [Optional] Prefix for output file. Default as Capy')

args = parser.parse_args()
threads = args.threads if args.threads else 8

# Function to find high-confidence genomes based on mash distance
def find_hc1030(genome):
    with tempfile.TemporaryDirectory() as temp_dir:
        # Temporary files setup
        input_msh = os.path.join(temp_dir, "query.msh")
        input_msh_list = os.path.join(temp_dir, "mash.list")
        output_dist = os.path.join(temp_dir, "dist.out")

        # Writing the list of genome mash files
        with open(input_msh_list, 'w') as iml:
            for i in genomes_msh_list:
                print(i, file=iml)

        # Step 1: Sketching the query genome with mash
        step1 = f"{mash} sketch {genome} -o {input_msh}"
        result = subprocess.run(step1, shell=True, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error executing '{step1}': {result.stderr}")
            return

        # Step 2: Calculating distance between sketched genomes
        step2 = f"{mash} dist {input_msh} -l {input_msh_list} > {output_dist}"
        result = subprocess.run(step2, shell=True, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error executing '{step2}': {result.stderr}")
            return

        # Finding the closest genome
        closest = ['', '', 1, 0.05]  # Default values for comparison
        with open(output_dist, 'r') as dists:
            for i in dists:
                part = i.rstrip().split('\t')
                genome, mash_dist, p = part[1], float(part[2]), float(part[3])
                if mash_dist < closest[2] and p < closest[3]:
                    closest = [genome, genome_cc_dict[genome.replace('.fna.gz', '')], mash_dist, p]
        return closest

# Function to parse VCF files for identified SNPs
def parse_vcf(vcf_file):
    vcf_dict = {}
    with open(vcf_file, 'rt') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                parts = line.rstrip().split('\t')
                pos, ref_base, alt_base, quality = int(parts[1]), parts[3], parts[4], parts[5]
                vcf_dict[pos] = [ref_base, alt_base, quality]
    return vcf_dict

# Function to perform SNP calling for ESL markers
def esl_snp(genome):
    with tempfile.TemporaryDirectory() as temp_dir:
        # Setup temporary files for sequence alignment and SNP calling
        mp_sam = os.path.join(temp_dir, "mp.sam")
        mp_bam = os.path.join(temp_dir, "mp.bam")
        sorted_bam = os.path.join(temp_dir, "sorted.bam")
        vcf_file = os.path.join(temp_dir, "vcf.file")

        # Sequence alignment and conversion to binary format
        step1 = f"{minimap2_exe} -a {esl_ref} {genome} -o {mp_sam} -t {threads}"
        step2 = f"{samtools_exe} view -bS -o {mp_bam} {mp_sam}"
        step3 = f"{samtools_exe} sort {mp_bam} -o {sorted_bam}"
        step4 = f"{samtools_exe} index {sorted_bam}"
        step5 = f"{bcftools_exe} mpileup -Ou -f {esl_ref} {sorted_bam} | {bcftools_exe} call -mv -Ov -o {vcf_file}"

        # Execute each step and check for errors
        for step in [step1, step2, step3, step4, step5]:
            result = subprocess.run(step, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"Error executing command: {step}")
                print(result.stderr)
                return

        # Analyze VCF for SNP matching
        vcf_dict = parse_vcf(vcf_file)
        mapped_esl = []
        for site in vcf_dict:
            if site in genotype_snp:
                mapped_esl.append((genotype_snp[site][0], vcf_dict[site][2]))
        return mapped_esl

# Main execution flow
def main():
    output_rows = []
    header = ['query', 'ESL', 'Clade', 'Coverage']
    esl_marker = 0

    # Determine the source of input files
    files_to_process = [args.query] if args.query else open(args.list, 'r').read().splitlines()

    # Process each file for SNP mapping and variant detection
    for query in files_to_process:
        mapped_hcs = find_hc1030(query)
        print(mapped_hcs)
        
        GC1_count = 1 if mapped_hcs[1] == 'Clonal_complex_of_ST1' else 0
        GC2_count = 1 if mapped_hcs[1] == 'Clonal_complex_of_ST2' else 0

        esl_detected = GC1_count > 0 or GC2_count > 0
        esl_marker += esl_detected

        # Collect results for detected ESL markers
        if esl_detected:
            for j in esl_snp(query):
                lineage, variant = j
                output_rows.append([query, str(esl_detected), lineage, variant])
        else:
            output_rows.append([query, str(esl_detected), '-', '-'])

    # Write results to a TSV file
    with open(f'{args.prefix}.tsv', 'w') as capout:
        print('\t'.join(header), file=capout)
        for row in output_rows:
            print('\t'.join(row), file=capout)

    print(f'--Report generated: {args.prefix}.tsv--')
    if esl_marker:
        print('--ESL detected--')
    else:
        print('--No ESL detected--')

if __name__ == '__main__':
    main()
