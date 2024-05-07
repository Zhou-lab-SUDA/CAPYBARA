"""Author: Naclist
Contact: Zxzzzhu@gmail.com
Update time: 2024/5/6
Version: 1.1 Beta"""

from shutil import which
import os

# Define paths to external executable tools
exe = {
    'minimap2': which('minimap2'),
    'mash': which('mash'),
    'samtools': which('samtools'),
    'bcftools': which('bcftools')
}

# Get the directory of the current script
dirname = os.path.dirname(os.path.abspath(__file__))

# Define database folders
db_folder = os.path.join(dirname, 'capydb')
hc_ref = os.path.join(db_folder, 'hc1030_ref.fna')

# Define ESL related paths
esl_folder = os.path.join(db_folder, 'esl')
esl_list = {
    'esl_ref': os.path.join(esl_folder, 'esl_ref.fna'),
    'lineage_snp': os.path.join(esl_folder, 'lineage.SNP'),
    'variant_snp': os.path.join(esl_folder, 'variant.SNP')
}

# Define mash folder and related lists
msh_folder = os.path.join(db_folder, 'msh')
cc_mash_list = os.path.join(msh_folder, 'cc.mash.list')

# Create dictionaries to hold genome data
genome_cc_dict = {}
genomes_msh_list = []

# Populate the genome_cc_dict and genomes_msh_list from cc.mash.list
with open(cc_mash_list) as cc_mash:
    for line in cc_mash:
        part = line.rstrip().split('\t')
        genome, cc = part
        genome_cc_dict[genome] = cc
        genomes_msh_list.append(os.path.join(msh_folder, '.'.join([genome, 'msh'])))

# A series of SNPs for identifying GC1 and GC2

hc1030_snp = {1324677: ['GC1', 'T', 'C'],
              1330326: ['GC1', 'A', 'G'],
              1504581: ['GC1', 'T', 'C'],
              1674628: ['GC1', 'C', 'T'],
              2097573: ['GC1', 'A', 'T'],
              2142755: ['GC1', 'G', 'A'],
              2193024: ['GC1', 'C', 'T'],
              2198442: ['GC1', 'G', 'T'],
              2217841: ['GC1', 'T', 'G'],
              2307186: ['GC1', 'T', 'A'],
              2355160: ['GC1', 'G', 'A'],
              2802506: ['GC1', 'A', 'T'],
              2922385: ['GC1', 'T', 'A'],
              2974472: ['GC1', 'A', 'C'],
              3487156: ['GC1', 'C', 'A'],
              991694: ['GC2', 'T', 'C'],
              1786844: ['GC2', 'T', 'A'],
              1787249: ['GC2', 'C', 'T'],
              2016304: ['GC2', 'G', 'A'],
              2016630: ['GC2', 'G', 'T'],
              2035410: ['GC2', 'G', 'A'],
              2035459: ['GC2', 'T', 'A'],
              }

# A series of SNPs for identifying lineages

lineage_snp = {29867: ['1.1', 'G', 'T'],
               1417459: ['1.2', 'C', 'T'],
               1410148: ['1.3', 'G', 'A'],
               266654: ['1.4', 'A', 'G'],
               947162: ['2.1', 'T', 'A'],
               592890: ['2.2', 'A', 'G'],
               266845: ['2.3', 'G', 'T'],
               2110750: ['2.5', 'T', 'A'],
               }

# A series of SNPs for identifying variants

variant_snp = {161777: ['1.1.1', 'C', 'T'],
               2028825: ['1.2.1', 'T', 'A'],
               2312839: ['1.2.2', 'G', 'T'],
               2223946: ['1.3.10', 'G', 'C'],
               123384: ['1.3.2', 'C', 'T'],
               623806: ['1.3.3', 'C', 'T'],
               807437: ['1.3.4', 'C', 'T'],
               1042347: ['1.3.5', 'A', 'G'],
               210560: ['1.3.6', 'G', 'A'],
               399: ['1.3.7', 'C', 'T'],
               3649: ['1.3.8', 'C', 'A'],
               1759764: ['1.3.9', 'A', 'C'],
               1761833: ['1.4.1', 'C', 'T'],
               417396: ['2.1.1', 'A', 'G'],
               738524: ['2.1.2', 'C', 'T'],
               480007: ['2.1.3', 'C', 'T'],
               2724138: ['2.1.4', 'T', 'G'],
               2048415: ['2.2.1', 'T', 'A'],
               7075: ['2.2.2', 'G', 'A'],
               482779: ['2.2.3', 'C', 'T'],
               1109249: ['2.2.4', 'G', 'T'],
               1671760: ['2.2.5', 'G', 'A'],
               1319866: ['2.2.6', 'C', 'T'],
               277039: ['2.3.1', 'C', 'A'],
               1327555: ['2.3.2', 'C', 'T'],
               2235691: ['2.3.3', 'G', 'A'],
               2113909: ['2.5.1', 'C', 'A'],
               2235965: ['2.5.2', 'G', 'A'],
               629861: ['2.5.3', 'C', 'T'],
               169108: ['2.5.4', 'A', 'C'],
               2246114: ['2.5.5', 'T', 'C'],
               1733723: ['2.5.6', 'T', 'G'],
               1097028: ['2.5.7', 'C', 'T'],
               2860886: ['2.5.8', 'T', 'C'],
               42866: ['2.5.9', 'A', 'C'],
               2906145: ['2.4.10', 'C', 'T'],
               2089945: ['2.4.11', 'C', 'A'],
               299105: ['2.4.12', 'C', 'T'],
               976877: ['2.4.13', 'A', 'G'],
               2172149: ['2.4.14', 'G', 'A'],
               1460331: ['2.4.15', 'T', 'C'],
               1290543: ['2.4.16', 'C', 'T'],
               189843: ['2.4.17', 'G', 'T'],
               1041006: ['2.4.1', 'T', 'C'],
               745014: ['2.4.2', 'G', 'A'],
               351931: ['2.4.4', 'A', 'G'],
               2901430: ['2.4.5', 'G', 'A'],
               454398: ['2.4.6', 'C', 'T'],
               3259086: ['2.4.7', 'T', 'G'],
               1040967: ['2.4.9', 'C', 'A']
               }
