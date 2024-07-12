# CAPYBARA

Capybara is a Core-snp Assignment PYthon tool for <i>Acinetobacter baumannii</i>. It screens either raw reads or assemblies for :

* Identifying whether a query belongs to the epidemic super-lineage (ESL), a super-set of two predominant international clones: IC1 and IC2.
* Assignment of query strain into one of the lineages, clusters, and clades in the ESL based on a pre-curated set of SNPs. 

## Citation
Shengkai Li, Heng Li, Guilai Jiang, Shengke Wang, Min Wang, Yilei Wu, Xiao Liu, Ling Zhong, Shichang Xie, Yi Ren, Yongliang Lou, Jimei Du, Zhemin Zhou, 2024, Emergence and Global Spread of a Dominant Multidrug-Resistant Variant in *Acinetobacter baumannii*, https://doi.org/10.21203/rs.3.rs-4224555/v1

------

## Installation:
Capybara was devoloped and tested in Python 3.9.0, and requires several modules:

~~~~~~~~~~
minimap2
mash
samtools
bcftools
~~~~~~~~~~

You can easily install these packages using command below:

~~~~~~~~~~
conda install -c bio-conda samtools bcftools minimap2 mash
~~~~~~~~~~

Then you can use git to clone Capybara into your PC.

~~~~~~~~~~
git clone git@github.com:Zhou-lab-SUDA/CAPYBARA.git
~~~~~~~~~~

## Quick Start (with examples)

~~~~~~~~~~
$ cd /path/to/Capybara/

$ capy.py -i Examples/2.5.6.fna
~~~~~~~~~~

It will generate a report file for Examples/2.5.6.fna about its population.

A single run for an assembled genome will finish <3 minutes for a 4 CPUs laptop (>10 minutes for short reads).

## Usage

~~~~~~~~~~
$ Usage: capy.py [OPTIONS]

Options:

  -i, --query TEXT  [Required] Input data, both assembled genome or short reads are acceptable.

  -p, --prefix TEXT [Optional] Prefix for output file. Default as Capy.

  -t, --threads INTEGER [Optional] Number of process to use. default: 8

  -l, --list TEXT   [Optional] A file containing list of query files, one per line.

  --help    Show this message and exit.

~~~~~~~~~~

Capybara generates a report file in format below:

| query | ESL | Clusters | Clades |
| ---- | ---- | ---- | ---- |
| 2.5.6.fna | True | 2.5 | 2.5.6 |
| IC7.fna | False | - | - |

## Work flow and Reproduction Instructions

### Work flow

A basic run for Capybara is as follows: 
* **ESL identification:**
	* We pre-sketched all 5,824 representative genomes. Genetic distance between query data and pre-sketched data will be evaluated to find the most closed genomes.
	* If query data does not contains any sequential information related to ESL genomes, it will be classified as non-ESL. Otherwise, it will be analyzed as follows.
* **Sequential alignment:**
	* Query data will be aligned onto ESL's reference genome (MDR-TJ:GCF_000187205.2) to generate a BAM file.
* **SNP calling:**
	* A series SNPs will be called from BAM and then generate an VCF file.
* **Population assignment:**
	* Using a pre-built SNP scheme to assign hierarchical population of query data.


Workflow chart:

<img src="https://github.com/Zhou-lab-SUDA/CAPYBARA/blob/main/workflow.png" width="300px">

### Reproduction Instructions
All data required for reproduction of the analysis were distributed in this repository under CAPYBARA/capydb/

which included:
* esl/esl.fna
~~~
Reference genome for ESL.
~~~
* msh/*.msh
~~~
5,824 pre-sketched files by Mash sketch.
~~~
### Our pulished release
You may also be interested [KleTy](https://github.com/Zhou-lab-SUDA/KleTy), a pipeline for analysis of Klebsiella and can also genotype plasmids from short-reads file.


<a href="https://zenodo.org/doi/10.5281/zenodo.11349698"><img src="https://zenodo.org/badge/803854575.svg" alt="DOI"></a>
