# CAPYBARA

Capybara, a Core-snp Assignment PYthon tool for <i>Acinetobacter baumannii</i>


Capybara enables you to identify hierarchical populations in epidemic super-lineage (ESL) of <i>Acinetobacter baumannii</i> using a set of core-genome SNPs. For ESL or citation of Capybara, see DOI:[[10.21203/rs.3.rs-4129268/v1]]

## Installation:

Capybara was devoloped and tested in Python 3.9.0, and requires a several modules:

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

## Quick Start (with examples)

~~~~~~~~~~
$ cd /path/to/Capybara/

$ Capy.py -i Examples/2.5.6.fna
~~~~~~~~~~

It will generate a report file for Examples/2.5.6.fna about its population.

A single run for an assembled genome will finish <3 minutes for a 4 CPUs laptop (>10 minutes for short reads).

## Usage

~~~~~~~~~~
$ Usage: Capy.py [OPTIONS]

Options:

  -i, --query TEXT  [Required] Input data, both assembled genome or short reads are acceptable.

  -p, --prefix TEXT [Optional] Prefix for output file. Default as Capy.

  -t, --threads INTEGER [Optional] Number of process to use. default: 8

  -l, --list TEXT   [Optional] A file containing list of query files, one per line.

  --help    Show this message and exit.

~~~~~~~~~~

Capybara generates a report file in format below:

| query | ESL | Lineage | Variant |
| ---- | ---- | ---- | ---- |
| 2.5.6.fna | True | 2.5 | 2.5.6 |
| IC7.fna | False | - | - |

## Work flow and Reproduction Instructions

### Work flow

A basic run for Capybara is as follows: 
* (1) ESL identification: 
	* We pre-sketched all 5,824 representative genomes. Genetic distance between query data and pre-sketched data will be evaluated to find the most closed genomes.
	* If query data does not contains any sequential information related to ESL genomes, it will be classified as non-ESL. Otherwise, it will be analyzed as follows.
* (2) Sequential alignment: 
	* Query data will be aligned onto ESL's reference genome (MDR-TJ:GCF_000187205.2) to generate a BAM file.
* (3) SNP calling:
	* A series SNPs will be called from BAM and then generate an VCF file.
* (4) Population assignment:
	* Using a pre-built SNP scheme to assign hierarchical population of query data.


Workflow chart:

<img src="https://github.com/Zhou-lab-SUDA/CAPYBARA/blob/main/workflow.png" width="300px">

### Reproduction Instructions
All data required for reproduction of the analysis were distributed in this repository under xxx/x/xx/

which included:
* esl/esl.fna
~~~
Reference genome for ESL.
~~~
* msh/*.msh
~~~
5,824 pre-sketched files by Mash sketch.
