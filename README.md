# VCF Stats

VCF Stats is a collection of VCF parser scripts used to sumamrise and explore variant call format (VCF) files.
This tool consists of standalone python, R scripts and a nextflow pipeline framework. The scripts are able to do the following:

- Calculate variant specific statistics such as MAF, HWE, call rate and allele frequncy. 
- Recount the amount of alternative and reference calls for each genotype.
- Find variant locations using a reference GTF.
- Create visualisations showing the allele frequencies and call rates.

The statistics are either written in a compressed TSV format or PNG format in the case of the R script.

## Installation
To install and use the scripts, simply make sure you have the Docker container built.
You can either build it from the Dockerfile or clone it from the following repository on dockerhub like this:

```bash
docker pull ogkourlias/vcf_compare:latest
```
Once the docker container is installed and running, you are ready to use the individual scripts.

## Usage
There are three distinct scripts you are able to call upon for the primary VCF analysis.

1. vcf_stats.py
2. get_region.py
3.  summary.R

### vcf_Stats.py
This script takes a VCF.GZ input file and creates and output file containg summary statistics for each variant.
To run the script, some arguments need to be passed first. Fill in the blanks below and you can run the script.

```bash
./vcf_stats.py -i <path_to_input_vcf_gz> -o <path_to_output_tsv.gz>
```

### get_regions.py
The get_regions script takes the tsv.gz from the vcf_stats script append the new variant regions to a new .tsv.gz file.
It is therefore requried that the vcf_stats script has been run prior.
You are required to have a GTF reference file. The -all argument option is required when you are running this is a standalone script.

```bash
./get_regions.py -all -i <path_to_input_tsv_gz> -o <path_to_output_tsv.gz> -g <path_to_reference_gtf_file.gz>
```
### summary.R
This R script takes the prior tsv.gz file from get_regions.py to create PNG images containing plots of allele frequencie and call rates.
It creates a new directory in the working dirctory called ''figures'' with all the generated images.

```bash
Rscript summary.R <path_to_input_tsv_gz>
```

## License

[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)