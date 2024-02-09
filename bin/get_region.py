#!/usr/bin/env python3

"""
    usage:
        ./vcf_rompare.py -f [INPUT VCF FILE] -r [VCF COMPARISON FILE]
        -ih [INPUT HEADERS TEXT FILE] -ch [COMPARISON HEADER TEXT FILE]
        -chr [CHROMOSOME TEXT FILE] -n [CHUNKSIZE (Variant amount per run)]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import csv
import gzip
import re
from GTFAnnotation import GTFAnnotation
from features.chromosome import Chromosome

# FUNCTIONS
def arg_parse():

    """
    Function to parse arguments.

    :file: Input VCF file.
    :reference: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file named as input chromosome.
    :return: Argparse namespace object.
    """

    # Parser init
    parser = argparse.ArgumentParser(prog="VcfCompare")
    parser.add_argument("-i", "--tsv", type=str)
    parser.add_argument("-g", "--gtf", type=str)
    parser.add_argument("-c", "--chr", type=str)
    parser.add_argument("-n", "--chunk_size", type=int, default=10000)
    parser.add_argument("-o", "--output", type=str)
    return parser.parse_args()

class get_region:
    def __init__(self, tsv, gtf, chr, chunk_size, output):
        self.tsv = tsv
        self.chunk_size = chunk_size
        self.chr = int(chr[3:])
        self.output = output
        self.annotation = GTFAnnotation(gtf)
    
    def get_genes(self):
        self.genes_by_chr = self.annotation.getGenesByChromosome()

    def get_transcripts(self):
        ...

    def get_exons(self):
        ...
    
    def get_regions(self):
        genes = self.annotation.getGenesByChromosome()
        chrom = Chromosome.parse(str(self.chr))
        genes = genes[chrom]
        variant_list = []
        with open(self.output, "w+") as fo:
            fo.write(
            "var_id\t"
            + "is_indel\t"
            + "sample_size\t"
            + "nrhoma_i\t"
            + "nrhets_i\t"
            + "nrhomb_i\t"
            + "nrhoma_r\t"
            + "nrhets_r\t" 
            + "nrhomb_r\t"
            + "matching\t"
            + "match_ratio\t"
            + "maf_i\t"
            + "maf_r\t"
            + "hwe_i\t"
            + "hwe_r\t"
            + "p_corr\t"
            + "s_corr\t"
            + "missing_i\t"
            + "missing_r\t"
            + "gt_total\t"
            + "maf_type_i\t"
            + "maf_type_r\t"
            + "pass_r\t"
            + "gq_i\t"
            + "gq_r\t"
            + "match\t"
            + "gene_types\t"
            + "transcripts\t"
            + "exons\t"
            + "var_status\n"
            )
            with open(self.tsv, "r+") as f:
                spamreader = csv.reader(f, delimiter='\t')
                next(spamreader, None)
                for row in spamreader:
                    variant_list.append(row)
                    if len(variant_list) == self.chunk_size:
                        self.compare(variant_list, genes, fo)
                        variant_list = []
                self.compare(variant_list, genes, fo)
                variant_list = []

    def compare(self, variant_list, genes, fo):
        for variant in variant_list:
            pos = int(variant[0].split(":")[1])
            overlapGenes = []
            region_status = "nan"
            gene_types = []
            overlapTranscripts = []
            overlapExons = []
            for gene in genes:
                if pos >= gene.start and pos <= gene.stop:
                    overlapGenes.append(gene)
            
            if len(overlapGenes) > 0:
                for gene in overlapGenes:
                    gene_types.append(gene.type)
                    transcripts = gene.transcripts
                    for transcript in transcripts:
                        if pos >= transcript.start and pos <= transcript.stop:
                            overlapTranscripts.append(transcript.name)
                            if region_status != "exonic":
                                region_status = "intronic"
                            exons = transcript.exons
                            for exon in exons:
                                if pos >= exon.start and pos <= exon.stop:
                                    overlapExons.append(exon.name)
                                    region_status = "exonic"
            
            overlapTranscripts = ":".join(overlapTranscripts)
            overlapExons = ":".join(overlapExons)
            gene_types = ":".join(gene_types)
                    
            
            variant.append(gene_types)
            variant.append(overlapTranscripts)
            variant.append(overlapExons)
            variant.append(region_status)
            variant_str = "\t".join(variant) + "\n"
            fo.write(variant_str)
    
# MAIN
def main(args):
    """Main function"""
    # Get args.
    args = arg_parse()
    getter = get_region(args.tsv, args.gtf, args.chr, args.chunk_size, args.output)
    getter.get_regions()
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
