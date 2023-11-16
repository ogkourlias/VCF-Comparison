#!/usr/bin/env python3

"""
    usage:
        python3 orfeas_gkourlias_deelopdracht01.py
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import tabix
import numpy as np
import pandas as pd
import subprocess
import vcf
import gzip

# FUNCTIONS
def arg_parse():
    """
    Arg parser.

    :file_path: Tabix file
    :return: Argparse namespace object.
    """
    parser = argparse.ArgumentParser(
    prog="TabixParser"
    )
    parser.add_argument('-f', '--file', type=str)
    parser.add_argument('-c', '--compare', type=str)
    return parser.parse_args()

def open_vcf(input_file):
    return vcf.Reader(filename=input_file)

def get_records(vcf_reader):

    """
    Handle tabix file and parse lines to line_parser.

    :file_path: Tabix file
    :return: None
    """
    chr_dict = {}
    for record in vcf_reader:
        if record.CHROM in chr_dict:
            chr_dict[record.CHROM].append(record)
        else:
            chr_dict[record.CHROM] = [record]

    return chr_dict

def isec_for_chr(chr_dict, comp_f):
    tbi_f = tabix.open(comp_f)
    rec_list = []
    rec_list_c = []
    for chr in chr_dict.values():
        selection = get_selection_chrom(chr, tbi_f)
        selection = [rec for rec in selection]
        for rec in chr:
            print(rec.samples[0])
            rec = [rec.CHROM, rec.POS, rec.ID, rec.REF, rec.ALT, rec.QUAL, rec.FILTER, rec.INFO, rec.FORMAT, rec.samples]
            for rec_c in selection:
                #print(rec_c[9:])
                #print(" ".join(rec_c[6:]))
                if rec[1] == int(rec_c[1]):
                    rec_list.append(rec)
                    rec_list_c.append(rec_c)

        #chr_isec_dict[chr[0].CHROM] = [rec for rec in chr if rec.POS in [int(rec_c[1]) for rec_c in selection]]
    
    return rec_list, rec_list_c   

def get_selection_chrom(chr_entry, tb_f):
    chr = chr_entry[0].CHROM
    min, max = chr_entry[0].POS, chr_entry[-1].POS
    return tb_f.query(chr, min, max)

def make_df(isecs):
    df = pd.DataFrame(columns=["1","2"])
    for chr in isecs.values():
        for file in chr:
            ...
    

def line_parser(self, records, gene):

    """
    Parses the line from tabix and gathers all the information for the array.

    :param records: Tabix line
    :param gene: Corresponding gene
    :return: None
    """

    for element in records:
        rs_id = element[2].split(":")[2]

        # Check if combination exists in the test dataframe
        combination_exists = self.test_df[(self.test_df['Gene'] == gene) & (self.test_df['SNP'] == rs_id)]

        if combination_exists.empty:
            continue

        a1, a2 = element[3], element[4]

        # Extract genotypes
        genotypes = [genome.split(":")[0].count("1") if not genome.startswith("./") else np.nan for genome in element[9:]]

        # Extract expression data
        expression = []
        for name in self.names_vcf:
            try:
                expression_name = self.GTE_dict[name]
                expression.append(self.expression.loc[gene,expression_name])
            except KeyError:
                expression.append(np.nan)

        # Extract disease status
        disease_status = [self.diseaselabels_dict.get(name, np.nan) for name in self.names_vcf]

        # Create data array
        data_array = np.array(list(zip(expression, genotypes, disease_status)))

        # Filter out rows with NaN values
        masked_array = data_array[~np.isnan(data_array).any(axis=1)]

        if masked_array.shape[0] == 0:
            continue

# MAIN
def main(args):
    """ Main function """
    args = arg_parse()
    input_f = open_vcf(args.file)
    records = get_records(input_f)
    rec_list, rec_list_c = isec_for_chr(records, args.compare)
    #make_df(rec_list, rec_list_c)
    #min_pos, max_pos = get_range(vcf_reader)
    #selected_comp = get_selection(min_pos, max_pos, vcf_reader_comp)
    #selected = get_selection(min_pos, max_pos, vcf_reader)
    #for x, y in zip(selected, selected_comp):
       #print(x,y)
    # FINISH
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
