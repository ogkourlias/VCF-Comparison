#!/usr/bin/env python3

"""
    usage:
        ./vcf_prase.py -f [INPUT VCF FILE] -c [VCF COMPARISON FILE]
        -ih [INPUT HEADERS TEXT FILE] -ch [COMPARISON HEADER TEXT FILE]
        -chr [CHROMOSOME TEXT FILE]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import numpy as np
import pandas as pd

# FUNCTIONS
def arg_parse():

    """
    Function to parse arguments.

    :file: Input VCF file.
    :compare: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file named as input chromosome.
    :return: Argparse namespace object.
    """
    
    # Parser init
    parser = argparse.ArgumentParser(
    prog="TabixParser"
    )
    parser.add_argument('-f', '--file', type=str)
    parser.add_argument('-c', '--compare', type=str)
    return parser.parse_args()

def compare(chr, i_handle, c_handle, input_headers_f, comp_headers_f):

    """
    Takes chromosome entry from nextflow process and compare record entries between two files.
    i in var name = VCF File containing all records for which you want to find hits.
    c in var name = vCF File to be compard with.
    :chr: chromosome number in (chr1) format.
    :i_handle: Input VCF file handle.
    :c_handle: Comparison VCF file handle.
    :input_headers_f: Input VCF header file handle.
    :comp_headers_f: Comparison VCF header file handle.
    :return: Pandas dataframes containing records with intersectioned positions.
    """

    # Create pandas DF's with vcf headers as columns.
    i_df = pd.read_csv(input_headers_f, sep='\t')
    c_df = pd.read_csv(comp_headers_f, sep='\t')

    return i_df, c_df

def write_output(i_df, c_df):

    """ 
    Take dataframes with intersections and difference counters to write them to output file.
    
    :i_df: Input dataframe.
    :c_df: Comparison file dataframe.
    :return: None
    """
    file = i_df.iat[0,0]
    print(file)
    i_df.to_csv(f"{file}_1.csv", sep="\t", index=False)
    c_df.to_csv(f"{file}_2.csv", sep="\t", index=False)
    

# MAIN
def main(args):
    """ Main function """
    # Get args.
    args = arg_parse()
    # FINISH
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
