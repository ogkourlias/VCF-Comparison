#!/usr/bin/env python3

"""
    usage:
        ./vcf_compare.py -f [INPUT VCF FILE] -c [VCF COMPARISON FILE]
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
import tabix
import gzip
import numpy as np
import pandas as pd
from itertools import zip_longest

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
    parser.add_argument('-ih', '--input_header', type=str)
    parser.add_argument('-ch', '--comp_header', type=str)
    parser.add_argument('-chr', '--chr', type=str)
    parser.add_argument('-o', '--output', type=str)
    return parser.parse_args()

def open_tbi(input_file):
    
    """
    Tabix file opener function.
    :input_file: Input VCF file.
    :return: Opened tabix file.
    """
    return tabix.open(input_file)

def compare2(chr, i_handle, c_handle, input_headers_f, comp_headers_f):

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

    # out_f = f"{chr}.txt"
    
    # Get generator for input vcf file.
    rec_gen = i_handle.query(chr, 1, 999999999)
    # Iterate over generator to create list and find min/max positions of records.
    recs = [rec for rec in rec_gen]
    min, max = recs[0][1], recs[-1][1]

    # Query on chromosome in range of min, max values.
    c_rec_gen = c_handle.query(chr, int(min), int(max))
    recs_c = [rec for rec in c_rec_gen]

    # Initiate alt difference counter.
    alt_same = []

    # For each record in the input file.
    for rec in recs:
        # For each record within selected region of comparison file.
        for rec_c in recs_c:
            # If overlapping positions are found between records.
            if rec[1] == rec_c[1]:
                # Insert row at end of dataframes.
                i_df.loc[len(i_df)] = rec
                c_df.loc[len(c_df)] = rec_c
                
                # Check if SNP's are equal.
                # If not, append 1 for counted difference.
                if rec[4] == rec_c[4]:
                    alt_same.append(1)
                else:
                    alt_same.append(0)
    
    # Insert new column with the alt difference marker.
    i_df['alt_same'] = alt_same
    c_df['alt_same'] = alt_same

    return i_df, c_df

def compare(chr, i_handle, c_handle, headers, output, chunk_size = 100):

    """
    Takes chromosome entry from nextflow process and compare record entries between two files.
    i in var name = VCF File containing all records for which you want to find hits.
    c in var name = vCF File to be compard with.
    :chr: chromosome number in (chr1) format.
    :i_handle: Input VCF file handle.
    :c_handle: Comparison VCF file handle.
    :input_headers_f: Input VCF header file handle.
    :comp_headers_f: Comparison VCF header file handle.
    :re
    turn: Pandas dataframes containing records with intersectioned positions.
    """

    # Get start and stop positions of chr selection.
    start, stop = get_range(chr, i_handle)

    # Open and query tabix comparison file.
    comp_f = open_tbi(c_handle)
    comp_i = open_tbi(i_handle)

    #comp_gen = comp_f.query(chr, int(start), int(stop))start

    for chunk in range(start, stop, chunk_size):
        i_records = comp_i.query(chr, chunk, chunk+chunk_size)
        i_records = [rec for rec in i_records]

        c_records = comp_f.query(chr, chunk, chunk+chunk_size)

        with open(output, 'w') as f:
            f.write("\t".join(headers.keys()) + "\n")
            for c_record in c_records:
                for i_record in i_records:
                    if i_record[1] == c_record[1] and i_record[3] == c_record[3] and i_record[4] == c_record[4]:
                        record_extract(i_record, c_record)
                        
                    
        

        # Comparison



    return 0

def header_subset(input_headers_f, comp_headers_f):
    head_idx = {}
    for header_i in np.genfromtxt(input_headers_f, delimiter='\t', dtype=str):
        for i, header_c in enumerate(np.genfromtxt(comp_headers_f, delimiter='\t', dtype=str)):
            if header_i == header_c:
                head_idx[header_i] = i
                break

    return head_idx

def get_range(chr, i_handle):
    # Get start and stop positions of chr selection.
    with gzip.open(i_handle, 'rt') as f:
        i = 0
        for line in f:
            if line.startswith(chr):
                if i == 0:
                    start = int(line.split("\t")[1])
                    i =+ 1

                else:
                    stop = int(line.split("\t")[1])
    return start, stop

def record_extract(c_record, i_record):
    c_record[0\]
    return entry

def write_output(i_df, c_df):

    """
    Take dataframes with intersections and difference counters to write them to output file.
    
    :i_df: Input dataframe.
    :c_df: Comparison file dataframe.
    :return: None
    """
    file = i_df.iat[0,0]
    i_df.to_csv(f"{file}_1.csv", sep="\t", index=False)
    c_df.to_csv(f"{file}_2.csv", sep="\t", index=False)
    

# MAIN
def main(args):
    """ Main function """
    # Get args.
    args = arg_parse()
    # Get dataframes with compared and intersecting entries.
    headers = header_subset(args.input_header, args.comp_header)
    compare(args.chr, args.file, args.compare, headers, args.output)
    # Write to tsv file for next nextflow process.
    #write_output(i_df, c_df)
    # FINISH
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
