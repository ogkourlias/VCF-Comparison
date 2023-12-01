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
    parser = argparse.ArgumentParser(prog="VcfCompare")
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
    with open(output, 'w') as f:
        for chunk in range(start, stop, chunk_size):
            print(f"Comparing between {chunk} - {chunk+chunk_size}")
            print(f"{(stop-chunk)/chunk_size} iterations left.")

            i_records = comp_i.query(chr, chunk, chunk+chunk_size)
            i_records = [rec for rec in i_records]

            c_records = comp_f.query(chr, chunk, chunk+chunk_size)

            f.write("varID\t" + "samplesize\t" + "isIndel\t" + "mafI\t" +"mafC\t" + "hweI\t" + "hweC\t" "pcc\n")
            for c_record in c_records:
                for i_record in i_records:
                    if i_record[1] == c_record[1] and i_record[3] == c_record[3] and i_record[4] == c_record[4]:
                        f.write(record_extract(i_record, c_record, list(headers.values())))
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
                    i += 1

                else:
                    stop = int(line.split("\t")[1])

    return start, stop

def record_extract(i_record, c_record, head_idx):
    # Init vars
    c_record = [entry for i, entry in enumerate(c_record) if i in head_idx]
    sampleSize = 0
    isIndel = 0

    # Extract
    varID = f"{i_record[0]}:{i_record[1]}:{i_record[3]}_{i_record[4]}" 
    
    for i, (entry_i, entry_c) in enumerate(zip(i_record, c_record)):
        match i:
            case 0 | 1 | 2 | 5 | 6 | 7:
                continue

            case 3 | 4:
                if len(entry_i) > 1:
                    isIndel = 1
            
            case 8:
                val_dict_i = {}
                val_dict_c = {}

                for ind_i, ind_c in zip(entry_i.split(":"), entry_c.split(":")):
                    val_dict_i[ind_i] = "./."
                    val_dict_c[ind_c] = "./."

            case other:
                if entry_i == "./.":
                    entry_i = "./." + (len(val_dict_c) -1) * ":./."
                for key_i, key_c, val_i, val_c in zip(val_dict_i, val_dict_c, entry_i.split(":"), entry_c.split(":")):
                    val_dict_i[key_i] = val_i
                    val_dict_c[key_c] = val_c

                if val_dict_i["GT"] != "./." and val_dict_c["GT"] != "./.":
                    vals = val_dict_i["GT"].split("/")
                    sampleSize += (int(vals[0]) +  int(vals[1]))
        
    return f"{varID}\t{sampleSize}\t{isIndel}\n"

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
