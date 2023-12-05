#!/usr/bin/env python3

"""
    usage:
        ./vcf_rompare.py -f [INPUT VCF FILE] -c [VCF COMPARISON FILE]
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
from scipy.stats import pearsonr

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
    parser.add_argument("-f", "--file", type=str)
    parser.add_argument("-r", "--reference", type=str)
    parser.add_argument("-ih", "--input_header", type=str)
    parser.add_argument("-ch", "--comp_header", type=str)
    parser.add_argument("-chr", "--chr", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-n", "--chunk", type=int)
    return parser.parse_args()


def open_tbi(input_rile):

    """
    Tabix file opener function.
    :input_file: Input VCF file.
    :return: Opened tabix file.
    """
    return tabix.open(input_rile)


def compare(chr, start, stop, recs_i, r_handle, headers, chunksize):

    """
    Takes chromosome entry from nextflow process and compare record entries between two files.
    i in var name = VCF File containing all records for which you want to find hits.
    c in var name = vCF File to be compard with.
    :chr: chromosome number in (chr1) format.
    :i_handle: Input VCF file handle.
    :r_handle: Comparison VCF file handle.
    :input_headers_r: Input VCF header file handle.
    :comp_headers_r: Comparison VCF header file handle.
    :re
    turn: Pandas dataframes containing records with intersectioned positions.
    """
    # Open and query tabix files.
    comp_r = open_tbi(r_handle)
    i_records = recs_i
    results = []

    # CLI feedback.
    # print(f"Comparing between {chunk} - {chunk+chunk_size}")
    # print(f"{(stop-chunk)/chunk_size} iterations left.")

    # Query the chunk in the comparison file for a generator.
    r_records = comp_r.query(chr, start, stop)
    # For comparison record in the generator.
    for r_record in r_records:
        # For each retrieved record check whether there's a match in the input record list.
        for i_record in i_records:
            if (
                i_record[1] == r_record[1]  # If CHR is the same.
                and i_record[3] == r_record[3]  # If ref is the same.
                and i_record[4] == r_record[4]  # If alt is the same.
            ):
                results.append(
                    record_extract(i_record, r_record, list(headers.values()))
                )
    # Comparison

    return results


def header_subset(input_headers_r, comp_headers_r):
    """Looks for header intersections between files and retrieves
    the relevant indexes.

    Args:
        input_headers_r (path): File path for input headers.
        comp_headers_r (path): File path for input headers.

    Returns:
        dictionary: Dictionary of indexes corresponding to intersecting headers between vcf files.
    """

    # Initialise an empty dictionary.
    head_idx = {}

    # For each input header, check whether there is a match in comparison headers.
    for header_i in np.genfromtxt(input_headers_r, delimiter="\t", dtype=str):
        for i, header_r in enumerate(
            np.genfromtxt(comp_headers_r, delimiter="\t", dtype=str)
        ):  # If match is found, save the index in the dictionary.
            if header_i == header_r:
                head_idx[header_i] = i
                break

    return head_idx


def chunk_iterate(chr, i_handle, r_handle, headers, output, chunksize=2):
    """Retrieves the start and stop positions of relevant vcf
    selection.

    Args:
        chr (string): Chromosome string.
        i_handle (path): Input vcf file handle.

    Returns:
        string, string: number strings of start and stop range.
    """
    # Get start and stop positions of chr selection.
    with open(output, "w") as fo:

        # Write headers to output file.
        fo.write(
            "var_id\t"
            + "samplesize\t"
            + "is_indel\t"
            + "maf_i\t"
            + "maf_r\t"
            + "hwe_i\t"
            + "hwe_r\t"
            "pcc\n"
        )

        with gzip.open(i_handle, "rt") as f:
            recs = [0 for _ in range(chunksize)]
            i = -1
            for line in f:
                if line.startswith(chr):
                    recs[i] = line.split("\t")
                    if i == 0:
                        start = recs[i][1]
                        stop = int(recs[i][1]) + 1
                    else:
                        stop = recs[i][1]

                    i += 1

                    if len([rec for rec in recs if rec != 0]) == chunksize:
                        results = compare(
                            chr,
                            int(start),
                            int(stop),
                            recs,
                            r_handle,
                            headers,
                            chunksize,
                        )
                        for res in results:
                            fo.write(res)
                        recs = [0 for _ in range(chunksize)]
                        i = 0

            # process trailing records
            if i != 0:
                recs = [rec for rec in recs if rec != 0]
                results = compare(
                    chr,
                    int(start),
                    # Adding 10 in case of n = 1. Tabix cannot deal with single position queries.
                    int(stop) + 10,
                    recs,
                    r_handle,
                    headers,
                    chunksize,
                )
                
                # Write trailing results in 
                for res in results:
                    fo.write(res)

    return 0


def record_extract(i_record, r_record, head_idx):
    """Extracts the relevant intersection information from two tabix records.

    Args:
        i_record (list): list containing record entries for input file.
        r_record (list): list containing record entries for comparison file.
        head_idx (dict): header indexes

    Returns:
        string: string to be written to output file.
    """
    # Get a column value for each row if header/column index is present in both files.
    r_record = [entry for i, entry in enumerate(r_record) if i in head_idx]
    # Init new vars.
    sampleSize = 0
    isIndel = 0

    nrhoma_i = 0
    nrhets_i = 0
    nrhomb_i = 0

    nrhoma_r = 0
    nrhets_r = 0
    nrhomb_r = 0

    gt_list_r = []
    gt_list_i = []

    # Construct identifier string for the record/row.
    var_id = f"{i_record[0]}:{i_record[1]}:{i_record[3]}_{i_record[4]}"

    # Selection is same length because of column index selection.
    # Zip to iterate over both record values and track index.
    for i, (entry_i, entry_r) in enumerate(zip(i_record, r_record)):
        match i:  # Match the index value.
            case 0 | 1 | 2 | 5 | 6 | 7:  # Values on these indexes are not relevant.
                continue

            case 3 | 4:  # Check whether ref or alt contains an indel.
                if len(entry_i) > 1:
                    isIndel = 1

            case 8:  # Initialse the format column indicators/headers.
                val_dict_i = {}
                val_dict_r = {}

                # Store format indication/header into a dictionary with empty values.
                for ind_i, ind_r in zip(entry_i.split(":"), entry_r.split(":")):
                    val_dict_i[ind_i] = "./."
                    val_dict_r[ind_r] = "./."

            case other:  # Sample rows.
                # Assign values in format order.
                for key_i, key_r, val_i, val_r in zip(
                    val_dict_i, val_dict_r, entry_i.split(":"), entry_r.split(":")
                ):
                    val_dict_i[key_i] = val_i
                    val_dict_r[key_r] = val_r

                gt_i = val_dict_i["GT"]
                gt_r = val_dict_r["GT"]

                vals_i = gt_i.split("/")
                vals_r = gt_r.split("/")

                if (
                    vals_i[0] != "."
                    and vals_r[0] != "."
                    and int(vals_i[0]) + int(vals_i[1])
                    == int(vals_r[0]) + int(vals_r[1])
                ):
                    sampleSize += 1

                # If genotype is not missing, sum alternative alele counts.
                nrhoma_i, nrhets_i, nrhomb_i, gt_list_i = al_counts(
                    vals_i, nrhoma_i, nrhets_i, nrhomb_i, gt_list_i
                )
                nrhoma_r, nrhets_r, nrhomb_r, gt_list_r = al_counts(
                    vals_r, nrhoma_r, nrhets_r, nrhomb_r, gt_list_r
                )

    # Get the MAF
    maf_i = get_maf(nrhoma_i, nrhets_i, nrhoma_i)
    maf_r = get_maf(nrhoma_r, nrhets_r, nrhoma_r)

    # Get HWE
    hwe_i = get_hwe(nrhoma_i, nrhets_i, nrhomb_i)
    hwe_r = get_hwe(nrhoma_r, nrhets_r, nrhomb_r)

    # Lazy pearson corr implementation. Maybe change nr vars to be lists.
    corr, pval = pearsonr(gt_list_r, gt_list_r)
    # Return a constructed row string for immediate output writing.
    return f"{var_id}\t{sampleSize}\t{isIndel}\t{maf_i}\t{maf_r}\t{hwe_i}\t{hwe_r}\t{corr}\n"


def al_counts(vals, nrhoma, nrhets, nrhomb, gt_list):
    """Counts alele frequencies.

    Args:
        vals (_type_): _description_
        nrhoma (_type_): _description_
        nrhets (_type_): _description_
        nrhomb (_type_): _description_

    Returns:
        _type_: _description_
    """

    # Make List 

    if vals[0] != ".":
        num = int(vals[0]) + int(vals[1])
        gt_list.append(num)
        match num:
            case 0:
                nrhoma += 1

            case 1:
                nrhets += 1

            case 2:
                nrhomb += 1

    return nrhoma, nrhets, nrhomb, gt_list


def get_maf(nrhoma, nrhets, nrhomb):
    # MAF Calculation
    sampleSize = nrhoma + nrhomb + nrhets
    maf = min(nrhoma, nrhomb, nrhets) / sampleSize
    return maf


def get_hwe(nrhoma, nrhet, nrhomb):
    # Init vars.
    # AA, Aa, aa
    dom, het, rec = nrhomb + 1, nrhet + 1, nrhoma + 1
    pop = rec + dom + het

    # Alele freqs.
    p = (2 * dom + het) / (2 * (dom + het + rec))
    q = 1 - p

    # Expected HW.
    exp_dom = p**2 * pop
    exp_het = 2 * p * q * pop
    exp_rec = q**2 * pop

    # Chi-squared for deviation.
    hwe = (
        (((dom - exp_dom) ** 2) / exp_dom)
        + (((het - exp_het) ** 2) / exp_het)
        + (((rec - exp_rec) ** 2) / exp_rec)
    )
    return hwe


# MAIN
def main(args):
    """Main function"""
    # Get args.
    args = arg_parse()
    # Get headers
    headers = header_subset(args.input_header, args.comp_header)
    # Perform comparison and write output.
    chunk_iterate(args.chr, args.file, args.reference, headers, args.output, args.chunk)
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
