#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
import pandas as pd


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to convert a PLINK formatted phenotype file to a Pyseer formatted phenotype file.\n" \
                  "See https://zzz.bwh.harvard.edu/plink/data.shtml - Alternate phenotype files, for information on " \
                  "PLINK phenotype file. " \
                  "See https://pyseer.readthedocs.io/en/master/usage.html#phenotype-and-covariates for information on " \
                  "Pyseer phenotype files.\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-p", "--plink_phenotype_file", action="store", dest="plink_phenotype_file",
        help="Input PLINK-formatted alternate phenotype files",
        required=True, metavar="PLINK_PHEN"
    )
    group.add_argument(
        "-o", "--pyseer_phenotype_file", action="store", dest="pyseer_phenotype_file",
        help="Output Pyseer-formatted phenotype file",
        required=True, metavar="PYSEER_PHEN"
    )
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-s", "--samples_list", action="store", dest="samples_list",
        help="Samples to keep from input PLINK-formatted phenotype file.\nOne sample id per line expected.\n",
        required=False, metavar="SAMPLES_LIST"
    )
    group.add_argument(
        "-c", "--plink_phenotype_coding", action="store", dest="plink_phenotype_coding",
        help="Plink phenotype coding of missing/unaffected/affected (e.g. -9/0/1).\nDefault: -9/1/2.\n",
        required=False, default="-9/1/2", metavar="PHEN_CODING"
    )
    group.add_argument(
        "-n", "--pyseer_column_names", action="store", dest="pyseer_column_names",
        help="Column names of output Pyseer-formatted phenotype file. By default, they are named as: "
             "samples,phenotype1,phenotype2,...\n",
        required=False, metavar="COLUMN_NAMES"
    )

    return parser.parse_args()


def check_file_exists(my_file):
    if not os.path.isfile(my_file):
        logging.error(f'File {my_file} not found!')
        sys.exit(-1)


def count_nan(input_list):
    num_nan = input_list.isna().sum()
    return num_nan


# ------------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------------

def _main():
    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    # Making sure input files exist
    logging.info('Making sure input files exist...')
    input_files = [args.plink_phenotype_file]
    if args.samples_list is not None:
        input_files.append(args.samples_list)
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # 0. Parsing optional arguments
    plink_coding = args.plink_phenotype_coding
    if len(plink_coding.split('/')) != 3:
        logging.error(f'Chosen --plink_phenotype_coding \'{args.plink_phenotype_coding}\' does not contain three values')
        sys.exit(-1)

    # 1. Loading input PLINK-formatted phenotype files
    logging.info(f"Loading input PLINK-formatted phenotype files {args.plink_phenotype_file}...")
    pheno = pd.read_csv(args.plink_phenotype_file, sep=' ', header=None)
    print("Number of samples: " + str(pheno.shape[0]) + ". Number of phenotypes: " + str(pheno.shape[1]-2))
    logging.info(f"Loading input PLINK-formatted phenotype files {args.plink_phenotype_file}. DONE")

    # 2. Editing PLINK-formatted phenotype content (as dataframe)
    logging.info(f"Changing formats...")
    # 2.1 Removing columns where all elements are NaN
    columns_num_nan = pheno.apply(count_nan, axis=0).tolist()
    for idx, num_nan in enumerate(columns_num_nan):
        if num_nan == pheno.shape[0]:
            pheno.drop(pheno.columns[idx], axis=1, inplace=True)
            logging.warning(f"Column {idx+1}th removed as contained NaN values")

    # 2.2 Removing first Family ID column
    pheno.drop(pheno.columns[0], axis=1, inplace=True)

    # 2.3 Replacing column names, as default or as specified in input
    column_names_old = pheno.columns.values.tolist()
    if args.pyseer_column_names is not None:
        column_names_new = args.pyseer_column_names.strip().split(',')
        if len(column_names_new) != pheno.shape[1]:
            logging.error(f'Length of chosen --pyseer_column_names \'{args.pyseer_column_names}\' does match expected '
                          f'number of columns')
            sys.exit(-1)
    else:
        column_names_new = ['samples']
        for i in range(1, pheno.shape[1]):
            column_names_new.append('phenotype' + str(i))
    pheno.rename(columns=dict(zip(column_names_old, column_names_new)), inplace=True)

    # 2.4 Replacing coded phenotypes, from -9,1,2 to NA,0,1 for non-available,control,case
    [na_c, con_c, cas_c] = plink_coding.split('/')
    for column in pheno.columns[1:]:
        pheno[column].replace(int(na_c), 'NA', inplace=True)
        pheno[column].replace(int(con_c), '0', inplace=True)
        pheno[column].replace(int(cas_c), '1', inplace=True)

    # 3. Keeping samples specified in input. Reading input sample list if provided
    samples_to_keep = dict()
    if args.samples_list is not None:
        with open(args.samples_list, 'r') as lines:
            for line in lines:
                samples_to_keep[line.strip()] = 0
    if len(samples_to_keep) > 0:
        logging.info(f"Keeping samples specified in {args.samples_list}")
        rows_to_keep = list()
        for idx, sample in enumerate(pheno['samples']):
            if sample in samples_to_keep:
                rows_to_keep.append(idx)
        pheno = pheno.loc[rows_to_keep]

    # 4. Saving dataframe as table
    logging.info(f"Saving Pyseer-formatted file into {args.pyseer_phenotype_file}")
    pheno.to_csv(args.pyseer_phenotype_file, header=True, sep='\t', index=False)


if __name__ == "__main__":
    _main()



