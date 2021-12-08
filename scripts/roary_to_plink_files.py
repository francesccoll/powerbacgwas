#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
import pandas as pd


# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# Script tested on:
#   PLINK v1.90b6.17 64-bit (28 Apr 2020)

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to convert a Roary gene_presence_absence.Rtab file to input PLINK binary PED files " \
                  "(.fam, .bim and .bed )\n" \
                  "See https://www.cog-genomics.org/plink/2.0/input#bed for PLINK input information.\n" \
                  "See https://www.cog-genomics.org/plink/1.9/formats#bed for PLINK input format information.\n" \
                  "PLINK input-formatted files are needed to run GCTA phenotype simulations\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-g", "--gene_presence_absence", action="store", dest="gene_presence_absence",
        help="Roary gene_presence_absence.Rtab file",
        required=True, metavar="INPUT_TABLE"
    )
    group.add_argument(
        "-f", "--input_format", action="store", dest="input_format",
        help="Roary file format: (1) gene_presence_absence.csv or (2) gene_presence_absence.Rtab",
        required=True, metavar="INPUT_FORMAT"
    )
    group.add_argument(
        "-o", "--output_prefix", action="store", dest="output_prefix",
        help="output prefix used to name PLINK formatted files",
        required=True, metavar="OUTPUT_PREFIX"
    )
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-s", "--samples_list", action="store", dest="samples_list",
        help="Samples to keep from Roary gene_presence_absence.csv file.\nOne sample id per line expected.\n",
        required=False, metavar="SAMPLES_LIST"
    )
    group.add_argument(
        "-p", "--plink_path", action="store", dest="plink_path",
        help="Full path to PLINK executable", required=False, default='plink', metavar="PLINK_PATH"
    )

    return parser.parse_args()


def check_dependency(executable_name):
    """ Returns true if executable exists, else false """
    found = False
    try:
        output = subprocess.check_output(['which', executable_name]).strip()
        if output:
            found = True
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)
    return found


def check_file_exists(my_file):
    if not os.path.isfile(my_file):
        logging.error(f'File {my_file} not found!')
        sys.exit(-1)


def run_command_shell_string(command_line_string):
    """
    This function executes a command line, check for execution errors and but does not return stdout
    This is to be used when the stdout is not needed
    Note: shell=True needs to be set if I/O redirection operators are to be used (e.g. >) in the command line,
    otherwise they will have no special meaning, they are treated as ordinary arguments
    Note: if shell=True is used then the command line must be provided as a string, not a list
    :param command_line_string: it must be a string not a list
    """
    print('Running\n' + command_line_string)
    try:
        subprocess.run(command_line_string,
                       check=True,
                       shell=True,
                       )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)


def assign_allele(x):
    x.replace('"', '')
    allele = '0' if x == '' else '1'
    return allele


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

    # Making sure dependencies exist
    logging.info('Making sure dependencies exist...')
    dependencies = [args.plink_path]
    for dependency in dependencies:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)

    # Making sure input_format is either 1 or 2
    if args.input_format != '1':
        if args.input_format != '2':
            logging.error(f'Input format chosen (-f, --input_format) \'{args.input_format}\' must be either 1 or 2.')
            sys.exit(-1)

    # Making sure input files exist
    logging.info('Making sure input files exist...')
    input_files = [args.gene_presence_absence]
    if args.samples_list is not None:
        input_files.append(args.samples_list)
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # 1. Converting Roary file into made-up VCF

    # 1.1 Loading gene_presense_absense.Rtab as dataframe
    logging.info(f"Reading Roary file {args.gene_presence_absence}...")
    df = pd.read_csv(args.gene_presence_absence, sep='\t')
    logging.info(f"Reading Roary file {args.gene_presence_absence}. DONE")

    # 1.2 Inserting new VCF columns (CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT)
    logging.info(f"Converting Roary file {args.gene_presence_absence} into intermediate VCF...")
    if args.input_format == '1':
        gene_names = df.iloc[:, 0].tolist()
        gene_names = [g.replace('"', '') for g in gene_names]
        df2 = df.iloc[:, 14:int(df.shape[1])]
        df2.applymap(assign_allele)
        samples_names = list(df.columns)
        samples_names = samples_names[14:len(samples_names)]
    if args.input_format == '2':
        gene_names = df.iloc[:, 0].tolist()
        df2 = df.iloc[:, 1:int(df.shape[1])]
        samples_names = list(df.columns)
        samples_names.pop(0)
    format_list = ['GT'] * int(df.shape[0])
    df2.insert(0, 'FORMAT', format_list)  # FORMAT column
    dots_list = ['.'] * int(df.shape[0])
    df2.insert(0, 'INFO', dots_list)  # INFO column
    df2.insert(0, 'FILTER', dots_list)  # FILTER column
    df2.insert(0, 'QUAL', dots_list)  # QUAL column
    alt_list = ['1'] * int(df.shape[0])
    df2.insert(0, 'ALT', alt_list)
    ref_list = ['0'] * int(df.shape[0])
    df2.insert(0, 'REF', ref_list)
    id_list = gene_names
    df2.insert(0, 'ID', id_list)
    pos_list = list(range(1, int(df.shape[0])+1, 1))
    df2.insert(0, 'POS', pos_list)
    chr_list = ['1'] * int(df.shape[0])
    df2.insert(0, 'CHROM', chr_list)
    print("Number of samples: " + str(len(samples_names)) + ". Number of variants: " + str(df2.shape[0]))

    # 1.3 Creating made-up VCF header
    vcf_header_lines = ['##fileformat=VCFv4.1', '##FILTER=<ID=PASS,Description="All filters passed">',
                        '##contig=<ID=1,length='+str(df.shape[0]) +'>',
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']

    # 1.4 Printing made-up vcf
    tmp_vcf_file = open(args.output_prefix + '_tmp.vcf', 'w')
    for vcf_header_line in vcf_header_lines:
        tmp_vcf_file.write(vcf_header_line + '\n')
    vcf_header_line = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples_names)
    tmp_vcf_file.write(vcf_header_line + '\n')
    for index, row in df2.iterrows():
        new_vcf_line = row.tolist()
        new_vcf_line = [str(i) for i in new_vcf_line]
        tmp_vcf_file.write('\t'.join(new_vcf_line) + '\n')
    logging.info(f"Converting Roary file {args.gene_presence_absence} into intermediate VCF. DONE.")

    # 2. Running PLINK to convert made-up VCF file into BED files
    logging.info('Running PLINK to convert intermediate VCF file to BED format...')
    plink_command = ['plink', '--vcf', args.output_prefix + '_tmp.vcf', '--make-bed', '--double-id',
                     '--out', args.output_prefix]
    run_command_shell_string(' '.join(plink_command))
    logging.info('Running PLINK to convert intermediate VCF file to BED format. DONE.')

    # 3, Removing temporary files
    # logging.info('Removing temporary files.')
    # rm_command = ['rm', args.output_prefix + '_tmp.vcf']
    # run_command_shell_string(' '.join(rm_command))


if __name__ == "__main__":
    _main()
