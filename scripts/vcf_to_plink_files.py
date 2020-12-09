#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf
import subprocess


# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# Script tested on:
#   PLINK v1.90b6.17 64-bit (28 Apr 2020)


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to convert a multi-sample VCF file to input PLINK binary PED files (.fam, .bim and .bed )\n" \
                  "See https://www.cog-genomics.org/plink/2.0/input#bed for PLINK input information.\n" \
                  "See https://www.cog-genomics.org/plink/1.9/formats#bed for PLINK input format information.\n" \
                  "PLINK input-formatted files are needed to run GCTA phenotype simulations\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-v", "--input_vcf", action="store", dest="input_vcf",
        help="multi-sample VCF file",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-o", "--output_prefix", action="store", dest="output_prefix",
        help="output prefix used to name PLINK formatted files",
        required=True, metavar="OUTPUT_PREFIX"
    )
    group.add_argument(
        "-p", "--plink_path", action="store", dest="plink_path",
        help="Full path to PLINK executable or executable name", required=False, default='plink', metavar="PLINK_PATH"
    )
    group.add_argument(
        "-b", "--bcftools_path", action="store", dest="bcftools_path",
        help="Full path to bcftools executable or executable name", required=False, default='bcftools', metavar="BCFTOOLS_PATH"
    )

    return parser.parse_args()


def get_vcf_reader(my_vcf):
    if os.path.splitext(my_vcf)[-1].lower() == '.gz':
        return vcf.Reader(open(my_vcf, 'rb'))
    else:
        return vcf.Reader(open(my_vcf, 'r'))


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
    dependencies = [args.plink_path, args.bcftools_path]
    for dependency in dependencies:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)

    # Making sure input files exist
    logging.info('Making sure input files exist...')
    input_files = [args.input_vcf]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Editing original VCF to change chromosome ids and add variant ids
    logging.info(f"Editing input VCF to change chromosome id and add variant ids...\n")
    vcf_tmp = args.input_vcf + '_tmp'  # temporary edited VCF file
    vcf_reader = get_vcf_reader(args.input_vcf)
    chr_id = ''  # extract chromosome id
    for header_line in vcf_reader._header_lines:
        if header_line.startswith("##contig"):
            chr_id = header_line.strip().split(',')[0].replace('##contig=<ID=', '')
    replace_command = [args.bcftools_path, ' view ', args.input_vcf, ' | ',
                       args.bcftools_path, ' annotate -Ov -x ID -I +\'%CHROM\.%POS\.%REF\.%ALT\'', ' | ',
                       ' sed \'s/^', chr_id, '/1/g\'', ' | ',
                       ' sed \'s/##contig=<ID=', chr_id, '/##contig=<ID=1/g\'', ' > ', vcf_tmp]
    run_command_shell_string(''.join(replace_command))
    logging.info(f"Editing input VCF to change chromosome id and add variant ids. DONE\n")

    # Running PLINK
    logging.info('Running PLINK to change VCF format to BED...')
    plink_command = [args.plink_path, '--vcf', vcf_tmp, '--make-bed', '--double-id', '--out', args.output_prefix]
    run_command_shell_string(' '.join(plink_command))
    logging.info('Running PLINK to change VCF format to BED. DONE.')

    # Removing temporary files
    logging.info('Removing temporary files.')
    rm_command = ['rm', vcf_tmp]
    run_command_shell_string(' '.join(rm_command))


if __name__ == "__main__":
    _main()