#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import random
import vcf
import subprocess
import math


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to randomly subsample a subset of samples from a multi-sample VCF file, " \
                  "containing a specific set of variants.\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('I/O arguments')
    group.add_argument(
        "-v", "--input_vcf", action="store", dest="input_vcf",
        help="multi-sample VCF file",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-c", "--causal-loci", action="store", dest="causal_loci",
        help="VCF variants in GCTA causal.snplist",
        required=True, metavar="CAUSAL_LOCI"
    )
    group.add_argument(
        "-s", "--subset_size", action="store", dest="subset_size",
        help="Total number of samples to subsample from --input_vcf",
        required=True, metavar="SIZE"
    )
    group.add_argument(
        "-o", "--output_samples", action="store", dest="output_samples",
        help="List of output sub-sampled samples.",
        required=True, metavar="OUT_LIST"
    )
    return parser.parse_args()


def check_positive_integer(value, name, minim, maxim):
    """
    Function to check that a value is a positive integer
    :param value: value of variable to check
    :param name: message about variable to print
    :param minim: minimum value the variable should be within
    :param maxim: maximum value the variable should be within
    :return:
    """
    try:
        val = int(value)
        if val < 0:
            logging.error(f"{name} provided ''{value}'' must be a positive integer")
            sys.exit(-1)
        if minim is not None:
            if val < minim:
                logging.error(f"{name} provided ''{value}'' must be greater than {minim}")
                sys.exit(-1)
        if maxim is not None:
            if val > maxim:
                logging.error(f"{name} provided ''{value}'' must be smaller than {maxim}")
                sys.exit(-1)
    except ValueError:
        logging.error(f"{name} provided '{value}' is not an integer")
        sys.exit(-1)


def get_vcf_reader(my_vcf):
    if os.path.splitext(my_vcf)[-1].lower() == '.gz':
        return vcf.Reader(open(my_vcf, 'rb'))
    else:
        return vcf.Reader(open(my_vcf, 'r'))


def run_command_shell_string(command_line_string):
    """
    This function executes a command line, check for execution errors and but does not return stdout
    This is to be used when the stdout is not needed
    Note: shell=True needs to be set if I/O redirection operators are to be used (e.g. >) in the command line,
    otherwise they will have no special meaning, they are treated as ordinary arguments
    Note: if shell=True is used then the command line must be provided as a string, not a list
    :param command_line_string: it must be a string not a list
    """
    print('\tRunning: ' + command_line_string)
    try:
        subprocess.run(command_line_string,
                       check=True,
                       shell=True,
                       )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)


def run_command_string(command_line_string):
    """
    This function executes a command line, check for execution errors and returns stdout
    :param command_line_string: it must be a string
    :return: stdout
    """
    print('\tRunning: ' + command_line_string)
    try:
        process_completed = subprocess.run(
            command_line_string,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            shell=True,
        )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)
    return process_completed.stdout.decode('utf-8')


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

    # Making sure input file exists
    if not os.path.exists(args.input_vcf):
        logging.error(f'Input file in --input_vcf {args.input_vcf} not found!')
        sys.exit(-1)
    if not os.path.exists(args.causal_loci):
        logging.error(f'Input file in --causal-loci {args.causal_loci} not found!')
        sys.exit(-1)

    # Making sure subset_size is a positive intenger
    check_positive_integer(args.subset_size, '--subset_size', None, None)

    # Saving variants in causal variants file.
    logging.info(f"Reading causal variant file {args.causal_loci}...")
    causal_variants = dict()
    file = open(args.causal_loci, "r")
    for line in file:
        [var_id, *_] = line.strip().split('\t')
        causal_variants[var_id] = ""
    print("\tNumber of causal variants extracted: " + str(len(causal_variants)))

    # Reading VCF to extract samples ids
    logging.info(f"Reading VCF file {args.input_vcf} to extract samples with and without causal variant...")
    vcf_reader = get_vcf_reader(args.input_vcf)
    vcf_samples = vcf_reader.samples
    print("\tNumber of samples extracted from VCF: " + str(len(vcf_samples)))

    # Make sure size of subset (--subset_size) is not greater than number of samples in VCF
    if int(args.subset_size) > len(vcf_samples):
        logging.error(f'subset_size ({args.subset_size}) is greater than number of samples in VCF '
                      f'({str(len(vcf_samples))})!')
        sys.exit(-1)

    # Creating uncompressed temporary VCF file
    tmp_vcf_file = args.input_vcf + '_tmp.vcf'
    logging.info(f"Creating uncompressed temporary VCF file {tmp_vcf_file}...")
    if os.path.splitext(args.input_vcf)[-1].lower() == '.gz':
        run_command_shell_string(
            ' '.join(["gunzip", "-c", args.input_vcf, ">", tmp_vcf_file])
        )
    else:
        run_command_shell_string(
            ' '.join(["cp", args.input_vcf, tmp_vcf_file])
        )
    if not os.path.exists(tmp_vcf_file):
        logging.error(f'Uncompressed temporary VCF file {tmp_vcf_file} not found. Something went wrong.')
        sys.exit(-1)

    # Extracting VCF line with causal variant
    logging.info(f"Extracting VCF line(s) with causal variant from {tmp_vcf_file}...")
    vcf_lines = dict()
    for var_id in causal_variants:
        vcf_line = run_command_string(
            ''.join(["cat ", tmp_vcf_file, " | ", " awk -F'\t' '{ if ($3==\"", var_id, "\") print $0}'"])
        )
        if vcf_line.strip() == '':
            logging.error(f'{var_id} not found in {tmp_vcf_file}. Make sure variant IDs are included.')
            sys.exit(-1)
        vcf_lines[var_id] = vcf_line

    # Extracting samples with and without causal variant (mutated and wild-type)
    logging.info(f"Extracting samples with and without causal variant(s) (mutated and wild-type)...")
    vcf_samples_mut = []  # list to store samples with causal variant(s) (mutated)
    vcf_samples_wt = []  # list to store samples without causal variants(s) (wild-type)
    for var_id in causal_variants:
        variant_calls = vcf_lines[var_id].strip().split('\t')[9:]
        if len(vcf_samples) != len(variant_calls):
            logging.error(f'Number of sample ids extracted from VCF {str(len(vcf_samples))} does not match number of '
                          f'genotype calls {str(len(variant_calls))} for {var_id}. Check VCF format.')
            sys.exit(-1)
        for idx, sample in enumerate(vcf_samples):
            allele = str(variant_calls[idx])
            if allele == '1':
                vcf_samples_mut.append(sample)
    vcf_samples_mut = list(set(vcf_samples_mut))
    vcf_samples_wt = list(set(vcf_samples) - set(vcf_samples_mut))
    if len(vcf_samples_mut) == 0:
        logging.error(f'Number of samples with causal variant is zero. Check VCF format, and allele coding.')
        sys.exit(-1)
    if len(vcf_samples_wt) == 0:
        logging.error(f'Number of samples without causal variant is zero. Check VCF format, and allele coding.')
        sys.exit(-1)
    print("\t\tNumber of mutated samples: " + str(len(vcf_samples_mut)))
    print("\t\tNumber of wild-type samples: " + str(len(vcf_samples_wt)))
    allele_frequency = len(vcf_samples_mut)/(len(vcf_samples_mut)+len(vcf_samples_wt))
    print("\t\tObserved allele frequency: " + str(allele_frequency))

    # Extracting sub-sampled samples with and without causal variant
    # NOTE: calculating the number of samples to sub-sample, with (n1) and without (n2) causal variants, to
    # retain the same allele frequency in the sub-sample
    logging.info(f'Randomnly sampling {str(args.subset_size)} from {args.input_vcf}')
    n1 = allele_frequency * int(args.subset_size)
    n2 = (1 - allele_frequency) * int(args.subset_size)
    n1 = math.ceil(n1)
    n2 = math.ceil(n2)
    sampled_vcf_samples_mut = random.sample(vcf_samples_mut, n1)
    sampled_vcf_samples_wt = random.sample(vcf_samples_wt, n2)
    sampled_vcf_samples_all = sampled_vcf_samples_mut + sampled_vcf_samples_wt

    # Saving subsampled list
    logging.info(f'Writing sampled items into {args.output_samples}')
    output = open(args.output_samples, 'w')
    for sample in sampled_vcf_samples_all:
        output.write(sample + '\n')
    output.close()

    # Removing temporary files
    logging.info('Removing temporary files.')
    rm_command = ['rm', tmp_vcf_file]
    run_command_shell_string(' '.join(rm_command))


if __name__ == "__main__":
    _main()