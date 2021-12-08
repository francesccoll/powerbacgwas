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
        "-r", "--roary_table", action="store", dest="roary_table",
        help="Roary gene_presence_absence.Rtab file",
        required=True, metavar="ROARY_TABLE"
    )
    group.add_argument(
        "-f", "--input_format", action="store", dest="input_format",
        help="Roary file format: (1) gene_presence_absence.csv or (2) gene_presence_absence.Rtab",
        required=True, metavar="INPUT_FORMAT"
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
    if not os.path.exists(args.roary_table):
        logging.error(f'Input file in --roary_table {args.roary_table} not found!')
        sys.exit(-1)
    if not os.path.exists(args.causal_loci):
        logging.error(f'Input file in --causal-loci {args.causal_loci} not found!')
        sys.exit(-1)

    # Making sure subset_size is a positive intenger
    check_positive_integer(args.subset_size, '--subset_size', None, None)

    # Making sure input_format is either 1 or 2
    if args.input_format != '1':
        if args.input_format != '2':
            logging.error(f'Input format chosen (-f, --input_format) \'{args.input_format}\' must be either 1 or 2.')
            sys.exit(-1)

    # Saving variants in causal variants file.
    logging.info(f"Reading causal variant file {args.causal_loci}...")
    causal_variants = dict()
    file = open(args.causal_loci, "r")
    for line in file:
        [var_id, *_] = line.strip().split('\t')
        causal_variants[var_id] = ""
    print("\tNumber of causal variants extracted: " + str(len(causal_variants)))

    # Opening gene_presence_absence.csv file and saving gene frequencies
    logging.info(f"Opening {args.roary_table} file and calculating gene frequency...")
    gene_calls = dict()  # dictionary to save gene presence/absence calls: gene_calls[gene_id] = presabs
    num_genes = 0
    roary_samples = list()
    f = open(args.roary_table, 'r')
    if args.input_format == '1':
        fields = f.readline().strip().split('","')  # sample ids are extracted from first line
        samples = fields[14:len(fields)]
        samples = [s.replace('"', '') for s in samples]
        roary_samples = samples
        for line in f:
            num_genes = num_genes + 1
            fields = line.strip().split('","')
            gene_id = fields[0].replace('"', '')
            presabs = fields[14:len(fields)]
            presabs = [allele.replace('"', '') for allele in presabs]
            gene_calls[gene_id] = presabs
            # print("gene_calls[" + gene_id + "] --> ", gene_calls[gene_id])
    if args.input_format == '2':
        fields = f.readline().strip().split('\t')  # sample ids are extracted from first line
        samples = fields
        samples.pop(0)
        roary_samples = samples
        for line in f:
            num_genes = num_genes + 1
            fields = line.strip().split('\t')
            gene_id = fields[0]
            presabs = fields
            presabs.pop(0)
            gene_calls[gene_id] = presabs
            # print("gene_calls[" + gene_id + "] --> ", gene_calls[gene_id])
    f.close()
    print("\tNumber of genes: ", str(num_genes))
    print("\tNumber of genes saved: ", str(len(gene_calls)))
    print("\tNumber of samples: ", str(len(roary_samples)))
    logging.info(f"Opening gene_presence_absence file and calculating gene frequency. DONE.")

    # Make sure size of subset (--subset_size) is not greater than number of samples in VCF
    if int(args.subset_size) > len(roary_samples):
        logging.error(f'subset_size ({args.subset_size}) is greater than number of samples in Roary file '
                      f'({str(len(roary_samples))})!')
        sys.exit(-1)

    # Extracting samples with and without causal variant (mutated and wild-type)
    logging.info(f"Extracting samples with and without causal variant(s) (mutated and wild-type)...")
    roary_calls = dict()
    for var_id in causal_variants:
        variant_calls = list()
        if var_id in gene_calls:
            variant_calls = gene_calls[var_id]
        else:
            logging.error(f'Gene {var_id} in {args.causal_loci} not found in {args.roary_table}')
            sys.exit(-1)
        if len(roary_samples) != len(variant_calls):
            logging.error(f'Number of sample ids extracted from Roary {str(len(roary_samples))} does not match number of '
                          f'genotype calls {str(len(variant_calls))} for {var_id}. Check Roary format.')
            sys.exit(-1)
        roary_calls[var_id] = variant_calls
    roary_samples_mut = []  # list to store samples with causal gene (mutated)
    roary_samples_wt = []  # list to store samples without causal gene (wild-type)
    for var_id in causal_variants:
        variant_calls = roary_calls[var_id]
        for idx, sample in enumerate(roary_samples):
            allele = str(variant_calls[idx])
            # if allele is '0':
            #     roary_samples_wt.append(sample)
            if allele == '1':
                roary_samples_mut.append(sample)
    roary_samples_mut = list(set(roary_samples_mut))
    roary_samples_wt = list(set(roary_samples) - set(roary_samples_mut))
    if len(roary_samples_mut) == 0:
        logging.error(f'Number of samples with causal variant is zero. Check VCF format, and allele coding.')
        sys.exit(-1)
    if len(roary_samples_wt) == 0:
        logging.error(f'Number of samples without causal variant is zero. Check VCF format, and allele coding.')
        sys.exit(-1)
    print("\t\tNumber of mutated individuals available to sample from: " + str(len(roary_samples_mut)))
    print("\t\tNumber of wild-type individuals available to sample from: " + str(len(roary_samples_wt)))
    allele_frequency = len(roary_samples_mut)/(len(roary_samples_mut)+len(roary_samples_wt))
    print("\t\tObserved allele frequency: " + str(allele_frequency))

    # Extracting sub-sampled samples with and without causal variant
    # NOTE: calculating the number of samples to sub-sample, with (n1) and without (n2) causal variants, to
    # retain the same allele frequency in the sub-sample
    logging.info(f'Randomnly sampling {str(args.subset_size)} from {args.roary_table}')
    n1 = allele_frequency * int(args.subset_size)
    n2 = (1 - allele_frequency) * int(args.subset_size)
    n1 = math.ceil(n1)
    n2 = math.ceil(n2)
    sampled_roary_samples_mut = random.sample(roary_samples_mut, n1)
    sampled_roary_samples_wt = random.sample(roary_samples_wt, n2)
    sampled_roary_samples_all = sampled_roary_samples_mut + sampled_roary_samples_wt

    # Saving subsampled list
    logging.info(f'Writing sampled items into {args.output_samples}')
    output = open(args.output_samples, 'w')
    for sample in sampled_roary_samples_all:
        output.write(sample + '\n')
    output.close()


if __name__ == "__main__":
    _main()