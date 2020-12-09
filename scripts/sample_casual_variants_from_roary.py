#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import random
import gzip
import numpy


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to sample causal genes from a Roary gene_presence_absence.csv file.\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-v", "--input_roary", action="store", dest="input_roary",
        help="input Roary gene_presence_absence.csv file",
        required=True, metavar="INPUT_ROARY"
    )
    group.add_argument(
        "-f", "--input_format", action="store", dest="input_format",
        help="Roary file format: (1) gene_presence_absence.csv or (2) gene_presence_absence.Rtab",
        required=True, metavar="INPUT_FORMAT"
    )
    group.add_argument(
        "-p", "--pastml_steps_file", action="store", dest="pastml_steps_file",
        help="table with number of homoplasies per gene created by script ancestral_state_reconstruction_roary.py",
        required=True, metavar="PASTML_STEPS"
    )
    group.add_argument(
        "-n", "--number_variants", action="store", dest="number_variants",
        help="number of randomly sampled genes (integer)",
        required=True, metavar="NUMBER_VARIANTS"
    )
    group.add_argument(
        "-a", "--allele_frequency_range", action="store", dest="allele_frequency_range",
        help="range of frequency of selected genes, separated by comma (e.g. 0.05,0.10 for 5%% to 10%%)",
        required=True, metavar="ALLELE_FREQUENCY_RANGE"
    )
    group.add_argument(
        "-s", "--degree_homoplasy_range", action="store", dest="degree_homoplasy_range",
        help="range of degree of homoplasy of selected genes (number of steps), separated by comma "
             "(e.g. 2,5 for 2 to 5 steps)",
        required=True, metavar="DEGREE_HOMOPLASY_RANGE"
    )
    group.add_argument(
        "-e", "--effect_size", action="store", dest="effect_size",
        help="effect size of selected variants, in odd ratio units (e.g. 2 for OR of 2)",
        required=True, metavar="EFFECT_SIZE"
    )
    group.add_argument(
        "-o", "--output_gcta_file", action="store", dest="output_gcta_file",
        help="randomly sampled Roary genes in GCTA format",
        required=True, metavar="OUTPUT_TABLE"
    )

    return parser.parse_args()


def get_vcf_reader(my_vcf):
    if os.path.splitext(my_vcf)[-1].lower() == '.gz':
        return vcf.Reader(open(my_vcf, 'rb'))
    else:
        return vcf.Reader(open(my_vcf, 'r'))


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


def check_positive_float(value, name, minim, maxim):
    """
    Function to check that a value is a positive numberical value (float)
    :param value: value of variable to check
    :param name: message about variable to print
    :param minim: minimum value the variable should be within
    :param maxim: maximum value the variable should be within
    :return:
    """
    try:
        val = float(value)
        if val < 0:
            logging.error(f"{name} provided ''{value}'' must be a positive value")
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
        logging.error(f"{name} provided '{value}' is not a numerical value")
        sys.exit(-1)


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
    input_files = [args.input_roary, args.pastml_steps_file]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Making sure input_format is either 1 or 2
    if args.input_format is not '1':
        if args.input_format is not '2':
            logging.error(f'Input format chosen (-f, --input_format) \'{args.input_format}\' must be either 1 or 2.')
            sys.exit(-1)

    # Parsing parameters
    (af_from, af_to) = args.allele_frequency_range.split(',')
    (steps_from, steps_to) = args.degree_homoplasy_range.split(',')
    # (or_from, or_to) = args.effect_size_range.split(',')
    odds_ratio = args.effect_size
    number_variants = args.number_variants.strip()

    # Making sure input parameters meet expected
    int_param = [steps_from, steps_to, number_variants]
    int_param_names = ['start of range in number of steps (-s)', 'end of range in number of steps (-s)',
                       'number of variants (-n)']
    for idx, param in enumerate(int_param):
        check_positive_integer(param, int_param_names[idx], None, None)
    float_param = [af_from, af_to]
    float_param_names = ['start of range in allele frequency (-f)', 'end of range in allele frequency (-f)']
    for idx, param in enumerate(float_param):
        check_positive_float(param, float_param_names[idx], 0, 1)
    check_positive_float(odds_ratio, "effect size (--effect_size)", None, None)

    # Printing chosen parameters
    print("Chosen variants parameters: ")
    print("\tNumber of causal variants: " + number_variants)
    print("\tAllele frequency: " + af_from + "-" + af_to)
    print("\tDegree of homoplasy: " + steps_from + "-" + steps_to)
    print("\tEffect size: " + odds_ratio)

    # Opening PastML steps output
    logging.info(f"Saving number of homoplasies (steps) per gene...")
    allele_steps = dict()    # dictionary to save allele frequency: allele_frequency[var_id] = sample_allele
    file = open(args.pastml_steps_file, "r")
    file.readline()
    for line in file:
        [var_id, steps] = line.strip().split('\t')
        allele_steps[var_id] = steps
        # print("allele_steps["+var_id+"] "+allele_steps[var_id])

    # Opening gene_presence_absence.csv file and saving gene frequencies
    logging.info(f"Opening {args.input_roary} file and calculating gene frequency...")
    allele_frequency = dict()  # dictionary to save allele frequency: allele_frequency[gene_id] = sample_allele
    num_genes = 0
    f = open(args.input_roary, 'r')
    if args.input_format is '1':
        fields = f.readline().strip().split('","')  # sample ids are extracted from first line
        samples = fields[14:len(fields)]
        samples = [s.replace('"', '') for s in samples]
        for line in f:
            num_genes = num_genes + 1
            fields = line.strip().split('","')
            gene_id = fields[0].replace('"', '')
            presabs = fields[14:len(fields)]
            frequency = 0
            for allele in presabs:
                allele.replace('"', '')
                if allele is not '':
                    frequency = frequency + 1
            allele_frequency[gene_id] = frequency / len(samples)
            # print("allele_frequency[" + gene_id + "] --> ", allele_frequency[gene_id])
    if args.input_format is '2':
        fields = f.readline().strip().split('\t')  # sample ids are extracted from first line
        samples = fields
        samples.pop(0)
        for line in f:
            num_genes = num_genes + 1
            fields = line.strip().split('\t')
            gene_id = fields[0]
            presabs = fields
            presabs.pop(0)
            frequency = 0
            for allele in presabs:
                if allele is not '0':
                    frequency = frequency + 1
            allele_frequency[gene_id] = frequency / len(samples)
            # print("allele_frequency[" + gene_id + "] --> ", allele_frequency[gene_id])
    f.close()
    print("\tNumber of genes: ", str(num_genes))
    print("\tNumber of genes saved: ", str(len(allele_frequency)))
    logging.info(f"Opening gene_presence_absence file and calculating gene frequency. DONE.")

    # Selecting variants meeting criteria (allele frequency and homoplasy steps)
    logging.info(f"Selecting variants meeting criteria...")
    set1 = set()
    for var_id, freq in allele_frequency.items():
        if (float(freq) >= float(af_from)) and (float(freq) <= float(af_to)):
            set1.add(var_id)
    set2 = set()
    for var_id, steps in allele_steps.items():
        if (int(steps) >= int(steps_from)) and (int(steps) <= int(steps_to)):
            set2.add(var_id)
    union_set = set1.intersection(set2)
    print(str(len(union_set)) + "/" + str(len(allele_frequency)) + " variants meeting criteria")
    logging.info(f"Selecting variants meeting criteria. DONE.")

    # Randomly sampling of variants
    logging.info(f"Randomly sampling of variants meeting criteria.")
    selected_variants = random.sample(union_set, int(number_variants))

    # Saving sampled causal variants in GCTA format
    # Input file format: causal.snplist (columns are SNP ID and effect size)
    output_gcta = open(args.output_gcta_file, 'w')
    for var_id in selected_variants:
        newline = var_id + '\t' + str(odds_ratio) + '\n'
        output_gcta.write(newline)
        print("\tselected_variant " + var_id)
        print("\tallele frequency " + str(allele_frequency[var_id]))
        print("\thomoplasy steps " + str(allele_steps[var_id]))
        print("\teffect size " + str(odds_ratio))
    output_gcta.close()


if __name__ == "__main__":
    _main()