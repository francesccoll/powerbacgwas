#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf
import subprocess
from scipy.optimize import least_squares
import random
import pandas as pd
import math


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to simulate binary phenotypes from VCF variants.\n Designed to test the " \
                  "effect of changing odds ratio (effect size). A file with sampled causal variants obtained " \
                  "from sample_casual_variants_from_vcf.py is needed as input. Input VCF need to be formatted as " \
                  "a multi-sample VCF file. \n" \
                  "NOTE: Simulation of quantitative phenotypes is not supported, use " \
                  "simulate_phenotype_using_gcta.py script instead.\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('I/O arguments')
    group.add_argument(
        "-v", "--input_vcf", action="store", dest="input_vcf",
        help="multi-sample VCF file",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-o", "--out_prefix", action="store", dest="out_prefix",
        help="output prefix used to name simulated phenotype output file",
        required=True, metavar="OUTPUT_PREFIX"
    )
    group.add_argument(
        "-c", "--causal-loci", action="store", dest="causal_loci",
        help="Sampled VCF variants used as causal variants, generated by sample_casual_variants_from_vcf.py.",
        required=True, metavar="CAUSAL_LOCI"
    )
    group.add_argument(
        "-x", "--simu-cc", action="store", dest="simu_cc",
        help="Number of cases and controls, separated by comma (e.g. 100,200)",
        required=True, metavar="SIMU_CC"
    )
    group.add_argument(
        "-r", "--simu-rep", action="store", dest="simu_rep",
        help="Number of simulation repetitions. The default value is 1 if this option is not specified.",
        required=False, metavar="SIMU_REP", default="1"
    )
    group = parser.add_argument_group('Optional arguments: ')
    group.add_argument(
        "-a", "--allele-frequency", action="store", dest="allele_frequency",
        help="Allele frequency of causal variant. If not specified, observed allele frequency of causal variant will "
             "be used. If specified, it cannot be higher than the observed allele frequency",
        required=False, metavar="ALLELE_FREQ", default="NA"
    )
    group.add_argument(
        "-e", "--effect-size", action="store", dest="effect_size",
        help="Effect size of causal variant in odds ratio units. If specified, this value will be used instead of that "
             "extracted from --causal-loci file",
        required=False, metavar="EFFECT_SIZE", default="NA"
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


def solve_odds_ratio_function(variables, mut, wt, cases, controls, odds_ratio):
    """
    This function is used to calculate the number of cases and controls, with and without the causal variant,
    to achieve the chosen odds ratio, where:
        controls_mut is the number of controls with causal variant
        cases_mut is the number of cases with causal variant
        controls_wt is the number of controls without causal variant
        cases_wt is the number of cases without causal variant
    A system of functions (f1 to f5) need to minimised to zero
    :param variables: initial guess of variables
    :param mut: number of individuals/samples with causal variant
    :param wt: number of individuals/samples without causal variant
    :param cases: number of cases
    :param controls: number of controls
    :param odds_ratio: chosen odds ratio of causal variant
    :return:
    """
    controls_mut = variables[0]
    cases_mut = variables[1]
    cases_wt = variables[2]
    controls_wt = variables[3]
    f1 = controls_mut + cases_mut - mut
    f2 = controls_wt + cases_wt - wt
    f3 = cases_wt + cases_mut - cases
    f4 = controls_wt + controls_mut - controls
    f5 = odds_ratio * (cases_wt/controls_wt) - (cases_mut/controls_mut)
    return [f1, f2, f3, f4, f5]


def solve_odds_ratio_function_2(variables, num_mut, num_wt, allele_frequency, sample_size, odds_ratio):
    """
    Same as above, but where the exact case/control ratio or allele frequency may not be met. This is used in cases
    where the number of mutated or wildtype samples needed to achieve the desired allele frequency exceeds the number
    of mutated or wildtype samples available to sub-sample from.

    A system of functions (f1 to f5) need to minimised to zero, where each function states:
    f1: the number of mutated cases + controls cannot exceed # of samples with variant available to sub-sample from
    f2: the number of WT cases + WT controls cannot exceed the number of WT samples available to sub-sample from
    f3: the chosen/desired overall sample size must be met
    f4: the chosen/desired allele frequency must be met
    f5: the chosen/desired odds ratio must be met

    :param variables: initial guess of variables
    :param num_mut: number of samples with causal variant according to variant file
    :param num_wt: number of samples without causal variant according to variant file
    :param allele_frequency: chosen allele frequency
    :param sample_size: number of cases and controls
    :param odds_ratio: chosen odds ratio of causal variant
    :return:
    """
    controls_mut = variables[0]
    cases_mut = variables[1]
    cases_wt = variables[2]
    controls_wt = variables[3]
    f1 = (controls_mut + cases_mut) - num_mut*allele_frequency
    f2 = (controls_wt + cases_wt) - num_wt*(1-allele_frequency)
    f3 = controls_mut + controls_wt + cases_mut + controls_wt - sample_size
    f4 = float(((controls_mut + cases_mut)/sample_size)) - float(allele_frequency)
    f5 = odds_ratio * (cases_wt/controls_wt) - (cases_mut/controls_mut)
    return [f1, f2, f3, f4, f5]


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
    input_files = [args.input_vcf, args.causal_loci]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Parsing parameters
    (num_cases, num_controls) = args.simu_cc.split(',')
    num_rep = args.simu_rep

    # Making sure input parameters meet expected values
    int_param = [num_cases, num_controls, num_rep]
    int_param_names = ['number of cases in (--simu-cc)', 'number of controls in (--simu-cc)',
                       'number of simulation repetitions (--simu-rep)']
    for idx, param in enumerate(int_param):
        check_positive_integer(param, int_param_names[idx], None, None)
    sample_size = int(num_cases) + int(num_controls)

    if args.allele_frequency != 'NA':
        check_positive_float(args.allele_frequency, " allele frequency (--allele-frequency) ", 0, 1)
    if args.effect_size != 'NA':
        check_positive_float(args.effect_size, " effect size (--effect-size)", 0, None)

    # Reading causal variant file.
    logging.info(f"Reading causal variant file {args.causal_loci}...")
    causal_variants = dict()
    odds_ratio = ''
    file = open(args.causal_loci, "r")
    for line in file:
        [var_id, odds_ratio] = line.strip().split('\t')
        check_positive_float(odds_ratio, "odds ratio in " + args.causal_loci, 0, None)
        causal_variants[var_id] = odds_ratio

    # Printing chosen parameters
    print("Chosen parameters: ")
    print("\tCausal variant(s): " + ' '.join(causal_variants))
    print("\tOdds ratio (from causal variant file): " + str(odds_ratio))
    print("\tNumber of cases: " + str(num_cases))
    print("\tNumber of controls: " + str(num_controls))
    print("\tSample size: " + str(sample_size))

    # Reading VCF to extract samples ids
    logging.info(f"Reading VCF file {args.input_vcf} to extract samples with and without causal variant...")
    vcf_reader = get_vcf_reader(args.input_vcf)
    vcf_samples = vcf_reader.samples
    print("\tNumber of samples extracted from VCF: " + str(len(vcf_samples)))

    # Creating uncompressed temporary VCF file
    tmp_vcf_file = args.out_prefix + '_tmp.vcf'
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

    # Replacing allele frequency and odds ratio if options specified
    if args.allele_frequency != 'NA':
        if float(args.allele_frequency) > allele_frequency:
            logging.error(f'Chosen allele frequence (--allele-frequency) {args.allele_frequency} cannot be higher than '
                          f'the observed one {allele_frequency}')
            sys.exit(-1)
        allele_frequency = float(args.allele_frequency)
    if args.effect_size != 'NA':
        odds_ratio = float(args.effect_size)

    # Calculating number of cases and controls, with and without the causal variant, to achieve the chosen odds ratio
    logging.info(f"Calculating number of cases and controls, with and without the causal variant, to achieve the "
                 f"chosen odds ratio...")
    initialGuess = [1, 1, 1, 1]
    mut = allele_frequency * float(sample_size)
    mut = int(math.floor(mut))
    wt = int(sample_size) - mut
    cases = num_cases
    controls = num_controls
    print("\tInput variables:")
    print("\t\tNumber of mutated samples: " + str(mut))
    print("\t\tNumber of wild-type samples: " + str(wt))
    print("\t\tNumber of cases: " + str(cases))
    print("\t\tNumber of controls: " + str(controls))
    print("\t\tSample size: " + str(sample_size))
    print("\t\tOdds ratio: " + str(odds_ratio))
    if mut > len(vcf_samples_mut) or wt > len(vcf_samples_wt):
        result = least_squares(solve_odds_ratio_function_2, initialGuess,
                               args=(len(vcf_samples_mut), len(vcf_samples_wt), float(args.allele_frequency), int(sample_size), float(odds_ratio)),
                               bounds=(0, int(sample_size)))
    else:
        result = least_squares(solve_odds_ratio_function, initialGuess,
                               args=(int(mut), int(wt), int(cases), int(controls), float(odds_ratio)),
                               bounds=(0, int(sample_size)))
    [controls_mut, cases_mut, cases_wt, controls_wt] = result.x
    # Values are rounded to integers, 0 samples are not allowed.
    # NOT DONE, as it will cause errors when only one mutated sample if found
    controls_mut = int(math.floor(controls_mut))
    cases_mut = int(math.floor(cases_mut))
    cases_wt = int(math.floor(cases_wt))
    controls_wt = int(math.floor(controls_wt))
    print("\tSolved variables:")
    print("\t\tNumber of mutated controls: " + str(round(controls_mut)))
    print("\t\tNumber of mutated cases: " + str(cases_mut))
    print("\t\tNumber of wild-type controls: " + str(controls_wt))
    print("\t\tNumber of wild-type cases: " + str(cases_wt))
    # print("\t\tCalculated odds ratio: " + str(round((cases_mut/controls_mut)/(cases_wt/controls_wt))))

    # Randomly selecting samples across repetitions
    logging.info(f"Randomly selecting samples across repetitions, creating phenotype file...")
    phen = pd.DataFrame({'samples': vcf_samples})
    for r in range(1, int(args.simu_rep)+1):
        cases_mut_samples = random.sample(vcf_samples_mut, cases_mut)
        controls_mut_samples = random.sample(set(vcf_samples_mut) - set(cases_mut_samples), controls_mut)
        cases_wt_samples = random.sample(vcf_samples_wt, cases_wt)
        controls_wt_samples = random.sample(set(vcf_samples_wt)-set(cases_wt_samples), controls_wt)
        cases_samples = list(cases_mut_samples) + list(cases_wt_samples)
        controls_samples = list(controls_mut_samples) + list(controls_wt_samples)
        print("\t\tNumber of samples selected as mutated contols: " + str(len(controls_mut_samples)))
        print("\t\tNumber of samples selected as mutated cases: " + str(len(cases_mut_samples)))
        print("\t\tNumber of samples selected as wild-type controls: " + str(len(controls_wt_samples)))
        print("\t\tNumber of samples selected as wild-type cases: " + str(len(cases_wt_samples)))
        # Assigning phenotypes to cases and controls
        phenotype = ['NA'] * len(vcf_samples)
        for idx, sample in enumerate(vcf_samples):
            if sample in cases_samples:
                phenotype[idx] = '1'
            if sample in controls_samples:
                phenotype[idx] = '0'
        # Appending phenotype list to dataframe
        phen['phenotype' + str(r)] = phenotype
    phen.to_csv(args.out_prefix + '.phen', header=True, sep='\t', index=False)
    logging.info(f"Simulated phenotype file saved as {args.out_prefix + '.phen'}")

    # Removing temporary files
    logging.info('Removing temporary files.')
    rm_command = ['rm', tmp_vcf_file]
    run_command_shell_string(' '.join(rm_command))


if __name__ == "__main__":
    _main()





















