#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf
import subprocess
import math

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Wrapper script used to simulate case-control phenotypes using GCTA.\n"\
                  "A file of sampled causal variants obtained from sample_casual_variants_from_vcf.py or " \
                  "sample_casual_variants_from_roary.py is needed as input. Input variants need to be formatted in " \
                  "PLINK format (e.g. by using vcf_to_plink_files.py or roary_to_plink_files.py)" \
                  "NOTE: Simulation of quantitative phenotypes is currently not supported, use GCTA tool directly.\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('I/O arguments')
    group.add_argument(
        "-b", "--bfile", action="store", dest="bfile",
        help="Input prefix of PLINK formatted files",
        required=True, metavar="BFILE"
    )
    group.add_argument(
        "-o", "--out_prefix", action="store", dest="out_prefix",
        help="output prefix used for GTCA output files",
        required=True, metavar="OUTPUT_PREFIX"
    )
    group.add_argument(
        "-g", "--gcta_path", action="store", dest="gcta_path",
        help="path to GTCA software tool", required=False, default='gcta64', metavar="GCTA_PATH"
    )
    group.add_argument(
        "-p", "--plink_path", action="store", dest="plink_path",
        help="Full path to PLINK executable", required=False, default='plink', metavar="PLINK_PATH"
    )
    group.add_argument(
        "-c", "--causal-loci", action="store", dest="causal_loci",
        help="List of SNPs as causal variants.",
        required=True, metavar="CAUSAL_LOCI"
    )
    group.add_argument(
        "-u", "--effect_size_units", action="store", dest="effect_size_units",
        help="Units of effect sizes in causal variants file: choose between 'beta' or 'oddsratio'. Default: 'beta'",
        required=False, default='beta', metavar="EF_UNITS"
    )
    group = parser.add_argument_group('GCTA arguments')
    group.add_argument(
        "-x", "--simu-cc", action="store", dest="simu_cc",
        help="Number of cases and controls for binary phenotypes, separated by comma (e.g. 100,200)",
        required=False, metavar="SIMU_CC"
    )
    group.add_argument(
        "-q", "--simu-qt", action="store_true", dest="simu_qt",
        help="Simulate a quantitative trait.",
        required=False
    )
    group.add_argument(
        "-e", "--keep", action="store", dest="keep",
        help="Subset of individuals to keep for quantitative trait simulation",
        required=False, metavar="KEEP"
    )
    group.add_argument(
        "-s", "--simu-hsq", action="store", dest="simu_hsq",
        help="Specify the heritability (or heritability of liability), e.g. 0.8. "
             "The default value is 0.1 if this option is not specified.",
        required=False, default="0.1", metavar="SIMU_HSQ"
    )
    group.add_argument(
        "-k", "--simu-k", action="store", dest="simu_k",
        help="Specify the disease prevalence, e.g. 0.01. The default value is 0.1 if this option is not specified.",
        required=False, default="0.01", metavar="SIMU_K"
    )
    group.add_argument(
        "-r", "--simu-rep", action="store", dest="simu_rep",
        help="Number of simulation replicates. The default value is 1 if this option is not specified.",
        required=False, metavar="SIMU_REP"
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


def run_command_shell_string(command_line_string):
    """
    This function executes a command line, check for execution errors and but does not return stdout
    This is to be used when the stdout is not needed
    Note: shell=True needs to be set if I/O redirection operators are to be used (e.g. >) in the command line,
    otherwise they will have no special meaning, they are treated as ordinary arguments
    Note: if shell=True is used then the command line must be provided as a string, not a list
    :param command_line_string: it must be a string not a list
    """
    # print('Running\n' + command_line_string)
    try:
        subprocess.run(command_line_string,
                       check=True,
                       shell=True,
                       )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)


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
    dependencies = [args.gcta_path, args.plink_path]
    for dependency in dependencies:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)

    # Making sure input files exist
    input_files = [args.bfile + '.bed', args.bfile + '.bim', args.bfile + '.fam', args.causal_loci]
    if args.keep is not None:
        input_files.append(args.keep)
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Making sure chosen parameters are not incompatible
    if args.simu_cc is not None:
        if args.simu_qt is True:
            logging.error(f'Choose either --simu-qt or --simu-cc, but not both.')
            sys.exit(-1)
    if args.simu_qt is True:
        if args.keep is None:
            logging.error(f'Argument --keep must be provided if --simu-qt chosen.')
            sys.exit(-1)

    # Parsing parameters
    (num_cases, num_controls, prevalence) = (None, None, None)
    if args.simu_cc is not None:
        (num_cases, num_controls) = args.simu_cc.split(',')
        prevalence = args.simu_k
    heritability = args.simu_hsq
    replicates = args.simu_rep
    ef_units = args.effect_size_units.strip()

    # Making sure input parameters meet expected values
    check_positive_integer(replicates, 'number of replicates (--simu-rep)', None, None)
    check_positive_float(heritability, 'heritability (--simu-hsq)', 0, 1)
    if args.simu_cc is not None:
        check_positive_integer(num_cases, 'number of cases (--simu-cc)', None, None)
        check_positive_integer(num_controls, 'number of controls (--simu-cc)', None, None)
        check_positive_float(prevalence, 'disease prevalence (--simu-k)', 0, 1)
    if ef_units != 'beta':
        if ef_units != 'oddsratio':
            logging.error(f'--effect_size_units specified {ef_units} must be either \'beta\' or \'oddsratio\'.')
            sys.exit(-1)

    # Making sure number of cases and controls do not exceed sample size
    plink_samples = dict()
    with open(args.bfile + '.fam') as input_file:
        for line in input_file:
            sample_id = line.strip().split(' ')[0]
            plink_samples[sample_id] = ''
    sample_size = len(plink_samples)

    # If binary phenotype chosen, make sure number of cases and controls do not exceed sample size
    if args.simu_cc is not None:
        if (int(num_cases)+int(num_controls)) > int(sample_size):
            logging.error(f'Sum of cases and controls {str(int(num_cases)+int(num_controls))} larger than sample size '
                          f'{str(sample_size)}')
            sys.exit(-1)
    # If quantitative phenotype chosen, make sure samples in --keep are found in PLINK files
    keep_samples = dict()
    if args.simu_qt is True:
        with open(args.keep) as input_file:
            for line in input_file:
                sample_id = line.strip()
                keep_samples[sample_id] = ''
                if sample_id not in plink_samples:
                    logging.error(f'individual {sample_id} in {args.keep} not found in PLINK input files.')
                    sys.exit(-1)
        if len(keep_samples) > len(plink_samples):
            logging.error(f'Number of individuals in {args.keep} exceeds those in PLINK input files.')
            sys.exit(-1)

    # Reading all variant ids from .bim file
    all_variants = dict()
    with open(args.bfile + '.bim', 'r') as input_file:
        for line in input_file:
            var_id = line.strip().split('\t')[1]
            all_variants[var_id] = 0

    # Reading causal variants. Changing effect size units from odds ratio (specified in parameters file) to beta
    # (expected by GCTA) if chosen (--effect_size_units)
    causal_variants = dict()
    causal_variants_new_text = ""
    with open(args.causal_loci, 'r') as input_file:
        for line in input_file:
            var_id = line.strip().split('\t')[0]
            causal_variants[var_id] = 0
            # Changing effect size units from odds ratio (specified in parameters file) to beta (expected by GCTA)
            odds_ratio = line.strip().split('\t')[1]
            beta = math.log(float(odds_ratio))
            causal_variants_new_text += var_id + '\t' + str(beta) + '\n'
    if ef_units == 'oddsratio':
        output_causal_loci = open(args.causal_loci, 'w')
        output_causal_loci.write(causal_variants_new_text)
        output_causal_loci.close()

    # Making sure causal variants exist in PLINK files
    for var_id in causal_variants:
        if var_id not in all_variants:
            logging.error(f'Causal variant \'{var_id}\' in {args.causal_loci} not found in PLINK file {args.bfile}.bim\n')
            sys.exit(-1)

    # Running PLINK to generate a list of minor allele frequencies (MAF) for each SNP
    # Making sure causal variant MAF does not exceed disease prevalence (applies to case-control phenotype)
    if args.simu_cc is not None:
        logging.info('Running PLINK to calculate minor allele frequencies.')
        plink_command = [args.plink_path, '--bfile', args.bfile, '--freq', '--out', args.bfile]
        if not os.path.exists(args.bfile + '.frq'):
            run_command_shell_string(' '.join(plink_command))
        maf_variants = dict()
        if os.path.exists(args.bfile + '.frq'):
            with open(args.bfile + '.frq', 'r') as input_file:
                for line in input_file:
                    fields = line.strip().split()
                    var_id = fields[1]
                    allele_frequency = fields[4]
                    maf_variants[var_id] = allele_frequency
        else:
            logging.error('Output .frq PLINK file does not exist. Something went wrong.')
            sys.exit(-1)

        for var_id in causal_variants:
            if float(maf_variants[var_id]) > float(args.simu_k):
                logging.error(f'Causal variant \'{var_id}\' in {args.causal_loci} with MAF greater '
                              f'({maf_variants[var_id]}) than selected disease prevalence {args.simu_k} (--simu-k)\n')
                sys.exit(-1)

    # Printing chosen parameters
    print("Chosen GCTA parameters: ")
    print("\tNumber of individuals in PLINK files: " + str(sample_size))
    print("\tNumber of variants in PLINK files: " + str(len(all_variants)))
    print("\tNumber of causal variants: " + str(len(causal_variants)))
    print("\tHeritability: " + str(heritability))
    print("\tNumber of replicates: " + str(replicates))
    if args.simu_cc is not None:
        print("\tSimulation of binary phenotype chosen: ")
        print("\t\tNumber of cases/controls: " + str(num_cases) + '/' + str(num_controls))
        print("\t\tDisease prevalence: " + str(prevalence))
    if args.simu_qt is True:
        print("\tSimulation of quantitative phenotype chosen: ")
        print("\t\tNumber of selected samples: " + str(len(keep_samples)))

    # Running GCTA
    logging.info('Running GCTA to simulate phenotypes.')
    if args.simu_cc is not None:
        gtca_command = [args.gcta_path, '--bfile', args.bfile, '--simu-cc', num_cases, num_controls,
                        '--simu-causal-loci', args.causal_loci, '--simu-hsq', heritability, '--simu-k', prevalence,
                        '--simu-rep', replicates, '--out', args.out_prefix]
        run_command_shell_string(' '.join(gtca_command))
    if args.simu_qt is True:
        gtca_command = [args.gcta_path, '--bfile', args.bfile, '--simu-qt', '--keep', args.keep,
                        '--simu-causal-loci', args.causal_loci, '--simu-hsq', heritability,
                        '--simu-rep', replicates, '--out', args.out_prefix]
        run_command_shell_string(' '.join(gtca_command))


if __name__ == "__main__":
    _main()