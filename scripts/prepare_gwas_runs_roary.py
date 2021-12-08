#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
from random import choice
from string import ascii_uppercase


# ---------------------------------------------------------------------------------------------------------------------
# Developing notes
# ---------------------------------------------------------------------------------------------------------------------

# To do:
#   - Add Pyseer options that can be parsed to pyseer GWAS command


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to prepare power calculation GWAS runs across parameters for Roary. \nThis script reads a " \
                  "parameters file (--parameters_file) with 'Causal variant sampling parameters' and 'Phenotype " \
                  "simulation parameters' to create an individual GWAS bash script for each unique parameter " \
                  "combination. Only PySeer is supported. In addition, this script outputs two files: a table with " \
                  "GWAS parameters and unique identifiers for each parameter combination, and a bash script with all " \
                  "individual GWAS scripts used to execute all GWAS runs sequentially or to submitted them as jobs" \
                  " (see --job_submission_command option)"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required: Input files')
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
        "-x", "--plink_prefix", action="store", dest="plink_prefix",
        help="Prefix of PLINK input files (i.e. prefix.bed, prefix.bim and prefix.fam)",
        required=True, metavar="PLINK_PREFIX"
    )
    group.add_argument(
        "-a", "--pastml_steps_file", action="store", dest="pastml_steps_file",
        help="table with number of homoplasies per variant created by script ancestral_state_reconstruction.py",
        required=True, metavar="PASTML_STEPS"
    )
    group.add_argument(
        "-p", "--parameters_file", action="store", dest="parameters_file",
        help="File with parameters to run power calculations",
        required=True, metavar="PARAMETERS"
    )
    group = parser.add_argument_group('Optional: Pyseer Options')
    group.add_argument(
        "-m", "--pyseer_model", action="store", dest="pyseer_model",
        help="Pyseer model to run. To choose between: \"LMM\" (Linear Mixed Models) or \"ENET\" (Elastic Net). "
             "Default: LMM",
        required=False, default='LMM', metavar="MODEL"
    )
    group.add_argument(
        "-w", "--phenotype_type", action="store", dest="phenotype_type",
        help="Type of phenotype to be simulated. To choose between: \"binary\" or \"quantitative\". [Default: binary]",
        required=False, default='binary', metavar="PHEN_TYPE"
    )
    group.add_argument(
        "-d", "--similarity", action="store", dest="similarity",
        help="Pyseer input strains similarity square matrix (for --lmm)",
        required=False, metavar="SIMILARITY"
    )
    group.add_argument(
        "-u", "--cpu", action="store", dest="cpu",
        help="Processes used by Pyseer [Default: 1]",
        required=False, default=1, metavar="CPU"
    )
    group = parser.add_argument_group('Optional: executable paths')
    group.add_argument(
        "-c", "--code_directory", action="store", dest="code_directory",
        help="Full path to directory containing power calculations scripts",
        required=False, default='', metavar="CODE_DIR"
    )
    group.add_argument(
        "-s", "--pyseer_path", action="store", dest="pyseer_path",
        help="Full path to Pyseer script", required=False, default='pyseer', metavar="PYSEER_PATH"
    )
    group.add_argument(
        "-k", "--plink_path", action="store", dest="plink_path",
        help="Full path to PLINK executable", required=False, default='plink', metavar="PLINK_PATH"
    )
    group.add_argument(
        "-g", "--gcta_path", action="store", dest="gcta_path",
        help="Full path to Pyseer script", required=False, default='gcta64', metavar="GCTA_PATH"
    )
    group = parser.add_argument_group('Optional: output')
    group.add_argument(
        "-o", "--output_dir", action="store", dest="output_dir",
        help="Output directory where pipeline output files will be stored",
        required=False, default='', metavar="OUTPUT_DIR"
    )
    group.add_argument(
        "-y", "--output_prefix", action="store", dest="output_prefix",
        help="Output prefix used to name output files (output_prefix.gwas_runs.csv and output_prefix.gwas_runs.sh). "
             "Where output_prefix.gwas_runs.csv contains GWAS parameters and unique identifiers for each parameter "
             "combination (GWAS run), and output_prefix.gwas_runs.sh is the bash script used to run GWAS scripts.",
        required=False, default='output', metavar="OUTPUT_PREFIX"
    )
    group.add_argument(
        "-j", "--job_submission_command", action="store", dest="job_submission_command",
        help="Job submission command to include in output_prefix.gwas_runs.sh. Place --job_submission_command option "
             "as the last argument when calling prepare_gwas_runs.py, to avoid an 'error: unrecognized arguments' caused"
             " by dashes in the argument string. Using the string 'job_id' in the argument provided will allow to "
             "assign a unique identifier (combination_id) to each job.",
        required=False, default='', nargs=argparse.REMAINDER, metavar="JOB_SUBMISSION"
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


def check_path_exists(my_dir):
    if not os.path.exists(my_dir):
        logging.error(f'Directory {my_dir} not found!')
        sys.exit(-1)


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


def check_list_not_empty(list, list_name):
    """
    Function used to check a list is not empty
    :param list: list to check length for
    :param list_name: name of variable
    :return:
    """
    if len(list) == 0:
        logging.error(f"{list_name} is not found in parameters file.")
        sys.exit(-1)


def check_value_not_zero(value, value_name):
    """
    Function to check value (int or float) is not zero
    :param value:
    :param value_name:
    :return:
    """
    if float(value) == 0:
        logging.error(f"{value_name} not found in parameters file.")
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

    # Making sure pyseer_model choosen option is correct
    if args.pyseer_model != "LMM":
        if args.pyseer_model != "ENET":
            logging.error(f'{args.pyseer_model} must be either \"LMM\" or \"ENET\"')
            sys.exit(-1)

    # Making sure phenotype type chosen option is correct
    if args.phenotype_type != "binary":
        if args.phenotype_type != "quantitative":
            logging.error(f'--phenotype_type must be either \"binary\" or \"quantitative\"')
            sys.exit(-1)

    # Making sure dependency executables exist
    logging.info('Making sure dependencies exist...')
    dependencies = [args.gcta_path, args.plink_path]
    for dependency in dependencies:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)
    logging.info('All required dependencies found.')

    # Making sure directories exist
    code_directory = ''
    if args.code_directory != '':
        code_directory = os.path.abspath(args.code_directory)
        check_path_exists(code_directory)
        code_directory += '/'
    output_dir = ''
    if args.output_dir != '':
        output_dir = os.path.abspath(args.output_dir)
        check_path_exists(output_dir)
        output_dir += '/'

    # Making sure python scripts exist
    logging.info('Making sure scripts exist...')
    scripts = [code_directory + 'sample_casual_variants_from_roary.py',
               code_directory + 'simulate_phenotype_using_gcta.py',
               code_directory + 'simulate_binary_phenotype_roary.py',
               code_directory + 'plink_to_pyseer_phenotype_file.py', args.pyseer_path]
    for script in scripts:
        check_dependency(script)
    logging.info('All required scripts found.')

    # Making sure input files exist
    logging.info('Making sure input files exist...')
    input_files = [args.roary_table, args.parameters_file, args.pastml_steps_file]
    if args.pyseer_model == "LMM":
        input_files.append(args.similarity)
    for input_file in input_files:
        check_file_exists(input_file)
    logging.info('All required input files found.')

    # Making sure PLINK input files exist
    logging.info('Making sure PLINK input files exist...')
    input_files = [args.plink_prefix + '.bed', args.plink_prefix + '.bim', args.plink_prefix + '.fam']
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found. Run roary_to_plink_files.py first.')
            sys.exit(-1)
    logging.info('All required PLINK input files found.')

    # Saving parameters from parameters file, and checking values format
    allele_frequency_from = list()
    allele_frequency_to = list()
    homoplasy_steps_from = list()
    homoplasy_steps_to = list()
    number_causal_variants = list()
    sampling_repetitions = 0
    sample_size = list()
    case_control_ratio = "NA"
    heritability = list()
    effect_size = list()
    prevalence = 0
    simulation_repetitions = 0
    with open(args.parameters_file, 'r') as input_file:
        for line in input_file:
            if line.startswith('allele_frequency_from'):
                allele_frequency_from = line.strip().split(' ')[1].split(',')
                for a in allele_frequency_from:
                    check_positive_float(a, 'allele_frequency_from', 0, 1)
            if line.startswith('allele_frequency_to'):
                allele_frequency_to = line.strip().split(' ')[1].split(',')
                for a in allele_frequency_to:
                    check_positive_float(a, 'allele_frequency_to', 0, 1)
            if line.startswith('effect_size'):
                effect_size = line.strip().split(' ')[1].split(',')
                for e in effect_size:
                    check_positive_float(e, 'effect_size', 0, None)
            if line.startswith('homoplasy_steps_from'):
                homoplasy_steps_from = line.strip().split(' ')[1].split(',')
                for h in homoplasy_steps_from:
                    check_positive_integer(h, 'homoplasy_steps_from', None, None)
            if line.startswith('homoplasy_steps_to'):
                homoplasy_steps_to = line.strip().split(' ')[1].split(',')
                for h in homoplasy_steps_to:
                    check_positive_integer(h, 'homoplasy_steps_to', None, None)
            if line.startswith('number_causal_variants'):
                number_causal_variants = line.strip().split(' ')[1].split(',')
                for h in number_causal_variants:
                    check_positive_integer(h, 'number_causal_variants', None, None)
            if line.startswith('sampling_repetitions'):
                sampling_repetitions = line.strip().split(' ')[1]
                check_positive_integer(sampling_repetitions, 'sampling_repetitions', None, None)
            if line.startswith('sample_size'):
                sample_size = line.strip().split(' ')[1].split(',')
                for s in sample_size:
                    check_positive_integer(s, 'sample_size', 0, None)
            if line.startswith('case_control_ratio'):
                case_control_ratio = line.strip().split(' ')[1]
                check_positive_float(case_control_ratio, 'case_control_ratio', 0, 1)
            if line.startswith('heritability'):
                heritability = line.strip().split(' ')[1].split(',')
                for h in heritability:
                    check_positive_float(h, 'heritability', 0, 1)
                # heritability = line.strip().split(' ')[1]
                # check_positive_float(heritability, 'heritability', 0, 1)
            if line.startswith('prevalence'):
                prevalence = line.strip().split(' ')[1]
                check_positive_float(prevalence, 'prevalence', 0, 1)
            if line.startswith('simulation_repetitions'):
                simulation_repetitions = line.strip().split(' ')[1]
                check_positive_integer(simulation_repetitions, 'simulation_repetitions', None, None)

    # check all required parameters found in parameters file
    check_list_not_empty(allele_frequency_from, 'allele_frequency')
    check_list_not_empty(allele_frequency_to, 'allele_frequency')
    if len(allele_frequency_from) != len(allele_frequency_to):
        logging.error(f'Make sure allele_frequency_from and allele_frequency_to have the same length in {args.parameters_file}.')
        sys.exit(-1)
    check_list_not_empty(effect_size, 'effect_size')
    check_list_not_empty(homoplasy_steps_from, 'homoplasy_steps')
    check_list_not_empty(homoplasy_steps_to, 'homoplasy_steps')
    if len(homoplasy_steps_from) != len(homoplasy_steps_to):
        logging.error(f'Make sure homoplasy_steps_from and homoplasy_steps_to have the same length in {args.parameters_file}.')
        sys.exit(-1)
    check_list_not_empty(number_causal_variants, 'number_causal_variants')
    check_list_not_empty(sample_size, 'sample_size')
    check_value_not_zero(sampling_repetitions, 'sampling_repetitions')
    check_list_not_empty(heritability, 'heritability')
    check_value_not_zero(prevalence, 'prevalence')
    check_value_not_zero(simulation_repetitions, 'simulation_repetitions')
    if args.phenotype_type == "binary":
        check_value_not_zero(case_control_ratio, 'case_control_ratio')

    # Make sure both heritability and effect_size vectors are not greater than 1.
    # If a range of heritability values is selected, then a single effect_size must be chosen, and vice versa.
    if len(heritability) > 1:
        if len(effect_size) > 1:
            logging.error(f'Range of values selected for both heritability and effect_size parameters. '
                          f'Range only supported for one of them.')
            sys.exit(-1)

    # Printing chosen parameters
    print("Causal variant sampling parameters: ")
    print("\tAllele frequency from: " + ' '.join(allele_frequency_from))
    print("\tAllele frequency to: " + ' '.join(allele_frequency_to))
    print("\tEffect sizes: " + ' '.join(effect_size))
    print("\tHomoplasy steps from: " + ' '.join(homoplasy_steps_from))
    print("\tHomoplasy steps to: " + ' '.join(homoplasy_steps_to))
    print("\tNumber causal variants: " + ' '.join(number_causal_variants))
    print("\tSampling repetitions: " + str(sampling_repetitions))
    print("Phenotype simulation parameters: ")
    print("\tSample sizes: " + ' '.join(sample_size))
    print("\tCase control ratio: " + str(case_control_ratio))
    print("\tHeritability: " + ' '.join(heritability))
    print("\tPrevalence: " + str(prevalence))
    print("\tSimulation repetitions: " + str(simulation_repetitions))

    # Total number of combinations
    num_combinations = len(allele_frequency_from) * len(effect_size) * len(homoplasy_steps_from)\
                       * len(number_causal_variants) * int(sampling_repetitions) * len(sample_size)\
                       * int(simulation_repetitions)
    print("Total number unique parameter combinations (GWAS runs): " + str(num_combinations))

    # Creating GWAS bash script for each unique parameter combination, and output file
    output_table_header = 'combination_id\tallele_frequency\teffect_size\thomoplasy_steps\tnumber_causal_variants\t' \
                          'sample_repetition\tsample_size\tcase_control_ratio\theritability\tprevalence\t' \
                          'simulation_repetition\n'
    output_table = open(args.output_prefix + '.gwas_runs.csv', 'w')
    output_table.write(output_table_header)
    output_bash = open(args.output_prefix + '.gwas_runs.sh', 'w')
    output_id_bash = open(args.output_prefix + '.ids_gwas_runs.csv', 'w')
    for idx, aff in enumerate(allele_frequency_from):
        aft = allele_frequency_to[idx]
        af = str(aff) + ',' + str(aft)
        for es in effect_size:
            for idx, hsf in enumerate(homoplasy_steps_from):
                hst = homoplasy_steps_to[idx]
                hs = str(hsf) + ',' + str(hst)
                for het in heritability:
                    for ncv in number_causal_variants:
                        for sr in list(range(1, int(sampling_repetitions)+1, 1)):
                            for ss in sample_size:
                                # Number of cases and controls
                                nca = "NA"
                                nco = "NA"
                                if args.phenotype_type == "binary":
                                    nca = int(int(ss) * float(case_control_ratio))
                                    nco = int(ss) - int(nca)
                                    check_positive_integer(nca, 'number of cases', None, None)
                                    check_positive_integer(nco, 'number of controls', None, None)

                                # creating combination id
                                combination_id = ''.join(choice(ascii_uppercase) for i in range(12))

                                # files common to all phenotype replicates
                                causal_variants = output_dir + combination_id + '.causal_variants.txt'
                                output_prefix = output_dir + combination_id
                                plink_phen_file = output_dir + combination_id + '.phen'
                                pyseer_phen_file = output_dir + combination_id + '.pyseer.phen'
                                gwas_script_file = output_dir + combination_id + '.gwas_run.sh'
                                sub_samples = output_dir + combination_id + '.subsample.txt'

                                # Creating one GWAS script per parameter combination (and for all phenotype replicates)
                                gwas_script = open(gwas_script_file, 'w')
                                command_line = '#!/bin/bash\n#GWAS run: ' + combination_id + '\n'
                                command_line += 'if [ ! -f ' + causal_variants + ' ]\n'
                                command_line += 'then\n'
                                command_line += code_directory + 'sample_casual_variants_from_roary.py' \
                                                + ' -v ' + args.roary_table + ' -n ' + str(ncv) + ' -a ' + str(af) \
                                                + ' -s ' + str(hs) + ' -e ' + str(es) + ' -o ' + causal_variants \
                                                + ' -p ' + args.pastml_steps_file + ' -f ' + args.input_format + '\n'
                                command_line += 'fi\n'
                                command_line += 'if [ -f ' + causal_variants + ' ]\n'
                                command_line += 'then\n'
                                command_line += 'if [ ! -f ' + pyseer_phen_file + ' ]\n'
                                command_line += 'then\n'
                                # If there is a range of heritability values use simulate_phenotype_using_gcta.py,
                                # else use simulate_binary_phenotype_vcf.py
                                if args.phenotype_type == "binary":
                                    if len(heritability) > 1:
                                        command_line += code_directory + 'simulate_phenotype_using_gcta.py' \
                                                        + ' --bfile ' + args.plink_prefix \
                                                        + ' --out ' + output_prefix \
                                                        + ' --plink_path ' + args.plink_path \
                                                        + ' --gcta_path ' + args.gcta_path \
                                                        + ' --causal-loci ' + causal_variants \
                                                        + ' --simu-cc ' + str(nca) + ',' + str(nco) \
                                                        + ' --simu-hsq ' + str(het) \
                                                        + ' --simu-k ' + str(prevalence) \
                                                        + ' --simu-rep ' + str(simulation_repetitions) + '\n'
                                        command_line += code_directory + 'plink_to_pyseer_phenotype_file.py' \
                                                        + ' --plink_phenotype_file ' + plink_phen_file \
                                                        + ' --pyseer_phenotype_file ' + pyseer_phen_file + '\n'
                                    else:
                                        command_line += code_directory + 'simulate_binary_phenotype_roary.py' \
                                                        + ' --input_roary ' + args.roary_table \
                                                        + ' --input_format ' + args.input_format \
                                                        + ' --out_prefix ' + output_prefix \
                                                        + ' --causal-loci ' + causal_variants \
                                                        + ' --simu-cc ' + str(nca) + ',' + str(nco) \
                                                        + ' --simu-rep ' + str(simulation_repetitions) + '\n'
                                        command_line += 'mv ' + output_prefix + '.phen ' + pyseer_phen_file + '\n'
                                if args.phenotype_type == "quantitative":
                                    command_line += code_directory + 'subsample_pangenome_samples.py' \
                                                    + ' --roary_table ' + args.roary_table \
                                                    + ' --input_format ' + args.input_format \
                                                    + ' --causal-loci ' + causal_variants \
                                                    + ' --subset_size ' + str(ss) \
                                                    + ' --output_samples ' + sub_samples + '\n'
                                    command_line += code_directory + 'simulate_phenotype_using_gcta.py' \
                                                    + ' --bfile ' + args.plink_prefix \
                                                    + ' --out ' + output_prefix \
                                                    + ' --plink_path ' + args.plink_path \
                                                    + ' --gcta_path ' + args.gcta_path \
                                                    + ' --causal-loci ' + causal_variants \
                                                    + ' --simu-qt --keep ' + sub_samples \
                                                    + ' --simu-hsq ' + str(het) \
                                                    + ' --simu-rep ' + str(simulation_repetitions) + '\n'
                                    command_line += code_directory + 'plink_to_pyseer_phenotype_file.py' \
                                                    + ' --plink_phenotype_file ' + plink_phen_file \
                                                    + ' --pyseer_phenotype_file ' + pyseer_phen_file + '\n'
                                command_line += 'fi\n'
                                # For each phenotype replicate
                                for ir in list(range(1, int(simulation_repetitions) + 1, 1)):
                                    replicate_id = combination_id + '_' + str(ir)
                                    # saving parameter combination
                                    newline_table = '\t'.join([str(replicate_id), str(af), str(es), str(hs), str(ncv),
                                                               str(sr), str(ss), str(case_control_ratio), str(het),
                                                               str(prevalence), str(ir)])
                                    output_table.write(newline_table + '\n')
                                    # creating one GWAS script per parameter combination
                                    # Phenotype file column to use in pyseer_phen_file
                                    phenotype_column = 'phenotype' + str(ir)
                                    pyseer_output = output_dir + replicate_id + '.pyseer_output.csv'
                                    command_line += 'if [ ! -f ' + pyseer_output + ' ]\n'
                                    command_line += 'then\n'
                                    command_line += args.pyseer_path
                                    if args.pyseer_model == "LMM":
                                        command_line += ' --lmm ' + ' --similarity ' + args.similarity
                                    if args.pyseer_model == "ENET":
                                        command_line += ' --wg enet ' + ' --alpha 0.01'
                                    if args.phenotype_type == "quantitative":
                                        command_line += ' --continuous'
                                    command_line += ' --pres ' + args.roary_table \
                                                    + ' --phenotypes ' + pyseer_phen_file \
                                                    + ' --phenotype-column ' + str(phenotype_column) \
                                                    + ' --min-af 0.001 --print-filtered --cpu ' + str(args.cpu) \
                                                    + ' > ' + pyseer_output + '\n'
                                    # command_line += args.pyseer_path + ' --lmm ' \
                                    #                 + ' --pres ' + args.roary_table \
                                    #                 + ' --phenotypes ' + pyseer_phen_file \
                                    #                 + ' --phenotype-column ' + str(phenotype_column) \
                                    #                 + ' --similarity ' + args.similarity \
                                    #                 + ' --min-af 0.001 --print-filtered ' \
                                    #                 + ' > ' + pyseer_output + '\n'
                                    command_line += 'fi\n'
                                # command_line += 'rm ' + causal_variants + '\n'
                                # command_line += 'rm ' + plink_phen_file + '\n'
                                # command_line += 'rm ' + pyseer_phen_file + '\n'
                                # command_line += 'rm ' + output_dir + combination_id + '.log' + '\n'
                                # command_line += 'rm ' + output_dir + combination_id + '.par' + '\n'
                                command_line += 'fi\n'
                                gwas_script.write(command_line)
                                gwas_script.close()
                                job_submission_command = ' '.join(args.job_submission_command)
                                # replacing string job_id with combination_id
                                job_submission_command = job_submission_command.replace('job_id', str(combination_id))
                                output_bash.write(job_submission_command + ' "bash ' + gwas_script_file + '"\n')
                                output_id_bash.write(str(combination_id) + '\t' + gwas_script_file + '\n')
                                print(command_line)
    output_table.close()
    output_bash.close()
    output_id_bash.close()


if __name__ == "__main__":
    _main()