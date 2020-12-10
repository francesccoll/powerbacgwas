#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf
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
    description = "Script to prepare power calculation GWAS runs for sub-sampling approach (VCF). \nThis script reads a " \
                  "parameters file (--parameters_file) with 'Sub-sampling phenotype simulation parameters' " \
                  "to create an individual GWAS bash script for each unique parameter combination. " \
                  "Only PySeer is supported. In addition, this script outputs two files: a table with " \
                  "GWAS parameters and unique identifiers for each parameter combination; and a bash script with all " \
                  "individual GWAS scripts used to run them all sequentially or to submit them as jobs" \
                  " (see --job_submission_command option)"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required: Input files')
    group.add_argument(
        "-v", "--input_vcf", action="store", dest="input_vcf",
        help="multi-sample VCF file (zipped and indexed)",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-d", "--similarity", action="store", dest="similarity",
        help="Pyseer input strains similarity square matrix (for --lmm)",
        required=True, metavar="SIMILARITY"
    )
    group.add_argument(
        "-p", "--parameters_file", action="store", dest="parameters_file",
        help="File with parameters to run power calculations",
        required=True, metavar="PARAMETERS"
    )
    group.add_argument(
        "-l", "--causal-loci", action="store", dest="causal_loci",
        help="GCTA-formatted file with causal loci causing the phenotype specified --phenotype",
        required=True, metavar="CAUSAL_LOCI"
    )
    group.add_argument(
        "-t", "--phenotype", action="store", dest="phenotype",
        help="Pyseer-formatted phenotype file caused by the causal loci specified in --causal-loci",
        required=True, metavar="PHENOTYPE"
    )
    group = parser.add_argument_group('Optional: burden testing')
    group.add_argument(
        "-b", "--burden", action="store", dest="burden",
        help="VCF regions to group variants by for burden testing. The regions file can be used to specify gene "
             "coordinates and must contain one region per line, with their name and the bcftools style region "
             "co-ordinates delimited by tab (locus_name\tchr_name:start-end).",
        required=False, metavar="REGIONS_FILE"
    )
    group = parser.add_argument_group('Optional: executable paths')
    group.add_argument(
        "-c", "--code_directory", action="store", dest="code_directory",
        help="Full path to directory containing power calculations scripts",
        required=False, default='./', metavar="CODE_DIR"
    )
    group.add_argument(
        "-s", "--pyseer_path", action="store", dest="pyseer_path",
        help="Full path to Pyseer script", required=False, default='pyseer', metavar="PYSEER_PATH"
    )
    group = parser.add_argument_group('Optional: output')
    group.add_argument(
        "-o", "--output_dir", action="store", dest="output_dir",
        help="Output directory where pipeline output files will be stored",
        required=False, default='./', metavar="OUTPUT_DIR"
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
             "as the last argument when calling prepare_gwas_runs_subsampling.py, to avoid an 'error: unrecognized arguments' caused"
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


def var_not_zero(variable):
    variable_return = variable
    if variable == 0:
        variable_return = 0.01
    return variable_return


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

    # Making sure input directories exist
    code_directory = os.path.abspath(args.code_directory)
    output_dir = os.path.abspath(args.output_dir)
    check_path_exists(code_directory)
    check_path_exists(output_dir)
    code_directory += '/'
    output_dir += '/'

    # Making sure python scripts exist
    logging.info('Making sure scripts exist...')
    scripts = [code_directory + 'simulate_binary_phenotype_vcf.py', args.pyseer_path]
    for script in scripts:
        check_file_exists(script)
    logging.info('All required scripts found.')

    # Making sure input files exist
    logging.info('Making sure input files exist...')
    input_files = [args.input_vcf, args.similarity, args.parameters_file, args.causal_loci, args.phenotype]
    if args.burden is not None:
        input_files.append(args.burden)
    for input_file in input_files:
        check_file_exists(input_file)
    logging.info('All required input files found.')

    # Saving parameters from parameters file, and checking values format
    allele_frequency = list()
    sample_size = list()
    case_control_ratio = 0
    effect_size = list()
    simulation_repetitions = 0
    with open(args.parameters_file, 'r') as input_file:
        for line in input_file:
            if line.startswith('allele_frequency'):
                allele_frequency = line.strip().split(' ')[1].split(',')
                for a in allele_frequency:
                    check_positive_float(a, 'allele_frequency', 0, 1)
            if line.startswith('effect_size'):
                effect_size = line.strip().split(' ')[1].split(',')
                for e in effect_size:
                    check_positive_float(e, 'effect_size', 0, None)
            if line.startswith('sample_size'):
                sample_size = line.strip().split(' ')[1].split(',')
                for s in sample_size:
                    check_positive_integer(s, 'sample_size', 0, None)
            if line.startswith('case_control_ratio'):
                case_control_ratio = line.strip().split(' ')[1]
                check_positive_float(case_control_ratio, 'case_control_ratio', 0, 1)
            if line.startswith('simulation_repetitions'):
                simulation_repetitions = line.strip().split(' ')[1]
                check_positive_integer(simulation_repetitions, 'simulation_repetitions', None, None)

    # check all required parameters found in parameters file
    check_list_not_empty(allele_frequency, 'allele_frequency')
    check_list_not_empty(effect_size, 'effect_size')
    check_list_not_empty(sample_size, 'sample_size')
    check_value_not_zero(case_control_ratio, 'case_control_ratio')
    check_value_not_zero(simulation_repetitions, 'simulation_repetitions')

    # Printing chosen parameters
    print("Sub-sampling phenotype simulation parameters: ")
    print("\tAllele frequency : " + ' '.join(allele_frequency))
    print("\tEffect sizes: " + ' '.join(effect_size))
    print("\tSample sizes: " + ' '.join(sample_size))
    print("\tCase control ratio: " + str(case_control_ratio))
    print("\tSimulation repetitions: " + str(simulation_repetitions))

    # Causal variant allele frequency needs to be calculated
    # Reading causal variant file. And making sure only one variant exists.
    logging.info(f"Reading causal variant file {args.causal_loci}...")
    causal_variants = dict()
    file = open(args.causal_loci, "r")
    for line in file:
        var_id = line.strip().split('\t')[0]
        causal_variants[var_id] = 0

    # Reading VCF to extract samples ids
    logging.info(f"Reading VCF file {args.input_vcf} to extract samples with and without causal variant...")
    vcf_reader = get_vcf_reader(args.input_vcf)
    vcf_samples = vcf_reader.samples
    print("\tNumber of samples extracted from VCF: " + str(len(vcf_samples)))

    # Creating uncompressed temporary VCF file
    tmp_vcf_file = args.output_prefix + '_tmp.vcf'
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

    # Removing temporary VCF file
    run_command_shell_string(' '.join(['rm', tmp_vcf_file]))

    # Calculating causal variant allele frequency
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
    allele_frequency_obs = len(vcf_samples_mut)/(len(vcf_samples_mut)+len(vcf_samples_wt))
    print("\t\tObserved allele frequency: " + str(allele_frequency_obs))

    # Reading phenotype file
    logging.info(f"Calculating causal variant odds ratio...")
    file = open(args.phenotype, "r")
    samples_cases = []  # list to cases extracted from phenotype file
    samples_controls = []  # list to controls extracted from phenotype file
    phenotype_samples = []
    for line in file:
        [sample, phenotype] = line.strip().split('\t')
        phenotype_samples.append(sample)
        if phenotype == '0':
            samples_controls.append(sample)
        if phenotype == '1':
            samples_cases.append(sample)
    if len(samples_cases) == 0:
        logging.error(f'Number of cases from phenotype file {args.phenotype} is zero. Check format and coding.')
        sys.exit(-1)
    if len(samples_controls) == 0:
        logging.error(f'Number of controls from phenotype file {args.phenotype} is zero. Check format and coding.')
        sys.exit(-1)
    print("\t\tNumber of cases: " + str(len(samples_cases)))
    print("\t\tNumber of controls: " + str(len(samples_controls)))
    print("\t\tNumber of phenotyped samples: " + str((len(samples_cases)+len(samples_controls))))

    # Calculating odds ratio
    controls_mut = len(set(samples_controls).intersection(set(vcf_samples_mut)))
    cases_mut = len(set(samples_cases).intersection(set(vcf_samples_mut)))
    cases_wt = len(set(samples_cases).intersection(set(vcf_samples_wt)))
    controls_wt = len(set(samples_controls).intersection(set(vcf_samples_wt)))
    controls_mut = var_not_zero(controls_mut)
    cases_mut = var_not_zero(cases_mut)
    cases_wt = var_not_zero(cases_wt)
    controls_wt = var_not_zero(controls_wt)
    odds_ratio_obs = ((cases_mut / controls_mut) / (cases_wt / controls_wt))
    print("\t\tCausal variant odds ratio: " + str(odds_ratio_obs))

    # Checks
    # Making sure highest OR chosen in parameters file not higher than observed
    # Making sure highest allele frequency chosen in parameters file not higher than observed
    # Making sure highest sample size in parameters file not higher than observed
    # Making sure samples in both VCF and phenotype files match
    logging.info(f"Making a few checks before continuing...")
    if float(allele_frequency[-1]) > float(allele_frequency_obs):
        logging.error(f'Observed causal variant allele frequency {allele_frequency_obs} higher than that chosen in '
                      f'{args.parameters_file} ({allele_frequency[-1]}).')
        sys.exit(-1)
    if float(effect_size[-1]) > float(odds_ratio_obs):
        logging.error(f'Observed causal variant odds ratio {odds_ratio_obs} higher than that chosen in '
                      f'{args.parameters_file} ({effect_size[-1]}).')
        sys.exit(-1)
    matched_samples = set(vcf_samples).intersection(set(phenotype_samples))
    if len(matched_samples) != len(vcf_samples):
        logging.error(f'Samples in phenotype file {args.phenotype} do not match those in VCF file {args.input_vcf}.')
        sys.exit(-1)
    if float(sample_size[-1]) > len(vcf_samples):
        logging.error(f'Observed sample size {len(vcf_samples)} higher than that chosen in {args.parameters_file} '
                      f'({sample_size[-1]}).')
        sys.exit(-1)
    logging.info(f"Making a few checks before continuing. DONE.")

    # Total number of combinations
    num_combinations = len(allele_frequency) * len(effect_size) * len(sample_size) * int(simulation_repetitions)
    print("Total number unique parameter combinations (GWAS runs): " + str(num_combinations))

    # Creating GWAS bash script for each unique parameter combination, and output file
    output_table_header = 'combination_id\tallele_frequency\teffect_size\tsample_size\tcase_control_ratio\t' \
                          'simulation_repetition\n'
    output_table = open(args.output_prefix + '.gwas_runs.csv', 'w')
    output_table.write(output_table_header)
    output_bash = open(args.output_prefix + '.gwas_runs.sh', 'w')
    for af in allele_frequency:
        for es in effect_size:
            for ss in sample_size:
                # Number of cases and controls
                nca = int(int(ss) * float(case_control_ratio))
                nco = int(ss) - int(nca)
                check_positive_integer(nca, 'number of cases', None, None)
                check_positive_integer(nco, 'number of controls', None, None)
                # creating combination id
                combination_id = ''.join(choice(ascii_uppercase) for i in range(12))
                # files common to all phenotype replicates
                causal_variants = args.causal_loci
                output_prefix = output_dir + combination_id
                pyseer_phen_file = output_dir + combination_id + '.pyseer.phen'
                gwas_script_file = output_dir + combination_id + '.gwas_run.sh'
                # Creating one GWAS script per parameter combination (and for all phenotype replicates)
                gwas_script = open(gwas_script_file, 'w')
                command_line = '#!/bin/bash\n#GWAS run: ' + combination_id + '\n'
                command_line += 'python3 ' + code_directory + 'simulate_binary_phenotype_vcf.py' \
                                + ' --input_vcf ' + args.input_vcf \
                                + ' --out_prefix ' + output_prefix \
                                + ' --causal-loci ' + causal_variants \
                                + ' --simu-cc ' + str(nca) + ',' + str(nco) \
                                + ' --simu-rep ' + str(simulation_repetitions) \
                                + ' --allele-frequency ' + str(af) \
                                + ' --effect-size ' + str(es) \
                                + '\n'
                command_line += 'mv ' + output_prefix + '.phen ' + pyseer_phen_file + '\n'
                # For each phenotype replicate
                for ir in list(range(1, int(simulation_repetitions)+1, 1)):
                    replicate_id = combination_id + '_' + str(ir)
                    # saving parameter combination
                    newline_table = '\t'.join([str(replicate_id), str(af), str(es), str(ss), str(case_control_ratio),
                                               str(ir)])
                    output_table.write(newline_table + '\n')
                    # creating one GWAS script per parameter combination
                    # Phenotype file column to use in pyseer_phen_file
                    phenotype_column = 'phenotype' + str(ir)
                    pyseer_output = output_dir + replicate_id + '.pyseer_output.csv'
                    command_line += 'python3 ' + args.pyseer_path + ' --lmm ' \
                                    + ' --vcf ' + args.input_vcf \
                                    + ' --phenotypes ' + pyseer_phen_file \
                                    + ' --phenotype-column ' + str(phenotype_column) \
                                    + ' --similarity ' + args.similarity \
                                    + ' --min-af 0.001 --print-filtered '
                    if args.burden is not None:
                        command_line += ' --burden ' + args.burden
                    command_line += ' > ' + pyseer_output + '\n'
                # command_line += 'rm ' + causal_variants + '\n'
                # command_line += 'rm ' + pyseer_phen_file + '\n'
                command_line += 'rm ' + output_dir + combination_id + '.log' + '\n'
                command_line += 'rm ' + output_dir + combination_id + '.par' + '\n'
                command_line += 'fi\n'
                gwas_script.write(command_line)
                gwas_script.close()
                job_submission_command = ' '.join(args.job_submission_command)
                # replacing string job_id with combination_id
                job_submission_command = job_submission_command.replace('job_id', str(combination_id))
                output_bash.write(job_submission_command + ' "bash ' + gwas_script_file + '"\n')
                print(command_line)
    output_table.close()
    output_bash.close()


if __name__ == "__main__":
    _main()