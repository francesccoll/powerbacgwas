#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
import pandas as pd
from multiprocessing import Pool

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to calculate number of ancestral state changes (homoplasies) per region from PastML output"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-p", "--pastml_output_table", action="store", dest="pastml_output_table",
        help="PastML output table with variants and their ancestral states produced by script "
             "ancestral_state_reconstruction.py",
        required=True, metavar="PASTML_TABLE"
    )
    group.add_argument(
        "-t", "--input_tree", action="store", dest="input_tree",
        help="input phylogenetic tree",
        required=True, metavar="INPUT_TREE"
    )
    group.add_argument(
        "-r", "--regions_file", action="store", dest="regions_file",
        help="The regions file must contain one region per line, with their name and the bcftools style region "
             "co-ordinates delimited by tab (locus_name\tchr_name:start-end).",
        required=True, metavar="REGIONS_FILE"
    )
    group.add_argument(
        "-o", "--output_steps", action="store", dest="output_steps",
        help="output table with number of ancestral state changes (homoplasies) per region",
        required=True, metavar="OUTPUT_STEPS"
    )
    group.add_argument(
        "-f", "--calculate_changes_path", action="store", dest="calculate_changes_path",
        help="Full path to PastML calculate_changes.py script (found in pastml-master/pastml/utilities)",
        required=False, default="calculate_changes.py",  metavar="CC_PATH"
    )
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-c", "--processes", action="store", dest="processes",
        help="Number of processes (cores) to use in multiprocessing",
        required=False, metavar="PROCESSES", default=4
    )
    group.add_argument(
        "-d", "--tmp_dir", action="store", dest="tmp_dir", help="Directory to saved temporary files",
        required=False, metavar="TMP_DIR", default="./"
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


def check_path_exists(my_dir):
    if not os.path.exists(my_dir):
        logging.error(f'Directory {my_dir} not found!')
        sys.exit(-1)


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
    input_files = [args.pastml_output_table, args.input_tree, args.regions_file]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    check_dependency(args.calculate_changes_path)

    # Making sure input directories exist
    tmp_dir = os.path.abspath(args.tmp_dir)
    check_path_exists(tmp_dir)
    tmp_dir += '/'

    # Saving region ids
    region_ids = dict()  # dictionary to save regions' ids
    region_coordinates = dict()  # dictionary to save regions' chromosome id and coordinates
    if args.regions_file is not None:
        logging.info(f"Reading regions file {args.regions_file}...")
        file = open(args.regions_file, "r")
        for line in file:
            [region_id, rest] = line.strip().split('\t')
            [chr_id, coordinates] = rest.split(':')
            [start, end] = coordinates.split('-')
            check_positive_integer(start, 'start coordinate in line: '+line.strip(), None, None)
            check_positive_integer(end, 'end coordinate in line: '+line.strip(), None, None)
            region_ids[region_id] = 1
            for pos in range(int(start), int(end)):
                region_coordinates[chr_id + '.' + str(pos)] = region_id
        print("\tNumber of regions saved: ", str(len(region_ids)))
        logging.info(f"Reading regions file {args.regions_file}. DONE.")

    # Opening PastML output table to save variant ids, and assign them to their region_id
    logging.info(f"Reading {args.pastml_output_table}. This may take several minutes...")
    df = pd.read_csv(args.pastml_output_table, sep='\t')
    logging.info(f"Reading {args.pastml_output_table}. DONE")
    print("Number of samples: " + str(df.shape[0]) + " number of variants " + str(df.shape[1]))

    logging.info(f"Extracting variant ids from {args.pastml_output_table} and linking them to their region id...")
    region_variants = dict()
    variants_saved = 0
    column_names = df.columns.values.tolist()
    column_names.pop(0)
    print("Length of column_names: " + str(len(column_names)))
    for var_id in column_names:
        print(var_id)
        [chr, pos, *rest] = var_id.strip().split('.')  # expected var_id: chr.pos.ref.alt
        var_id_tmp = chr + '.' + pos
        if var_id_tmp in region_coordinates:
            variants_saved += 1
            region_id = region_coordinates[var_id_tmp]
            if region_id in region_variants:
                region_variants[region_id] = region_variants[region_id] + [var_id]
            else:
                region_variants[region_id] = [var_id]
    logging.info(f"Extracting variant ids from {args.pastml_output_table} and linking them to their region id. DONE")
    logging.info(f"{str(variants_saved)} variants saved in {str(len(region_variants))} regions.")

    # Slicing PastML output table into files containing variants from same region
    logging.info(f"Slicing PastML output table into files containing variants from same region...")
    for region_id in region_variants:
        variants_ids = ['sample'] + region_variants[region_id]
        df_tmp = df[variants_ids]
        df_tmp.to_csv(tmp_dir + region_id + '_pastml_table.csv', header=True, sep='\t', index=False)
    logging.info(f"Region-sliced PastML output files saved in {tmp_dir}")
    logging.info(f"Slicing PastML output table into files containing variants from same region. DONE")

    # Running PastML calculate_changes.py script for each region
    logging.info(f"Running PastML calculate_changes.py script for each region...")
    pool = Pool(processes=int(args.processes))
    command_lines = []
    for region_id in region_variants:
        region_pastml_table = tmp_dir + region_id + '_pastml_table.csv'
        region_cc_output = tmp_dir + region_id + '_pastml_cc.csv'
        # NOTE: script calculate_changes.py needs to be called with python3 as it lack top line interpreter
        command_line = ["python3", args.calculate_changes_path,
                        "--tree", args.input_tree,
                        "--acr", region_pastml_table,
                        "--out_log", region_cc_output]
        command_line_string = ' '.join(command_line)
        if not os.path.isfile(region_cc_output):  # if PastML output does not exist
            command_lines.append(command_line_string)
    pool.map(run_command_shell_string, command_lines)
    logging.info(f"Running PastML calculate_changes.py script for each region. DONE")

    # Savings PastML calculate_changes.py script results
    logging.info(f"Savings PastML calculate_changes.py script results...")
    region_steps = dict()
    for region_id in region_variants:
        region_cc_output = tmp_dir + region_id + '_pastml_cc.csv'
        file = open(region_cc_output, "r")
        file.readline()
        steps = 0
        for line in file:
            [fro, to, count, normalized_count] = line.strip().split('\t')
            steps = steps + int(float(count))
        region_steps[region_id] = str(steps)
        print("region_steps[", region_id, "] --> ", region_steps[region_id])
    logging.info(f"Savings PastML calculate_changes.py script results. DONE.")

    # Saving output table with number of ancestral state changes (homoplasies) per region
    logging.info(f"Writing output table with number of ancestral state changes per region in {args.output_steps}")
    pastml_output_steps = open(args.output_steps, 'w')
    pastml_output_steps.write('region_id\tsteps\n')
    for region_id in region_ids:
        steps = 0
        if region_id in region_steps:
            steps = region_steps[region_id]
        newline = region_id + '\t' + str(steps)
        pastml_output_steps.write(newline + '\n')
    pastml_output_steps.close()
    logging.info(f"Writing output table with number of ancestral state changes per region in {args.output_steps}. DONE.")

    # Delete temporary files
    logging.info(f"Deleting temporary files...")
    for region_id in region_variants:
        region_pastml_table = tmp_dir + region_id + '_pastml_table.csv'
        region_cc_output = tmp_dir + region_id + '_pastml_cc.csv'
        run_command_shell_string(' '.join(["rm", region_pastml_table, region_cc_output]))


if __name__ == "__main__":
    _main()