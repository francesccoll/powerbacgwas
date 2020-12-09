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
    description = "Script to perform ancestral state reconstruction for Roary pan-genome using PastML"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-v", "--input_pastml_table", action="store", dest="input_pastml_table",
        help="matrix of genes obtained from a gene_presence_absence.csv file using the script roary_to_pastml_matrix.py",
        required=True, metavar="INPUT_TABLE"
    )
    group.add_argument(
        "-t", "--input_tree", action="store", dest="input_tree",
        help="input phylogenetic tree",
        required=True, metavar="INPUT_TREE"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output table with gene alleles (0/1) and their ancestral nodes",
        required=True, metavar="OUTPUT_TABLE"
    )
    group.add_argument(
        "-s", "--output_steps", action="store", dest="output_steps",
        help="output table with number of homoplasies per gene",
        required=True, metavar="OUTPUT_STEPS"
    )
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-p", "--processes", action="store", dest="processes",
        help="Number of processes (cores) to use in multiprocessing",
        required=False, metavar="PROCESSES", default=4
    )
    group.add_argument(
        "-a", "--pastml_path", action="store", dest="pastml_path",
        help="Full path to PastML executable", required=False, default='pastml', metavar="PASTML_PATH"
    )

    return parser.parse_args()


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


def duplicated(input_list):
    """ Returns a list of indices where duplicated items in a list are, except for the first instance
    :param list: list of items
    :return: list of indices of duplicated items
    """
    list_dict = dict()
    duplicated_idx = list()
    for idx, item in enumerate(input_list):
        if item in list_dict:
            duplicated_idx.append(idx)
        list_dict[item] = 0
    return duplicated_idx


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
    dependencies = [args.pastml_path]
    for dependency in dependencies:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)

    # Making sure input files exist
    input_files = [args.input_pastml_table, args.input_tree]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Slicing csv table: to run each slice separately using PastML
    logging.info(f"Reading {args.input_pastml_table}. This may take several minutes...")
    df = pd.read_csv(args.input_pastml_table, sep='\t')
    logging.info(f"Reading {args.input_pastml_table}. DONE")
    print("Number of samples: " + str(df.shape[0]) + " number of variants " + str(df.shape[1]))

    logging.info(f"Slicing {args.input_pastml_table}. This may take several minutes...")
    by = 1000
    number_columns = df.shape[1]
    column_names = df.columns.values.tolist()
    input_table_slices = []
    for c in range(1, number_columns, by):
        start = c
        end = c + by - 1
        end = df.shape[1] - 1 if end > number_columns else end
        slice_name = '_slice_' + str(start) + '_' + str(end)
        slice_file_name = args.input_pastml_table + slice_name
        input_table_slices.append(slice_file_name)
        if not os.path.isfile(slice_file_name):  # if slice file does not exist
            column_range = [0] + list(range(start, end + 1))
            df_slice = df.iloc[:, column_range]
            logging.info(f"Writing {slice_file_name}")
            df_slice.to_csv(slice_file_name, header=True, sep='\t', index=False)
    logging.info(f"Slicing {args.input_pastml_table}. DONE")

    # Running PastML (Ancestral reconstruction) for each slice - multi-processing
    pastml_dir = 'tmp_pastml_dir'
    logging.info(f"Running PastML (Ancestral reconstruction) for each slice - multi-processing...")
    pool = Pool(processes=int(args.processes))
    command_lines = []
    for input_table_slice in input_table_slices:
        command_line = ["pastml", "--tree", args.input_tree,
                        "--data", input_table_slice,
                        "--out_data", input_table_slice + '_pastml',
                        "--work_dir", pastml_dir,
                        "--prediction_method", "ACCTRAN"]
        command_line_string = ' '.join(command_line)
        if not os.path.isfile(input_table_slice + '_pastml'):  # if PastML output does not exist
            command_lines.append(command_line_string)
    pool.map(run_command_shell_string, command_lines)
    logging.info(f"Running PastML (Ancestral reconstruction) for each slice - DONE")

    # Processing PastML output files to keep one ACR solution per node
    # NOTE: PastML may output multiple lines (presumably multiple ancestral character
    # reconstruction (ACR) solutions) per node. Only the first occurrence/line per node is kept.
    logging.info(f"Processing PastML output files...")
    for input_table_slice in input_table_slices:
        if not os.path.isfile(input_table_slice + '_pastml_processed'):
            df_pastml_slice = pd.read_csv(input_table_slice + '_pastml', sep='\t', dtype=object)
            # extracting VCF table columns names
            df_roary_slice = pd.read_csv(input_table_slice, sep='\t', dtype=object)
            roary_column_names_prv = df_roary_slice.columns.values.tolist()
            roary_column_names_new = list()
            for cn in roary_column_names_prv:
                roary_column_names_new.append(cn.replace('.', ''))
            roary_column_names_new[0] = 'node'
            # matching order of columns between VCF table and PASTML output
            df_pastml_slice = df_pastml_slice.reindex(columns=roary_column_names_new)
            # identifying repeated nodes (i.e. ACR solutions)
            rows_to_rm = duplicated(df_pastml_slice.iloc[:, 0].tolist())
            df_pastml_slice.drop(rows_to_rm, axis='index', inplace=True)
            df_pastml_slice.to_csv(input_table_slice + '_pastml_processed', header=True, sep='\t', index=False)
    logging.info(f"Processing PastML output files. DONE.")

    # Joining all slices back together
    logging.info(f"Joining all slices back together. This may take several minutes...")
    if not os.path.isfile(args.output_table + '_tmp'):
        paste_command = "paste -d'\t' "
        for input_table_slice in input_table_slices:
            paste_command += input_table_slice + '_pastml_processed '
        paste_command += '> ' + args.output_table + '_tmp'
        run_command_shell_string(paste_command)

    logging.info(f"Creating final PastML output file {args.output_table}. This may take several minutes...")
    if not os.path.isfile(args.output_table):
        df_pastml_tmp = pd.read_csv(args.output_table + '_tmp', sep='\t', dtype=object)
        # the first column id needs to be changed
        df_pastml_tmp.rename(columns={"node": "sample"}, errors="raise", inplace=True)
        # node (former sample) column ids need to be removed
        cols_to_rm = [col for col in df_pastml_tmp.columns if 'node' in col]
        df_pastml_tmp.drop(cols_to_rm, axis=1, inplace=True)
        # column ids need to be changed back to original VCF table column ids
        column_names_tmp = df_pastml_tmp.columns.values.tolist()
        df_pastml_tmp.rename(columns=dict(zip(column_names_tmp, column_names)), inplace=True)
        # writing final PastML output table
        df_pastml_tmp.to_csv(args.output_table, header=True, sep='\t', index=False)
        logging.info(f"Creating final PastML output file {args.output_table}. DONE.")
    else:
        print('\tOutput file ' + args.output_table + 'already exists')

    # Extracting number of homoplasies per character (i.e. steps)
    logging.info(f"Extracting number of homoplasies per character (i.e. steps)...")
    if not os.path.isfile(args.output_steps):
        pastml_output_steps = open(args.output_steps, 'w')
        pastml_output_steps.write('variant_id\tsteps\n')
        variant_ids = column_names
        variant_ids.pop(0)
        for variant_id in variant_ids:
            variant_id_alt = variant_id.replace('.', '').replace('-', '').replace('(', '').replace(')', '').replace('\'', '')
            params_file = pastml_dir + '/params.character_' + variant_id_alt + '.method_ACCTRAN.tab'
            if os.path.isfile(params_file):
                with open(params_file, 'r') as params_lines:
                    for line in params_lines:
                        fields = line.strip().split('\t')
                        if 'steps' in fields[0]:
                            newline = variant_id + '\t' + fields[1]
                            pastml_output_steps.write(newline + '\n')
            else:
                logging.error(f"{params_file} file not found")
        pastml_output_steps.close()
        logging.info(f"Extracting number of homoplasies per character (i.e. steps). DONE")
    else:
        print('\tOutput file ' + args.output_steps + 'already exists')

    # Delete temporary files
    # logging.info(f"Deleting temporary files...")
    # files_to_rm = list()
    # for input_table_slice in input_table_slices:
    #     files_to_rm.append(input_table_slice)
    #     files_to_rm.append(input_table_slice + '_pastml')
    #     files_to_rm.append(input_table_slice + '_pastml_processed')
    # files_to_rm.append(args.output_table + '_tmp')
    # files_to_rm.append('tmp_pastml_dir')
    # for file_to_rm in files_to_rm:
    #     run_command_shell_string(' '.join(["rm -r -f", file_to_rm]))


if __name__ == "__main__":
    _main()