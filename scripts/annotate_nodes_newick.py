#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import Phylo
import dendropy
import subprocess


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to extract ALL taxa ids from a phylogenetic tree in Newick format"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i", "--input_tree", action="store", dest="input_tree",
        help="Input phylogenetic tree file in Newick format",
        required=True, metavar="INPUT_TREE")
    parser.add_argument(
        "-n", "--node_prefix", action="store", dest="node_prefix",
        help="Node prefix used to label internal nodes (e.g. \"node\")",
        required=False, default="node", metavar="NODE_PREFIX")
    parser.add_argument(
        "-o", "--output_tree", action="store", dest="output_tree",
        help="Output phylogenetic tree file in Newick format",
        required=True, metavar="OUTPUT_TREE")

    return parser.parse_args()


def is_tree_valid(input_tree):
    try:
        Phylo.read(input_tree, 'newick')
        tree = dendropy.Tree.get_from_path(input_tree, 'newick')
    except:
        print("Error with the input starting tree: Is it a valid Newick file?")
        return 0
    return 1


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
    if not os.path.exists(args.input_tree):
        logging.error(f'Input phylogenetic tree file {args.input_tree} not found!')
        sys.exit(-1)

    # Loading phylogenetic tree
    logging.info(f"Opening input tree file {args.input_tree}")
    if is_tree_valid(args.input_tree) == 0:
        logging.error(f'Input phylogenetic tree {args.input_tree} is invalid!')
        sys.exit(-1)
    tree = dendropy.Tree.get_from_path(args.input_tree, 'newick', preserve_underscores=True)

    # Printing basic information on tree before removing taxa
    logging.info(f"Input tree information (before removing taxa)")
    print('Number of taxa in tree: ' + str(len(tree.taxon_namespace)))
    print('Number of leaf nodes in tree: ' + str(len(tree.leaf_nodes())))
    print('Number of internal nodes in tree: ' + str(len(tree.internal_nodes())))

    # Adding a name to all internal nodes
    logging.info(f"Adding a name to all internal nodes ...")
    for idx, node in enumerate(tree.internal_nodes()):
        node_label = args.node_prefix + str(idx)
        node.label = node_label
        # print(node.label)

    # Writing output file
    logging.info(f"Saving output tree file {args.output_tree}")
    tree.write(path=args.output_tree, schema="newick")

    # Removing single quotes introduced by tree.write function
    command_line = ["cat", args.output_tree, "|",
                    "sed \"s/'//g\" > ", args.output_tree + '_tmp', "&&",
                    'mv', args.output_tree + '_tmp', args.output_tree]
    command_line_string = ' '.join(command_line)
    run_command_shell_string(command_line_string)
    

if __name__ == "__main__":
    _main()
