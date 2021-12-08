#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import random

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to randomly subsample a subset of individuals in a PLINK-formatted .fam file\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('I/O arguments')
    group.add_argument(
        "-i", "--plink_fam_in", action="store", dest="plink_fam_in",
        help="Input PLINK .fam file",
        required=True, metavar="LIST"
    )
    group.add_argument(
        "-s", "--subset_size", action="store", dest="subset_size",
        help="Number of individuals to subsample from --plink_fam_in",
        required=True, metavar="OUTPUT_PREFIX"
    )
    group.add_argument(
        "-o", "--plink_fam_out", action="store", dest="plink_fam_out",
        help="Output PLINK .fam file with only subsample",
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
    if not os.path.exists(args.plink_fam_in):
        logging.error(f'Input file in --plink_fam_in {args.plink_fam_in} not found!')
        sys.exit(-1)

    # Making sure subset_size is a positive intenger
    check_positive_integer(args.subset_size, '--subset_size', None, None)

    # Saving items from input list
    all_individuals = list()
    all_fam_lines = dict()
    with open(args.plink_fam_in, 'r') as input_file:
        for line in input_file:
            ind = line.strip().split(' ')[0]
            all_individuals.append(ind)
            all_fam_lines[ind] = line.strip()
    logging.info(f'Number of individuals from --plink_fam_in {args.plink_fam_in} is {str(len(all_individuals))}')

    # Make sure size of subset (--subset_size) is not greater than number of items in input list
    if int(args.subset_size) > len(all_individuals):
        logging.error(f'subset_size ({args.subset_size}) is greater than number of items in input list '
                      f'({str(len(all_individuals))})!')
        sys.exit(-1)

    # Raise warning if non-unique/repeated items found in list
    all_unique_individuals = list(set(all_individuals))
    if len(all_unique_individuals) < len(all_individuals):
        repeated_individuals = list(set([x for x in all_individuals if all_individuals.count(x) > 1]))
        logging.warning(f'There are individuals repeated in input .fam file {args.plink_fam_in}')
        logging.warning(f"{'Repeated individuals: ' + ' '.join(repeated_individuals)}")

    # Random subsample from input list
    logging.info(f'Randomnly sampling {str(args.subset_size)} from {args.plink_fam_in}')
    sampled_individuals = random.sample(all_individuals, int(args.subset_size))

    # Saving subsampled list
    logging.info(f'Writing sampled items into {args.plink_fam_out}')
    output = open(args.plink_fam_out, 'w')
    for ind in sampled_individuals:
        output.write(all_fam_lines[ind] + '\n')
    output.close()


if __name__ == "__main__":
    _main()