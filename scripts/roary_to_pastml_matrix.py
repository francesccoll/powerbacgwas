#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to convert a Roary gene_presence_absence.csv file to matrix csv file compatible with PastML.\n"\
                  "Where each row represents a sample and each column a gene\n" \
                  "Gene ids (in header) are taken from the Gene gene_presence_absence.csv file\n" \
                  "Each cell contains either 0 (absence) or 1 (presence)\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-g", "--gene_presence_absence", action="store", dest="gene_presence_absence",
        help="Roary gene_presence_absence.csv file",
        required=True, metavar="INPUT_TABLE"
    )
    group.add_argument(
        "-f", "--input_format", action="store", dest="input_format",
        help="Roary file format: (1) gene_presence_absence.csv or (2) gene_presence_absence.Rtab",
        required=True, metavar="INPUT_FORMAT"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output table csv file compatible with PastML",
        required=True, metavar="OUTPUT_TABLE"
    )
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-s", "--samples_list", action="store", dest="samples_list",
        help="Samples to keep from Roary gene_presence_absence.csv file.\nOne sample id per line expected.\n",
        required=False, metavar="SAMPLES_LIST"
    )

    return parser.parse_args()


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
    input_files = [args.gene_presence_absence]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Making sure input_format is either 1 or 2
    if args.input_format != '1':
        if args.input_format != '2':
            logging.error(f'Input format chosen (-f, --input_format) \'{args.input_format}\' must be either 1 or 2.')
            sys.exit(-1)

    # Reading input sample list if provided
    samples_to_keep = dict()
    if args.samples_list is not None:
        if not os.path.exists(args.samples_list):
            logging.error(f'Input file {args.samples_list} not found!')
            sys.exit(-1)
        else:
            with open(args.samples_list, 'r') as lines:
                for line in lines:
                    samples_to_keep[line.strip()] = 0

    # Opening gene_presence_absence file and saving sample alleles
    if args.input_format == '1':
        logging.info(f"Opening gene_presence_absence.csv file and saving sample alleles...")
        alleles = dict()    # dictionary to save sample alleles: alleles[gene_id][sample_id] = sample_allele (0/1)
        samples = list()
        with open(args.gene_presence_absence, 'r') as lines:
            fields = lines.readline().strip().split('","')  # sample ids are extracted from first line
            samples = fields[14:len(fields)]
            samples = [s.replace('"', '') for s in samples]
            for line in lines:
                fields = line.strip().split('","')
                gene_id = fields[0].replace('"', '')
                presabs = fields[14:len(fields)]
                for idx, allele in enumerate(presabs):
                    allele.replace('"', '')
                    allele = '0' if allele == '' else '1'
                    sample_id = samples[idx]
                    if alleles.get(gene_id) is None:  # Saving sample allele in alleles dictionary
                        alleles[gene_id] = {}
                        alleles[gene_id][sample_id] = allele
                    else:
                        alleles[gene_id][sample_id] = allele
        logging.info(f"Opening gene_presence_absence.csv file and saving sample alleles. DONE!")

    if args.input_format == '2':
        logging.info(f"Opening gene_presence_absence.Rtab file and saving sample alleles...")
        alleles = dict()  # dictionary to save sample alleles: alleles[gene_id][sample_id] = sample_allele (0/1)
        samples = list()
        with open(args.gene_presence_absence, 'r') as lines:
            fields = lines.readline().strip().split('\t')  # sample ids are extracted from first line
            samples = fields
            samples.pop(0)
            for line in lines:
                fields = line.strip().split('\t')
                gene_id = fields[0]
                presabs = fields
                presabs.pop(0)
                for idx, allele in enumerate(presabs):
                    sample_id = samples[idx]
                    if alleles.get(gene_id) is None:  # Saving sample allele in alleles dictionary
                        alleles[gene_id] = {}
                        alleles[gene_id][sample_id] = allele
                    else:
                        alleles[gene_id][sample_id] = allele
        logging.info(f"Opening gene_presence_absence.Rtab file and saving sample alleles. DONE!")

    # Samples to include in output file
    if len(samples_to_keep) == 0:  # if no input sample list provided, keep all samples
        samples_to_keep = samples
    else:  # make sure all samples in input sample list are in gene_presence_absence.csv file
        for sample_id in samples_to_keep:
            if sample_id not in samples:
                logging.error(f"{sample_id} in {args.samples_list} not found in {args.gene_presence_absence}")

    logging.info(f"Saving alleles in PastML format...")
    output_table = open(args.output_table, 'w')
    header = 'sample'
    for gene_id in alleles:
        header += '\t' + gene_id
    output_table.write(header + '\n')
    for sample_id in samples_to_keep:
        newline = sample_id
        for gene_id in alleles:
            newline += '\t' + alleles[gene_id][sample_id]
        output_table.write(newline + '\n')
    logging.info(f"Saving alleles in PastML format. DONE!")


if __name__ == "__main__":
    _main()