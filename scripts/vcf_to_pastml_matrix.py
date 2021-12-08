#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from cyvcf2 import VCF
# import vcf

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to convert a multi-sample VCF file to matrix csv file compatible with PastML.\n"\
                  "Where each row represents a sample and each column a variant\n" \
                  "Variant ids (in header) are created as chr.pos.ref.alt\n" \
                  "Each cell contains the reference or alternative nucleotide alleles\n" \
                  "Multi-allelic sites not supported.\n" \
                  "Input VCF needs to have multi-allelic sites split beforehand. Use bcftools norm -m -\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-v", "--input_vcf", action="store", dest="input_vcf",
        help="multi-sample VCF file",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output table csv file compatible with PastML",
        required=True, metavar="OUTPUT_TABLE"
    )

    return parser.parse_args()


def parse_allele(allele):
    if '/' in allele:
        allele = allele.split('/')[0]
    return allele


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
    input_files = [args.input_vcf]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Opening VCF file and saving sample alleles
    logging.info(f"Opening VCF file and saving sample alleles... (this may take several minutes)\n")
    alleles = dict()    # dictionary to save sample alleles: alleles{sample_id} = sample_alleles
    delim = '\t'
    num_vcf_records = 0
    saved_var_ids = dict()
    vcf = VCF(args.input_vcf)
    sample_ids = vcf.samples
    for variant in vcf:
        var_id = str(variant.CHROM) + '.' + str(variant.POS) + '.' + str(variant.REF) + '.' + str(variant.ALT[0])
        print(var_id)
        num_vcf_records = num_vcf_records + 1
        for idx, allele in enumerate(variant.gt_bases):
            saved_var_ids[var_id] = 0
            sample_id = str(sample_ids[idx])
            allele = parse_allele(allele)
            # Saving sample allele in alleles dictionary
            if alleles.get(sample_id) is None:
                alleles[sample_id] = str(allele)
            else:
                alleles[sample_id] = alleles[sample_id] + delim + str(allele)
    logging.info("Number of VCF records: " + str(num_vcf_records))
    logging.info("Number of VCF records saved: " + str(len(saved_var_ids)))
    logging.info(f"Opening VCF file and saving sample alleles - DONE\n")

    # Saving sample alleles into table format
    logging.info(f"Saving sample alleles into table format... (this may take several minutes)\n")
    output_table = open(args.output_table, 'w')

    # Creating file header
    logging.info(f"Creating file header... (this may take several minutes)\n")
    table_header = 'sample'
    for var_id in saved_var_ids:
        table_header = table_header + delim + var_id
    logging.info("Number of variants in header " + str(len(table_header.split('\t'))-1))
    output_table.write(table_header + '\n')

    for sample in alleles:
        newline = sample + delim + alleles[sample]
        logging.info("number of variants for sample " + sample + " " + str(len(newline.split('\t'))-1))
        output_table.write(newline + '\n')
    output_table.close()
    logging.info(f"Saving sample alleles into table format - DONE\n")


if __name__ == "__main__":
    _main()