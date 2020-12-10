#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess


# ---------------------------------------------------------------------------------------------------------------------
# Developing notes
# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to process power calculation GWAS runs prepared by script prepare_gwas_runs.py"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required: I/O files')
    group.add_argument(
        "-i", "--gwas_runs_in_table", action="store", dest="gwas_runs_in_table",
        help="Table produced by prepare_gwas_runs.py, containing GWAS parameters and unique identifiers for "
             "each parameter combination (GWAS run)",
        required=True, metavar="INPUT_TABLE"
    )
    group.add_argument(
        "-t", "--variant_type", action="store", dest="variant_type",
        help="VCF (v), VCF burden (b) or Roary presence/absence (r)",
        required=False, default='v', metavar="VARIANT_TYPE"
    )
    group.add_argument(
        "-o", "--gwas_runs_out_table", action="store", dest="gwas_runs_out_table",
        help="Output table with same information as in gwas_runs_in_table, but extended with causal variants, causal "
             "variant frequency in cases and controls, and PySeer p-values. This table is used by the "
             "script plot_gwas_runs.R",
        required=True, metavar="OUTPUT_TABLE"
    )
    group.add_argument(
        "-l", "--causal-loci", action="store", dest="causal_loci",
        help="GCTA-formatted file with causal loci used in sub-sampling approach/scripts",
        required=False, metavar="CAUSAL_LOCI"
    )
    group = parser.add_argument_group('Optional arguments:')
    group.add_argument(
        "-d", "--output_dir", action="store", dest="output_dir",
        help="Output directory where GWAS output files are stored (same as --output_dir in prepare_gwas_runs.py)",
        required=False, default='./', metavar="OUTPUT_DIR"
    )

    return parser.parse_args()


def check_file_exists(my_file):
    if not os.path.isfile(my_file):
        logging.error(f'File {my_file} not found!')
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

    # Making sure input files exist
    logging.info('Making sure input files exist...')
    input_files = [args.gwas_runs_in_table]
    if args.causal_loci is not None:
        input_files.append(args.causal_loci)
    for input_file in input_files:
        check_file_exists(input_file)
    logging.info('All required input files found.')

    # Make sure variant_type is either 'v', 'b' or 'r'
    valid_options = ['v', 'b', 'r']
    if args.variant_type not in valid_options:
            logging.error(f'value provided in option --variant_type {args.variant_type} must be \'v\', \'b\' or \'r\'')
            sys.exit(-1)

    # Opening GWAS runs table
    output_table = open(args.gwas_runs_out_table, 'w')
    logging.info(f'Opening {args.gwas_runs_in_table}')
    with open(args.gwas_runs_in_table, 'r') as input_file:
        output_header_saved = False
        input_header = input_file.readline().strip()
        pyseer_header = 'variant\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tbeta-std-err\tvariant_h2\tnotes'
        output_header = input_header + '\tcausal_variant\tn_cases_sim\tn_controls_sim\tn_missing_sim\t' \
                       + pyseer_header + '\n'
        output_table.write(output_header)
        for input_line in input_file:
            print(input_line)
            replicate_id = input_line.strip().split('\t')[0]  # made up of combination_id + '_' + phen_rep_num
            combination_id = replicate_id.strip().split('_')[0]
            phen_rep_num = int(replicate_id.strip().split('_')[1])  # phenotype replicate number
            # defining prepare_gwas_runs.py output files to load
            causal_variants_file = args.output_dir + combination_id + '.causal_variants.txt'
            if args.variant_type == 'b':
                causal_variants_file = args.output_dir + combination_id + '.causal_regions.txt'
            # for sub-sampling approach/scripts, the known causal loci need to be used instead
            if args.causal_loci is not None:
                causal_variants_file = args.causal_loci
                # to do: add regions_file? or specify causal regions in args.causal_loci
            pyseer_phen_file = args.output_dir + combination_id + '.pyseer.phen'
            pyseer_output = args.output_dir + replicate_id + '.pyseer_output.csv'
            # There is a possibility that causal variants meeting a parameter combination could not be found by
            # sample_casual_variants_from_vcf.py script, in which case a GWAS could not be run (Pyseer output absent)
            if os.path.isfile(causal_variants_file):
                # Extracting causal variants
                causal_variants = dict()
                with open(causal_variants_file, 'r') as cv_file:
                    for cv_line in cv_file:
                        cv = cv_line.strip().split('\t')[0]
                        causal_variants[cv] = 1
                # Extracting number of simulated cases and controls (from simulated phenotype file)
                # to do: this will need to change for quantitative phenotypes
                (n_cases, n_controls, n_missing) = ('NA', 'NA', 'NA')
                if os.path.isfile(pyseer_phen_file):
                    (n_cases, n_controls, n_missing) = (0, 0, 0)
                    with open(pyseer_phen_file, 'r') as pp_file:
                        for pp_line in pp_file:
                            phenotype = pp_line.strip().split('\t')[phen_rep_num]
                            if phenotype == '0':
                                n_controls += 1
                            elif phenotype == '1':
                                n_cases += 1
                            else:
                                n_missing += 1
                # print('n_cases ' + str(n_cases) + ' n_controls ' + str(n_controls) + ' n_missing ' + str(n_missing))
                # Extracting PySeer output line for causal variants
                # check_file_exists(pyseer_output)
                # to do: this will need to change for multiple loci
                if os.path.isfile(pyseer_output):
                    with open(pyseer_output, 'r') as po_file:
                        po_file.readline()
                        for po_line in po_file:
                            # NOTE: VCF variant Ids created by Pyseer (CHROM_POS_REF_ALT) differ from those created by
                            # sample_casual_variants_from_vcf.py (CHROM.POS.REF.ALT) in the field delimiter
                            variant_id = ''
                            if args.variant_type == 'v':
                                variant_id = po_line.strip().split('\t')[0].replace('_', '.')
                            else:
                                variant_id = po_line.strip().split('\t')[0]
                            if variant_id in causal_variants:
                                print(po_line)
                                out_line = input_line.strip() + '\t' + variant_id + '\t' + str(n_cases) + '\t'\
                                           + str(n_controls) + '\t' + str(n_missing) + '\t' + po_line.strip() + '\n'
                                output_table.write(out_line)
                else:
                    logging.warning(f'{pyseer_output} file not found.')
            else:
                # if causal variants file not found, add 12 missing fields to output table (4 + 8 PySeer output fields)
                out_line = input_line.strip() + '\t' + '\t'.join(['NA']*12) + '\n'
                output_table.write(out_line)
                print(causal_variants_file + " file not found")
    output_table.close()


if __name__ == "__main__":
    _main()