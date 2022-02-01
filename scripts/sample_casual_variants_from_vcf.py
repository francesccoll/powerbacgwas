#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import random
import numpy
from cyvcf2 import VCF


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to sample causal variants from a multi-sample VCF file to then run GCTA.\n" \
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
        "-p", "--pastml_steps_file", action="store", dest="pastml_steps_file",
        help="table with number of homoplasies per variant created by script ancestral_state_reconstruction.py. If a "
             "--regions_file is provided, ensure the --pastml_steps_file provided corresponds to that created by "
             "calculate_changes_per_region.py",
        required=True, metavar="PASTML_STEPS"
    )
    group.add_argument(
        "-n", "--number_variants", action="store", dest="number_variants",
        help="number of randomly sampled variants (integer)",
        required=True, metavar="NUMBER_VARIANTS"
    )
    group.add_argument(
        "-f", "--allele_frequency_range", action="store", dest="allele_frequency_range",
        help="range of allele frequency of selected variants, separated by comma (e.g. 0.05,0.10 for 5%% to 10%%)",
        required=True, metavar="ALLELE_FREQUENCY_RANGE"
    )
    group.add_argument(
        "-s", "--degree_homoplasy_range", action="store", dest="degree_homoplasy_range",
        help="range of degree of homoplasy of selected variants (number of steps), separated by comma "
             "(e.g. 2,5 for 2 to 5 steps)",
        required=True, metavar="DEGREE_HOMOPLASY_RANGE"
    )
    group.add_argument(
        "-e", "--effect_size", action="store", dest="effect_size",
        help="effect size of selected variants, in odd ratio units (e.g. 2 for OR of 2)",
        required=True, metavar="EFFECT_SIZE"
    )
    group.add_argument(
        "-o", "--output_gcta_file", action="store", dest="output_gcta_file",
        help="file to save randomly sampled VCF variants in GCTA format",
        required=True, metavar="OUTPUT_TABLE"
    )
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-r", "--regions_file", action="store", dest="regions_file",
        help="A regions file can be provided if a burden test is to be run with script prepare_gwas_runs.py (--burden)."
             "The regions file can be used to specify gene coordinates and must contain one region per line, with "
             "their name and the bcftools style region co-ordinates delimited by tab (locus_name\tchr_name:start-end). "
             "If a regions file is provided, then allele frequencies and number of homoplasies will be calculated per "
             "region. The variants printed in the output GCTA-format file will correspond to all the variants observed "
             "in the regions meeting the chosen criteria.",
        required=False, metavar="REGIONS_FILE"
    )
    group.add_argument(
        "-u", "--output_regions_file", action="store", dest="output_regions_file",
        help="file to save randomly sampled regions",
        required=False, metavar="OUTPUT_REGIONS"
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
    input_files = [args.input_vcf, args.pastml_steps_file]
    if args.regions_file is not None:
        input_files.append(args.regions_file)
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Making sure both --output_regions_file and --regions_file are provided, if applicable
    if args.regions_file is not None:
        if args.output_regions_file is None:
            logging.error(f'Both a --regions_file and --output_regions_file file names must be provided.')
            sys.exit(-1)

    # Parsing parameters
    (af_from, af_to) = args.allele_frequency_range.split(',')
    (steps_from, steps_to) = args.degree_homoplasy_range.split(',')
    odds_ratio = args.effect_size
    number_variants = args.number_variants.strip()

    # Making sure input parameters meet expected values
    int_param = [steps_from, steps_to, number_variants]
    int_param_names = ['start of range in number of steps (-s)', 'end of range in number of steps (-s)',
                       'number of variants (-n)']
    for idx, param in enumerate(int_param):
        check_positive_integer(param, int_param_names[idx], None, None)
    float_param = [af_from, af_to]
    float_param_names = ['start of range in allele frequency (-f)', 'end of range in allele frequency (-f)']
    for idx, param in enumerate(float_param):
        check_positive_float(param, float_param_names[idx], 0, 1)
    check_positive_float(odds_ratio, "effect size (--effect_size)", None, None)

    # Printing chosen parameters
    print("Chosen variants parameters: ")
    print("\tNumber of causal variants: " + number_variants)
    print("\tAllele frequency: " + af_from + "-" + af_to)
    print("\tDegree of homoplasy: " + steps_from + "-" + steps_to)
    print("\tEffect size: " + odds_ratio)

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

    # Opening VCF file and calculating allele frequency
    allele_frequency = dict()    # dictionary to save allele frequency: allele_frequency[var_id] = sample_allele
    num_vcf_records = 0
    saved_var_ids = dict()
    region_variants = dict()
    vcf = VCF(args.input_vcf)
    missing_calls = ['N', '-', '.', 'NA']   # symbols used in VCFs to encode a missing genotype call
    if args.regions_file is None:
        logging.info(f"Opening VCF file and calculating allele frequency... (may take several minutes)")
        for variant in vcf:
            var_id = str(variant.CHROM) + '.' + str(variant.POS) + '.' + str(variant.REF) + '.' + str(variant.ALT[0])
            num_vcf_records = num_vcf_records + 1
            frequency = 0
            missing_calls_n = 0
            for allele in variant.gt_bases:
                saved_var_ids[var_id] = 0
                if str(allele) in missing_calls:
                    missing_calls_n = missing_calls_n + 1
                if str(allele) is str(variant.ALT[0]):
                    frequency = frequency + 1
            # allele_frequency[var_id] = frequency / len(vcf.samples)
            # NOTE: the allele frequency is calculated not including isolates with missing calls at the site
            allele_frequency[var_id] = frequency / (len(vcf.samples) - missing_calls_n)
        print("\tNumber of VCF records: ", str(num_vcf_records))
        print("\tNumber of VCF records saved: ", str(len(saved_var_ids)))
    else:
        logging.info(f"Opening VCF file and calculating allele frequency at region level... (may take several minutes)")
        region_samples = dict()  # dictionary to save samples with ALT alleles in each region
        sample_ids = vcf.samples
        for variant in vcf:
            var_id = str(variant.CHROM) + '.' + str(variant.POS)
            var_id_2 = str(variant.CHROM) + '.' + str(variant.POS) + '.' + str(variant.REF) + '.' + str(variant.ALT[0])
            num_vcf_records = num_vcf_records + 1
            var_samples = list()
            for idx, allele in enumerate(variant.gt_bases):
                saved_var_ids[var_id] = 0
                sample_id = str(sample_ids[idx])
                if str(allele) is not str(variant.REF):
                    var_samples.append(sample_id)
            if var_id in region_coordinates:
                region_id = region_coordinates[var_id]
                if region_id in region_variants:
                    region_variants[region_id] = region_variants[region_id] + [var_id_2]
                else:
                    region_variants[region_id] = [var_id_2]
                if region_id in region_samples:
                    region_samples[region_id] = region_samples[region_id] + var_samples
                else:
                    region_samples[region_id] = var_samples
        print("\tNumber of VCF records: ", str(num_vcf_records))
        print("\tNumber of VCF records saved: ", str(len(saved_var_ids)))
        print("\tCalculating allele frequency at region level...")
        for region_id in region_ids:
            frequency = 0
            if region_id in region_samples:
                var_samples_unique = list(set(region_samples[region_id]))
                frequency = len(var_samples_unique) / len(sample_ids)
            allele_frequency[region_id] = frequency
        print("\tCalculating allele frequency at region level. DONE")
    logging.info(f"Opening VCF file and calculating allele frequency. DONE")

    # Opening PastML steps output
    logging.info(f"Saving number of homoplasies (steps) per variant...")
    allele_steps = dict()  # dictionary to save allele frequency: allele_frequency[var_id] = sample_allele
    file = open(args.pastml_steps_file, "r")
    file.readline()
    for line in file:
        [var_id, steps] = line.strip().split('\t')
        allele_steps[var_id] = steps
        # print("allele_steps["+var_id+"] "+allele_steps[var_id])
        # checking var_id in pastml_steps_file match variant ids or region ids in VCF and region files, respectively
        if var_id not in allele_frequency:
            logging.error(f"{var_id} in {args.pastml_steps_file} not found in other input files. Check.")
    logging.info(f"Saving number of homoplasies (steps) per variant. DONE.")

    # Selecting variants meeting criteria (allele frequency and homoplasy steps)
    logging.info(f"Selecting variants meeting criteria...")
    set1 = set()
    for var_id, freq in allele_frequency.items():
        if (float(freq) >= float(af_from)) and (float(freq) <= float(af_to)):
            set1.add(var_id)
            # print('Frequency met: ' + var_id + ' ' + str(freq))
    set2 = set()
    for var_id, steps in allele_steps.items():
        if (int(steps) >= int(steps_from)) and (int(steps) <= int(steps_to)):
            set2.add(var_id)
            # print('Homoplasy met: ' + var_id + ' ' + str(steps))
    intersection_set = set1.intersection(set2)
    print(str(len(intersection_set)) + "/" + str(len(allele_frequency)) + " variants meeting criteria")
    logging.info(f"Selecting variants meeting criteria. DONE.")

    if len(intersection_set) < int(number_variants):
        logging.error(f"Number of variants meeting criteria {str(len(intersection_set))} lower than selected ones "
                      f"{str(number_variants)}.")
        sys.exit(-1)

    # Randomly sampling of variants
    logging.info(f"Randomly sampling of variants meeting criteria.")
    selected_variants = random.sample(intersection_set, int(number_variants))

    # Saving sampled causal variants in GCTA format
    # Input file format: causal.snplist (columns are SNP ID and effect size)
    output_gcta = open(args.output_gcta_file, 'w')
    if args.regions_file is None:
        for var_id in selected_variants:
            newline = var_id + '\t' + str(odds_ratio) + '\n'
            output_gcta.write(newline)
            print("\tselected_variant " + var_id)
            print("\tallele frequency " + str(allele_frequency[var_id]))
            print("\thomoplasy steps " + str(allele_steps[var_id]))
            print("\teffect size " + str(odds_ratio))
    else:
        output_regions = open(args.output_regions_file, 'w')
        for region_id in selected_variants:
            output_regions.write(region_id + '\n')
            print("\tselected_region " + region_id)
            print("\tallele frequency " + str(allele_frequency[region_id]))
            print("\thomoplasy steps " + str(allele_steps[region_id]))
            print("\teffect size " + str(odds_ratio))
            var_ids = region_variants[region_id]
            for var_id in var_ids:
                newline = var_id + '\t' + str(odds_ratio) + '\n'
                output_gcta.write(newline)
                print("\tselected_variant " + var_id)
        output_regions.close()
        logging.info(f"Randomly sampled causal regions written to {args.output_regions_file}")
    output_gcta.close()
    logging.info(f"Randomly sampled causal variants written to {args.output_gcta_file}")
    

if __name__ == "__main__":
    _main()