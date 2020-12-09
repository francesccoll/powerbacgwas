#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
from cyvcf2 import VCF


# tested with:
#   bcftools version: 1.9 (using htslib 1.9)

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to re-format an input VCF file for the power calculations pipeline. " \
                  "Specifically, it will:\nCheck VCF file is a multi-sample VCF format; split multi-allelic sites; " \
                  "make sure GT genotypes are in haploid format, if not, convert diploid to haploid; add variant ID " \
                  "made up of CHROM.POS.REF.ALT; and select subset of samples, if chosen.\n"
    parser = argparse.ArgumentParser(description=description)
    group = parser.add_argument_group('Required arguments')
    group.add_argument(
        "-v", "--input_vcf", action="store", dest="input_vcf",
        help="multi-sample VCF file",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-o", "--output_vcf", action="store", dest="output_vcf",
        help="output VCF file correctly formatted",
        required=True, metavar="OUTPUT_VCF"
    )
    group = parser.add_argument_group('Optional arguments')
    group.add_argument(
        "-b", "--bcftools_path", action="store", dest="bcftools_path",
        help="Full path to bcftools executable or executable name", required=False, default='bcftools',
        metavar="BCFTOOLS_PATH"
    )
    group.add_argument(
        "-g", "--bgzip_path", action="store", dest="bgzip_path",
        help="Full path to bgzip executable or executable name", required=False, default='bgzip',
        metavar="BGZIP_PATH"
    )
    group.add_argument(
        "-t", "--tabix_path", action="store", dest="tabix_path",
        help="Full path to tabix executable or executable name", required=False, default='tabix',
        metavar="TABIX_PATH"
    )
    group.add_argument(
        "-c", "--chr_id", action="store", dest="chr_id",
        help="Chromosome ID to be used in VCF.\n",
        required=False, metavar="CHR_ID"
    )
    group.add_argument(
        "-f", "--to_haploid_format", action="store_true", dest="to_haploid_format",
        help="Whether to convert VCF GT fields to haploid format (0/0 or 0:. to 0)\n", required=False,
    )
    group.add_argument(
        "-s", "--samples_list", action="store", dest="samples_list",
        help="Samples to keep from VCF file.\nOne sample id per line expected.\n",
        required=False, metavar="SAMPLES_LIST"
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

    # Making sure input VCF file exists
    if not os.path.exists(args.input_vcf):
        logging.error(f'Input VCF file (--input_vcf) {args.input_vcf} not found!')
        sys.exit(-1)

    # Making sure dependencies exist
    check_dependency(args.bcftools_path)

    # Check multi-VCF format
    logging.info(f'Extracting samples from input VCF file {args.input_vcf}')
    stout = run_command_string(''.join([args.bcftools_path, ' query -l ', args.input_vcf]))
    vcf_samples = dict()
    for sample in stout.strip().split('\n'):
        if sample != '':
            vcf_samples[sample] = 1
    if len(vcf_samples) < 2:
        logging.error(f'Only {str(len(vcf_samples))} samples could be extracted from {args.input_vcf}. '
                      f'Use a multi-sample VCF file.')
        sys.exit(-1)
    else:
        logging.info(f'A total of {str(len(vcf_samples))} samples extracted from {args.input_vcf}')

    # Create temporary VCF file
    tmp_vcf_file = args.input_vcf + '_tmp'
    tmp_vcf_file2 = tmp_vcf_file + '_tmp'
    logging.info(f'Creating temporary VCF file {tmp_vcf_file}')
    run_command_shell_string(''.join(['cp ', args.input_vcf, ' ', tmp_vcf_file]))

    # Change CHROM id, if chosen
    if args.chr_id is not None:
        logging.info(f'Replacing VCF chromosome Id')
        stout = run_command_string(''.join([args.bcftools_path, ' view ', args.input_vcf, ' | grep \"^##contig\"']))
        num_contigs = stout.strip().split('\n')
        if len(num_contigs) == 1:
            old_chr_id = stout.strip().split('=')[2].split(',')[0]
            new_chr_id = args.chr_id.strip()
            logging.info(f'Replacing CHROM id \'{old_chr_id}\' with \'{new_chr_id}\'')
            # run_command_shell_string(''.join(['sed -i \'\' \'s/^', old_chr_id, '/', new_chr_id, '/g\'', ' ',
            #                                   tmp_vcf_file]))
            # run_command_shell_string(''.join(['sed -i \'\' \'s/contig=<ID=', old_chr_id, '/contig=<ID=', new_chr_id,
            #                                   '/g\'', ' ', tmp_vcf_file]))
            run_command_shell_string(''.join(['sed -i \'s/^', old_chr_id, '/', new_chr_id, '/g\'', ' ',
                                              tmp_vcf_file]))
            run_command_shell_string(''.join(['sed -i \'s/contig=<ID=', old_chr_id, '/contig=<ID=', new_chr_id,
                                              '/g\'', ' ', tmp_vcf_file]))
        else:
            logging.error(f'Only one contig/chrom supported in input VCF file if option --chr_id is chosen. '
                          f'{str(len(num_contigs))} contigs found.\n')
            sys.exit(-1)

    # Split multi-allelic alleles
    logging.info(f'Split multi-allelic sites')
    run_command_shell_string(''.join([args.bcftools_path, ' norm -Ov -m -any ', tmp_vcf_file, ' > ', tmp_vcf_file2]))
    run_command_shell_string(''.join(['mv ', tmp_vcf_file2, ' ', tmp_vcf_file]))

    # Add variant IDs
    logging.info(f'Adding variant IDs')
    run_command_shell_string(''.join([args.bcftools_path, ' annotate -Ov -x ID -I +\'%CHROM\.%POS\.%REF\.%ALT\' ',
                                      tmp_vcf_file, ' > ', tmp_vcf_file2]))
    run_command_shell_string(''.join(['mv ', tmp_vcf_file2, ' ', tmp_vcf_file]))

    # Check if diploid GT genotypes used, change to haploid
    # logging.info(f'Finding out if diploid or haploid GT format used in {args.input_vcf}')
    # run_command_shell_string(''.join(['cp ', tmp_vcf_file, ' | head -n 1000 > ', tmp_vcf_file2]))
    # stout = run_command_string(''.join([args.bcftools_path, ' query -f \'[%GT ]\' ', tmp_vcf_file2,
    #                                     ' | tr \' \' \'\n\' | sort | uniq']))
    # print(stout)
    # run_command_shell_string(''.join(['rm ', tmp_vcf_file2]))
    # diploid_format = True if '/' in stout else False
    # if ':' in stout:
    #    diploid_format = True
    # if diploid_format:
    if args.to_haploid_format:
        logging.info(f'Converting to haploid format\n')
        # 1. extract columns #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
        vcf_part1 = tmp_vcf_file + '_part1'
        run_command_shell_string(''.join(['cat ', tmp_vcf_file, ' | grep -v "^#" | ',
                                          'cut -f 1-9  > ', vcf_part1]))
        # 2. extract header ##
        vcf_header = tmp_vcf_file + '_header'
        run_command_shell_string(''.join(['cat ', tmp_vcf_file, ' | grep "^#" > ', vcf_header]))
        # 3. extract GT columns
        vcf_part2 = tmp_vcf_file + '_part2'
        run_command_shell_string(''.join([args.bcftools_path, ' query -f \'[\t%GT]\n\' ', tmp_vcf_file, ' > ', vcf_part2]))
        # 4. create new VCF from parts
        vcf_part2_tmp = tmp_vcf_file + '_part2_tmp'
        run_command_shell_string(''.join(['cat ', vcf_part2, ' | sed \'s/0\/0/0/g\' > ', vcf_part2_tmp]))
        run_command_shell_string(''.join(['mv ', vcf_part2_tmp, ' ', vcf_part2]))
        run_command_shell_string(''.join(['cat ', vcf_part2, ' | sed \'s/0\:\./0/g\' > ', vcf_part2_tmp]))
        run_command_shell_string(''.join(['mv ', vcf_part2_tmp, ' ', vcf_part2]))
        run_command_shell_string(''.join(['cat ', vcf_part2, ' | sed \'s/1\/1/1/g\' > ', vcf_part2_tmp]))
        run_command_shell_string(''.join(['mv ', vcf_part2_tmp, ' ', vcf_part2]))
        run_command_shell_string(''.join(['cat ', vcf_part2, ' | sed \'s/1\:\./1/g\' > ', vcf_part2_tmp]))
        run_command_shell_string(''.join(['mv ', vcf_part2_tmp, ' ', vcf_part2]))
        run_command_shell_string(''.join(['cat ', vcf_part2, ' | sed \'s/\.\/\./\./g\' > ', vcf_part2_tmp]))
        run_command_shell_string(''.join(['mv ', vcf_part2_tmp, ' ', vcf_part2]))
        run_command_shell_string(''.join(['cat ', vcf_part2, ' | sed \'s/\.\:\./\./g\' > ', vcf_part2_tmp]))
        run_command_shell_string(''.join(['mv ', vcf_part2_tmp, ' ', vcf_part2]))
        # 4. join new VCF parts back together
        vcf_body = tmp_vcf_file + '_body'
        run_command_shell_string(''.join(['paste -d\'\\0\' ', vcf_part1, ' ', vcf_part2, ' > ', vcf_body]))
        run_command_shell_string(''.join(['cat ', vcf_header, ' ', vcf_body, ' > ', tmp_vcf_file]))
        # remove temporary VCF parts
        run_command_shell_string(''.join(['rm ', vcf_header]))
        run_command_shell_string(''.join(['rm ', vcf_part1]))
        run_command_shell_string(''.join(['rm ', vcf_part2]))
        run_command_shell_string(''.join(['rm ', vcf_body]))

    # gzip and tabix
    run_command_shell_string(''.join(['mv ', tmp_vcf_file, ' ', args.output_vcf]))
    run_command_shell_string(''.join([args.bgzip_path, ' -c ', args.output_vcf , ' > ', args.output_vcf + '.gz']))
    run_command_shell_string(''.join([args.tabix_path, ' -p vcf ', args.output_vcf + '.gz']))


if __name__ == "__main__":
    _main()