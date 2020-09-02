from argparse import ArgumentParser

from betterproto import Casing
from logzero import logger

from neofox.exceptions import NeofoxInputParametersException
from neofox.model.conversion import ModelConverter
from neofox.neofox import NeoFox
import os
import json


def neofox_cli():
    parser = ArgumentParser(description='adds patient information given in sample file of a cohort to merged icam file')
    parser.add_argument('--input-file', dest='input_file', help='input file with neoantigens')
    parser.add_argument('-i', '--icam-file', dest='icam_file',
                        help='iCaM file with neoantigens (this is an alternative input for iCaM integration)')
    # TODO: once we support the neofox from the models this parameter will not be required
    parser.add_argument('-p', '--patient-id', dest='patient_id', help='the patient id for the iCaM file',
                        required=True)
    parser.add_argument('-d', '--patients-data', dest='patients_data',
                        help='file with data for patients with columns: identifier, estimated_tumor_content, '
                             'is_rna_available, mhc_i_alleles, mhc_i_i_alleles, tissue',
                        required=True)
    parser.add_argument('--output-folder', dest='output_folder', help='output folder', required=True)
    parser.add_argument('--output-prefix', dest='output_prefix',
                        help='prefix to name output files in the output folder', default='neofox')
    parser.add_argument('--with-short-wide-table', dest='with_short_wide_table', action='store_true',
                        help='output results in a short wide tab-separated table '
                             '(if no format is specified this is the default)')
    parser.add_argument('--with-tall-skinny-table', dest='with_tall_skinny_table', action='store_true',
                        help='output results in a tall skinny tab-separated table')
    parser.add_argument('--with-json', dest='with_json', action='store_true',
                        help='output results in JSON format')
    parser.add_argument('--num_cpus', dest='num_cpus', default=1, help='number of CPUs for computation')
    args = parser.parse_args()

    input_file = args.input_file
    icam_file = args.icam_file
    patient_id = args.patient_id
    patients_data = args.patients_data
    output_folder = args.output_folder
    output_prefix = args.output_prefix
    with_sw = args.with_short_wide_table
    with_ts = args.with_tall_skinny_table
    with_json = args.with_json
    num_cpus = int(args.num_cpus)
    if input_file and icam_file:
        raise NeofoxInputParametersException(
            "Please, define either an iCaM file or a standard input file as input. Not both")
    if not input_file and not icam_file:
        raise NeofoxInputParametersException(
            "Please, define one input file, either an iCaM file or a standard input file")
    if not with_sw and not with_ts and not with_json:
        with_sw = True  # if none specified short wide is the default

    # parse the input data
    if icam_file is not None:
        neoantigens = ModelConverter.parse_icam_file(icam_file)
    else:
        neoantigens = ModelConverter.parse_neoantigens_file(input_file)
    patients = ModelConverter.parse_patients_file(patients_data)

    # run annotations
    annotations = NeoFox(
        neoantigens=neoantigens, patients=patients, patient_id=patient_id, num_cpus=num_cpus).get_annotations()

    # writes the output
    if with_sw:
        ModelConverter.annotations2short_wide_table(annotations, neoantigens).to_csv(
            os.path.join(output_folder, "{}.tsv".format(output_prefix)), sep='\t', index=False)
    if with_ts:
        ModelConverter.annotations2tall_skinny_table(annotations).to_csv(
            os.path.join(output_folder, "{}.annotations.tsv".format(output_prefix)), sep='\t', index=False)
        ModelConverter.objects2dataframe(neoantigens).to_csv(
            os.path.join(output_folder, "{}.neoantigens.tsv".format(output_prefix)), sep='\t', index=False)
    if with_json:
        ModelConverter.objects2json(
            annotations, os.path.join(output_folder, "{}.annotations.json".format(output_prefix)))
        ModelConverter.objects2json(
            neoantigens, os.path.join(output_folder, "{}.neoantigens.json".format(output_prefix)))

    logger.info("Finished NeoFox")
