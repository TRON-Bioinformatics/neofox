#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
from argparse import ArgumentParser
from logzero import logger
from neofox.exceptions import NeofoxInputParametersException
from neofox.neofox import NeoFox
import os


def neofox_cli():
    parser = ArgumentParser(description='adds patient information given in sample file of a cohort to merged icam file')
    parser.add_argument('--model-file', dest='model_file',
                        help='input tabular file with neoantigens')
    parser.add_argument('--icam-file', dest='icam_file',
                        help='input iCaM file with neoantigens (this is an alternative input for iCaM integration)')
    # TODO: once we support the neofox from the models this parameter will not be required
    parser.add_argument('--patient-id', dest='patient_id', help='the patient id for the iCaM file',
                        required=True)
    parser.add_argument('--patient-data', dest='patients_data',
                        help='file with data for patients with columns: identifier, estimated_tumor_content, '
                             'is_rna_available, mhc_i_alleles, mhc_ii_alleles, tissue',
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
    parser.add_argument('--num-cpus', dest='num_cpus', default=1, help='number of CPUs for computation')
    args = parser.parse_args()

    model_file = args.model_file
    icam_file = args.icam_file
    patient_id = args.patient_id
    patients_data = args.patients_data
    output_folder = args.output_folder
    output_prefix = args.output_prefix
    with_sw = args.with_short_wide_table
    with_ts = args.with_tall_skinny_table
    with_json = args.with_json
    num_cpus = int(args.num_cpus)
    if model_file and icam_file:
        raise NeofoxInputParametersException(
            "Please, define either an iCaM file or a standard input file as input. Not both")
    if not model_file and not icam_file:
        raise NeofoxInputParametersException(
            "Please, define one input file, either an iCaM file or a standard input file")
    if not with_sw and not with_ts and not with_json:
        with_sw = True  # if none specified short wide is the default

    # makes sure that the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    neoantigens, patients = _read_data(icam_file, model_file, patients_data)

    # run annotations

    annotations = NeoFox(neoantigens=neoantigens, patients=patients, patient_id=patient_id, work_folder=output_folder,
                         output_prefix = output_prefix, num_cpus=num_cpus
                         ).get_annotations()

    _write_results(annotations, neoantigens, output_folder, output_prefix, with_json, with_sw, with_ts)

    logger.info("Finished NeoFox")


def _read_data(icam_file, model_file, patients_data):
    # NOTE: this import here is a compromise solution so the help of the command line responds faster
    from neofox.model.conversion import ModelConverter
    # parse the input data
    if icam_file is not None:
        neoantigens = ModelConverter.parse_icam_file(icam_file)
    else:
        neoantigens = ModelConverter.parse_neoantigens_file(model_file)
    patients = ModelConverter.parse_patients_file(patients_data)
    return neoantigens, patients

def _write_results(annotations, neoantigens, output_folder, output_prefix, with_json, with_sw, with_ts):
    # NOTE: this import here is a compromise solution so the help of the command line responds faster
    from neofox.model.conversion import ModelConverter
    # writes the output
    if with_sw:
        ModelConverter.annotations2short_wide_table(annotations, neoantigens).to_csv(
            os.path.join(output_folder, "{}_neoantigens_features_short_wide.tsv".format(output_prefix)), sep='\t',
            index=False)
    if with_ts:
        ModelConverter.annotations2tall_skinny_table(annotations).to_csv(
            os.path.join(output_folder, "{}_features_tall_skinny.tsv".format(output_prefix)), sep='\t', index=False)
        ModelConverter.objects2dataframe(neoantigens).to_csv(
            os.path.join(output_folder, "{}_neoantigens.tsv".format(output_prefix)), sep='\t', index=False)
    if with_json:
        ModelConverter.objects2json(
            annotations, os.path.join(output_folder, "{}_features.json".format(output_prefix)))
        ModelConverter.objects2json(
            neoantigens, os.path.join(output_folder, "{}_neoantigens.json".format(output_prefix)))
