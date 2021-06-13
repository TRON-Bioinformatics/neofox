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
from typing import Tuple, List
import dotenv
from logzero import logger
import orjson as json
import neofox
from neofox.model.neoantigen import Neoantigen, Patient, NeoantigenAnnotations
from neofox.exceptions import NeofoxInputParametersException
from neofox.neofox import NeoFox, AFFINITY_THRESHOLD_DEFAULT
import os
from neofox.model.conversion import ModelConverter
from neofox.references.installer import NeofoxReferenceInstaller
from neofox.references.references import ReferenceFolder, HlaDatabase

epilog = "NeoFox (NEOantigen Feature toolbOX) {}. Copyright (c) 2020-2021 " \
         "TRON – Translational Oncology at the University Medical Center of the " \
         "Johannes Gutenberg University Mainz gGmbH, all rights reserved".format(neofox.VERSION)


def neofox_configure():
    parser = ArgumentParser(description="NeoFox references installer", epilog=epilog)
    parser.add_argument(
        "--reference-folder",
        dest="reference_folder",
        help="the folder with the references required for Neofox",
        required=True,
    )
    parser.add_argument(
        "--install-r-dependencies",
        dest="install_r_dependencies",
        action="store_true",
        help="install the R dependencies automatically",
    )

    args = parser.parse_args()
    reference_folder = args.reference_folder
    install_r_dependencies = args.install_r_dependencies

    # makes sure that the output folder exists
    os.makedirs(reference_folder, exist_ok=True)

    logger.info("Starting the installation of references")
    NeofoxReferenceInstaller(
        reference_folder=reference_folder, install_r_dependencies=install_r_dependencies
    ).install()
    logger.info("Finished the installation succesfully!")


def neofox_cli():
    parser = ArgumentParser(
        description="NeoFox {} annotates a given set of neoantigen candidate sequences "
                    "derived from point mutation with relevant neoantigen features".format(neofox.VERSION),
        epilog=epilog
    )
    parser.add_argument(
        "--candidate-file",
        dest="candidate_file",
        help="input file with neoantigens candidates represented by long mutated peptide sequences",
    )
    parser.add_argument(
        "--json-file",
        dest="json_file",
        help="input JSON file with neoantigens candidates represented by long mutated peptide sequences",
    )
    parser.add_argument(
        "--patient-data",
        dest="patients_data",
        help="file with data for patients with columns: identifier, estimated_tumor_content, "
        "mhc_i_alleles, mhc_ii_alleles, tissue",
        required=True,
    )
    parser.add_argument(
        "--output-folder", dest="output_folder", help="output folder", required=True
    )
    parser.add_argument(
        "--output-prefix",
        dest="output_prefix",
        help="prefix to name output files in the output folder",
        default="neofox",
    )
    parser.add_argument(
        "--with-short-wide-table",
        dest="with_short_wide_table",
        action="store_true",
        help="output results in a short wide tab-separated table "
        "(if no format is specified this is the default)",
    )
    parser.add_argument(
        "--with-tall-skinny-table",
        dest="with_tall_skinny_table",
        action="store_true",
        help="output results in a tall skinny tab-separated table",
    )
    parser.add_argument(
        "--with-json",
        dest="with_json",
        action="store_true",
        help="output results in JSON format",
    )
    parser.add_argument(
        "--patient-id",
        dest="patient_id",
        help="the patient id for the input file. This parameter is only required, "
        'if the column "patient" has not been added to the candidate file',
    )
    parser.add_argument(
        "--affinity-threshold",
        dest="affinity_threshold",
        help="neoantigen candidates with a best predicted affinity greater than or equal than this threshold will be "
             "not annotated with features that specifically model neoepitope recognition. A threshold that is commonly "
             "used is 500 nM",
        default=AFFINITY_THRESHOLD_DEFAULT
    )
    parser.add_argument(
        "--num-cpus", dest="num_cpus", default=1, help="number of CPUs for computation"
    )
    parser.add_argument(
        "--config",
        dest="config",
        help="an optional configuration file with all the environment variables",
    )
    args = parser.parse_args()

    candidate_file = args.candidate_file
    json_file = args.json_file
    patient_id = args.patient_id
    patients_data = args.patients_data
    output_folder = args.output_folder
    output_prefix = args.output_prefix
    with_sw = args.with_short_wide_table
    with_ts = args.with_tall_skinny_table
    with_json = args.with_json
    affinity_threshold = int(args.affinity_threshold)
    num_cpus = int(args.num_cpus)
    config = args.config

    try:
        # check parameters
        if bool(candidate_file) + bool(json_file) > 1:
            raise NeofoxInputParametersException(
                "Please, define either a candidate file, a standard input file or a JSON file as input. Not many of them"
            )
        if not candidate_file and not json_file:
            raise NeofoxInputParametersException(
                "Please, define one input file, either a candidate file, a standard input file or a JSON file"
            )
        if not with_sw and not with_ts and not with_json:
            with_sw = True  # if none specified short wide is the default

        # makes sure that the output folder exists
        os.makedirs(output_folder, exist_ok=True)

        # loads configuration
        if config:
            dotenv.load_dotenv(config, override=True)
        reference_folder = ReferenceFolder()

        # reads the input data
        neoantigens, patients, external_annotations = _read_data(
            candidate_file, json_file, patients_data, patient_id, reference_folder.get_hla_database()
        )

        # run annotations
        annotations = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            patient_id=patient_id,
            work_folder=output_folder,
            output_prefix=output_prefix,
            num_cpus=num_cpus,
            reference_folder=reference_folder,
            affinity_threshold=affinity_threshold
        ).get_annotations()

        # combine neoantigen feature annotations and potential user-specific external annotation
        neoantigen_annotations = _combine_features_with_external_annotations(annotations, external_annotations)

        _write_results(
            neoantigen_annotations,
            neoantigens,
            output_folder,
            output_prefix,
            with_json,
            with_sw,
            with_ts,
        )
    except Exception as e:
        logger.exception(e)  # logs every exception in the file
        raise e

    logger.info("Finished NeoFox")


def _read_data(
    candidate_file, json_file, patients_data, patient_id, hla_database: HlaDatabase
) -> Tuple[List[Neoantigen], List[Patient], List[NeoantigenAnnotations]]:
    # parse patient data
    patients = ModelConverter.parse_patients_file(patients_data, hla_database)
    # parse the neoantigen candidate data
    if candidate_file is not None:
        neoantigens, external_annotations = ModelConverter.parse_candidate_file(
            candidate_file, patient_id
        )
    else:
        neoantigens = ModelConverter.parse_neoantigens_json_file(json_file)
        external_annotations = []

    return neoantigens, patients, external_annotations


def _write_results(
    annotations,  neoantigens, output_folder, output_prefix, with_json, with_sw, with_ts
):
    # NOTE: this import here is a compromise solution so the help of the command line responds faster
    from neofox.model.conversion import ModelConverter

    # writes the output
    if with_sw:
        ModelConverter.annotations2short_wide_table(annotations, neoantigens).to_csv(
            os.path.join(
                output_folder,
                "{}_neoantigen_candidates_annotated.tsv".format(output_prefix),
            ),
            sep="\t",
            index=False,
        )
    if with_ts:
        ModelConverter.annotations2tall_skinny_table(annotations).to_csv(
            os.path.join(
                output_folder, "{}_neoantigen_features.tsv".format(output_prefix)
            ),
            sep="\t",
            index=False,
        )
        ModelConverter.neoantigens2table(neoantigens).to_csv(
            os.path.join(
                output_folder, "{}_neoantigen_candidates.tsv".format(output_prefix)
            ),
            sep="\t",
            index=False,
        )
    if with_json:
        output_features = os.path.join(output_folder, "{}_neoantigen_features.json".format(output_prefix))
        with open(output_features, "wb") as f:
            f.write(json.dumps(ModelConverter.objects2json(annotations)))

        output_neoantigens = os.path.join(output_folder, "{}_neoantigen_candidates.json".format(output_prefix))
        with open(output_neoantigens, "wb") as f:
            f.write(json.dumps(ModelConverter.objects2json(neoantigens)))


def _combine_features_with_external_annotations(
        annotations: List[NeoantigenAnnotations],
        external_annotations: List[NeoantigenAnnotations],
) -> List[NeoantigenAnnotations]:
    final_annotations = []
    for annotation in annotations:
        for annotation_extern in external_annotations:
            if (
                    annotation.neoantigen_identifier
                    == annotation_extern.neoantigen_identifier
            ):
                annotation.annotations = (
                        annotation.annotations + annotation_extern.annotations
                )
        final_annotations.append(annotation)
    return final_annotations
