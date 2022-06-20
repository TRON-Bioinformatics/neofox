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
import neofox.neofox
from neofox.model.neoantigen import Neoantigen, Patient
from neofox.exceptions import NeofoxInputParametersException
from neofox.neofox import NeoFox
import os
from neofox.model.conversion import ModelConverter
from neofox.references.installer import NeofoxReferenceInstaller
from neofox.references.references import ReferenceFolder, ORGANISM_HOMO_SAPIENS, ORGANISM_MUS_MUSCULUS, MhcDatabase

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
        "--output-folder", dest="output_folder", help="output folder", required=True,
    )
    parser.add_argument(
        "--output-prefix",
        dest="output_prefix",
        help="prefix to name output files in the output folder",
        default="neofox",
    )
    parser.add_argument(
        "--with-table",
        dest="with_table",
        action="store_true",
        help="output results in a short wide tab-separated table "
        "(if no format is specified this is the default)",
    )
    parser.add_argument(
        "--with-json",
        dest="with_json",
        action="store_true",
        help="output results in JSON format",
    )
    parser.add_argument(
        "--with-all-neoepitopes",
        dest="with_all_neoepitopes",
        action="store_true",
        help="output annotations for all MHC-I and MHC-II neoepitopes on all HLA alleles",
    )
    parser.add_argument(
        "--patient-id",
        dest="patient_id",
        help="the patient id for the input file. This parameter is only required, "
        'if the column "patient" has not been added to the candidate file',
    )
    parser.add_argument(
        "--rank-mhci-threshold",
        dest="rank_mhci_threshold",
        help="MHC-I epitopes with a netMHCpan predicted rank greater than or equal than this threshold will be "
             "filtered out",
        default=neofox.RANK_MHCI_THRESHOLD_DEFAULT
    )
    parser.add_argument(
        "--rank-mhcii-threshold",
        dest="rank_mhcii_threshold",
        help="MHC-II epitopes with a netMHCIIpan predicted rank greater than or equal than this threshold will be "
             "filtered out",
        default=neofox.RANK_MHCII_THRESHOLD_DEFAULT
    )
    parser.add_argument(
        "--num-cpus", dest="num_cpus", default=1, help="number of CPUs for computation"
    )
    parser.add_argument(
        "--config",
        dest="config",
        help="an optional configuration file with all the environment variables",
    )
    parser.add_argument(
        "--organism",
        dest="organism",
        choices=[ORGANISM_HOMO_SAPIENS, ORGANISM_MUS_MUSCULUS],
        help="the organism to which the data corresponds",
        default="human"
    )
    parser.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true",
        help="verbose logs",
    )
    args = parser.parse_args()

    candidate_file = args.candidate_file
    json_file = args.json_file
    patient_id = args.patient_id
    patients_data = args.patients_data
    output_folder = args.output_folder
    output_prefix = args.output_prefix
    with_table = args.with_table
    with_json = args.with_json
    with_all_neoepitopes = args.with_all_neoepitopes
    rank_mhci_threshold = float(args.rank_mhci_threshold)
    rank_mhcii_threshold = float(args.rank_mhcii_threshold)
    num_cpus = int(args.num_cpus)
    config = args.config
    organism = args.organism

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
        if not with_table and not with_json:
            with_table = True  # if none specified short wide is the default

        # makes sure that the output folder exists
        os.makedirs(output_folder, exist_ok=True)

        # initialise logs
        log_file_name = NeoFox.get_log_file_name(work_folder=output_folder, output_prefix=output_prefix)
        neofox.neofox.initialise_logs(log_file_name, verbose=args.verbose)

        logger.info("NeoFox v{}".format(neofox.VERSION))

        # loads configuration
        if config:
            dotenv.load_dotenv(config, override=True)
        reference_folder = ReferenceFolder(organism=organism)

        # reads the input data
        neoantigens, patients = _read_data(
            candidate_file,
            json_file,
            patients_data,
            patient_id,
            reference_folder.get_mhc_database())

        # run annotations
        annotated_neoantigens = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            patient_id=patient_id,
            log_file_name=log_file_name,
            num_cpus=num_cpus,
            reference_folder=reference_folder,
            rank_mhci_threshold=rank_mhci_threshold,
            rank_mhcii_threshold=rank_mhcii_threshold,
            with_all_neoepitopes=with_all_neoepitopes
        ).get_annotations()

        _write_results(
            neoantigens=annotated_neoantigens,
            output_folder=output_folder,
            output_prefix=output_prefix,
            with_json=with_json,
            with_table=with_table,
            with_all_neoepitopes=with_all_neoepitopes
        )
    except Exception as e:
        logger.exception(e)  # logs every exception in the file
        raise e

    logger.info("Finished NeoFox")


def _read_data(
    candidate_file, json_file, patients_data, patient_id, mhc_database: MhcDatabase
) -> Tuple[List[Neoantigen], List[Patient]]:
    # parse patient data
    logger.info("Parsing patients data from: {}".format(patients_data))
    patients = ModelConverter.parse_patients_file(patients_data, mhc_database)
    logger.info("Loaded {} patients".format(len(patients)))

    # parse the neoantigen candidate data
    if candidate_file is not None:
        logger.info("Parsing candidate neoantigens from: {}".format(candidate_file))
        neoantigens = ModelConverter.parse_candidate_file(candidate_file, patient_id)
        logger.info("Loaded {} candidate neoantigens".format(len(neoantigens)))
    else:
        logger.info("Parsing candidate neoantigens from: {}".format(json_file))
        neoantigens = ModelConverter.parse_neoantigens_json_file(json_file)
        logger.info("Loaded {} candidate neoantigens".format(len(neoantigens)))

    return neoantigens, patients


def _write_results(neoantigens, output_folder, output_prefix, with_json, with_table, with_all_neoepitopes):
    # NOTE: this import here is a compromise solution so the help of the command line responds faster
    from neofox.model.conversion import ModelConverter
    # writes the output
    if with_table:
        ModelConverter.annotations2neoantigens_table(neoantigens).to_csv(
            os.path.join(
                output_folder,
                "{}_neoantigen_candidates_annotated.tsv".format(output_prefix),
            ),
            sep="\t",
            index=False,
        )

    if with_all_neoepitopes and with_table:
        ModelConverter.annotations2epitopes_table(neoantigens, mhc=neofox.MHC_I).to_csv(
            os.path.join(
                output_folder,
                "{}_mhcI_epitope_candidates_annotated.tsv".format(output_prefix),
            ),
            sep="\t",
            index=False,
        )
        ModelConverter.annotations2epitopes_table(neoantigens, mhc=neofox.MHC_II).to_csv(
            os.path.join(
                output_folder,
                "{}_mhcII_epitope_candidates_annotated.tsv".format(output_prefix),
            ),
            sep="\t",
            index=False,
        )

    if with_json:
        output_features = os.path.join(output_folder, "{}_neoantigen_candidates_annotated.json".format(output_prefix))
        with open(output_features, "wb") as f:
            f.write(json.dumps(ModelConverter.objects2json(neoantigens)))
