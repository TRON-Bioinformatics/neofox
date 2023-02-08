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
from typing import Tuple, List, Dict
import dotenv
from logzero import logger
import orjson as json
import neofox
import neofox.neofox
from neofox.model.neoantigen import Neoantigen, Patient, PredictedEpitope
from neofox.model.validation import ModelValidator
from neofox.neofox import NeoFox
import os
from neofox.model.conversion import ModelConverter
from neofox.neofox_epitope import NeoFoxEpitope
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
        "--input-file",
        dest="input_file",
        help="Input file with neoantigens candidates represented by long mutated peptide sequences. "
             "Supported formats: tab-separated columns (extensions: .txt or .tsv) or JSON (extension: .json)",
        required=True,
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
        "--with-all-neoepitopes",
        dest="with_all_neoepitopes",
        action="store_true",
        help="output annotations for all MHC-I and MHC-II neoepitopes on all HLA alleles",
    )
    parser.add_argument(
        "--rank-mhci-threshold",
        dest="rank_mhci_threshold",
        help="MHC-I epitopes with a netMHCpan predicted rank greater than or equal than this threshold will be "
             "filtered out (default: {})".format(neofox.RANK_MHCI_THRESHOLD_DEFAULT),
        default=neofox.RANK_MHCI_THRESHOLD_DEFAULT
    )
    parser.add_argument(
        "--rank-mhcii-threshold",
        dest="rank_mhcii_threshold",
        help="MHC-II epitopes with a netMHCIIpan predicted rank greater than or equal than this threshold will be "
             "filtered out (default: {})".format(neofox.RANK_MHCII_THRESHOLD_DEFAULT),
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

    input_file = args.input_file
    patients_data = args.patients_data
    output_folder = args.output_folder
    output_prefix = args.output_prefix
    with_all_neoepitopes = args.with_all_neoepitopes
    rank_mhci_threshold = float(args.rank_mhci_threshold)
    rank_mhcii_threshold = float(args.rank_mhcii_threshold)
    num_cpus = int(args.num_cpus)
    config = args.config
    organism = args.organism

    try:
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
            input_file,
            patients_data,
            reference_folder.get_mhc_database())

        # run annotations
        annotated_neoantigens = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
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
            with_all_neoepitopes=with_all_neoepitopes
        )
    except Exception as e:
        logger.exception(e)  # logs every exception in the file
        raise e

    logger.info("Finished NeoFox")


def _read_data(input_file, patients_data, mhc_database: MhcDatabase) -> Tuple[List[Neoantigen], List[Patient]]:
    # parse patient data
    logger.info("Parsing patients data from: {}".format(patients_data))
    patients = ModelConverter.parse_patients_file(patients_data, mhc_database)
    logger.info("Loaded {} patients".format(len(patients)))

    # parse the neoantigen candidate data
    if input_file.endswith('.txt') or input_file.endswith('.tsv'):
        logger.info("Parsing candidate neoantigens from: {}".format(input_file))
        neoantigens = ModelConverter.parse_candidate_file(input_file)
        logger.info("Loaded {} candidate neoantigens".format(len(neoantigens)))
    elif input_file.endswith('.json')  :
        logger.info("Parsing candidate neoantigens from: {}".format(input_file))
        neoantigens = ModelConverter.parse_neoantigens_json_file(input_file)
        logger.info("Loaded {} candidate neoantigens".format(len(neoantigens)))
    else:
        raise ValueError('Not supported input file extension: {}'.format(input_file))

    return neoantigens, patients


def _write_results(neoantigens, output_folder, output_prefix, with_all_neoepitopes):
    # NOTE: this import here is a compromise solution so the help of the command line responds faster
    from neofox.model.conversion import ModelConverter
    # writes the output
    ModelConverter.annotations2neoantigens_table(neoantigens).to_csv(
        os.path.join(
            output_folder,
            "{}_neoantigen_candidates_annotated.tsv".format(output_prefix),
        ),
        sep="\t",
        index=False,
    )

    if with_all_neoepitopes:
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

    output_features = os.path.join(output_folder, "{}_neoantigen_candidates_annotated.json".format(output_prefix))
    with open(output_features, "wb") as f:
        f.write(json.dumps(ModelConverter.objects2json(neoantigens)))


def neofox_epitope_cli():
    parser = ArgumentParser(
        description="NeoFox {} epitope annotates a given set of neoepitope candidates "
                    "derived from point mutation with relevant immunogenic features".format(neofox.VERSION),
        epilog=epilog
    )
    parser.add_argument(
        "--input-file",
        dest="input_file",
        help="Input file with neoepitope candidates. "
             "Supported formats: tab-separated columns (extensions: .txt or .tsv) or JSON (extension: .json)",
        required=True,
    )
    parser.add_argument(
        "--patient-data",
        dest="patients_data",
        help="file with data for patients with columns: identifier, estimated_tumor_content, "
        "mhc_i_alleles, mhc_ii_alleles, tissue",
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

    input_file = args.input_file
    patients_data = args.patients_data
    output_folder = args.output_folder
    output_prefix = args.output_prefix
    num_cpus = int(args.num_cpus)
    config = args.config
    organism = args.organism

    try:
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
        neoepitopes, patients = _read_data_epitopes(
            input_file,
            patients_data,
            reference_folder.get_mhc_database())

        # run annotations
        annotated_neoepitopes = NeoFoxEpitope(
            neoepitopes=neoepitopes,
            patients=patients,
            log_file_name=log_file_name,
            num_cpus=num_cpus,
            reference_folder=reference_folder
        ).get_annotations()

        _write_results_epitopes(
            neoepitopes=annotated_neoepitopes,
            output_folder=output_folder,
            output_prefix=output_prefix
        )
    except Exception as e:
        logger.exception(e)  # logs every exception in the file
        raise e

    logger.info("Finished NeoFox epitopes")


def _read_data_epitopes(
    input_file, patients_data, mhc_database: MhcDatabase) -> Tuple[List[PredictedEpitope], List[Patient]]:

    # parse patient data
    patients = []
    if patients_data is not None:
        logger.info("Parsing patients data from: {}".format(patients_data))
        patients = ModelConverter.parse_patients_file(patients_data, mhc_database)
        logger.info("Loaded {} patients".format(len(patients)))

    # parse the neoantigen candidate data
    if input_file.endswith('.txt') or input_file.endswith('.tsv'):
        logger.info("Parsing candidate neoepitopes from: {}".format(input_file))
        neoepitopes = ModelConverter.parse_candidate_neoepitopes_file(input_file, mhc_database)
        logger.info("Loaded {} candidate neoepitopes".format(len(neoepitopes)))
    # TODO: add support for input in JSON format
    #elif input_file.endswith('.json')  :
    #    logger.info("Parsing candidate neoepitopes from: {}".format(input_file))
    #    neoepitopes = ModelConverter.parse_neoepitopes_json_file(input_file)
    #    logger.info("Loaded {} candidate neoepitopes".format(len(neoepitopes)))
    else:
        raise ValueError('Not supported input file extension: {}'.format(input_file))

    return neoepitopes, patients


def _write_results_epitopes(neoepitopes: List[PredictedEpitope], output_folder, output_prefix):
    # NOTE: this import here is a compromise solution so the help of the command line responds faster
    from neofox.model.conversion import ModelConverter

    mhci_neoepitopes = [n for n in neoepitopes if ModelValidator.is_mhci_epitope(n)]
    mhcii_neoepitopes = [n for n in neoepitopes if ModelValidator.is_mhcii_epitope(n)]

    if mhci_neoepitopes:
        ModelConverter.annotated_neoepitopes2epitopes_table(mhci_neoepitopes, mhc=neofox.MHC_I).to_csv(
            os.path.join(
                output_folder,
                "{}_mhcI_epitope_candidates_annotated.tsv".format(output_prefix),
            ),
            sep="\t",
            index=False,
        )
    if mhcii_neoepitopes:
        ModelConverter.annotated_neoepitopes2epitopes_table(mhcii_neoepitopes, mhc=neofox.MHC_II).to_csv(
            os.path.join(
                output_folder,
                "{}_mhcII_epitope_candidates_annotated.tsv".format(output_prefix),
            ),
            sep="\t",
            index=False,
        )

    if neoepitopes:
        output_features = os.path.join(output_folder, "{}_neoepitope_candidates_annotated.json".format(output_prefix))
        with open(output_features, "wb") as f:
            f.write(json.dumps(ModelConverter.objects2json(neoepitopes)))
