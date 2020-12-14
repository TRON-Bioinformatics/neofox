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
import os
from logzero import logger
import neofox
from neofox.exceptions import NeofoxConfigurationException

PREFIX_HOMO_SAPIENS = "homo_sapiens"
HOMO_SAPIENS_FASTA = "Homo_sapiens.fa"
IEDB_FASTA = "IEDB.fasta"
PROTEOME_DB_FOLDER = "proteome_db"
IEDB_FOLDER = "iedb"
IEDB_BLAST_PREFIX = "iedb_blast_db"
NETMHCPAN_AVAILABLE_ALLELES_FILE = "netmhcpan_available_alleles.txt"
NETMHC2PAN_AVAILABLE_ALLELES_FILE = "netmhc2pan_available_alleles.txt"


class AbstractDependenciesConfiguration:
    def _check_and_load_binary(self, variable_name, optional=False):
        variable_value = os.environ.get(variable_name)
        if not optional and variable_value is None:
            raise NeofoxConfigurationException(
                "Please, set the environment variable ${} pointing to the right binary!".format(
                    variable_name
                )
            )
        # checks that the file exists
        if (
            variable_value is not None
        ):  # only optional variables can be None at this stage
            if not os.path.exists(variable_value):
                raise NeofoxConfigurationException(
                    "The provided binary '{}' in ${} does not exist!".format(
                        variable_value, variable_name
                    )
                )
            # checks that it is executable
            if not os.access(variable_value, os.X_OK):
                raise NeofoxConfigurationException(
                    "The provided binary '{}' in ${} is not executable!".format(
                        variable_value, variable_name
                    )
                )
        return variable_value


class DependenciesConfiguration(AbstractDependenciesConfiguration):
    def __init__(self):
        self.blastp = self._check_and_load_binary(neofox.NEOFOX_BLASTP_ENV)
        self.mix_mhc2_pred = self._check_and_load_binary(
            neofox.NEOFOX_MIXMHC2PRED_ENV, optional=True
        )
        if self.mix_mhc2_pred is not None:
            self.mix_mhc2_pred_alleles_list = os.path.join(
                os.path.dirname(self.mix_mhc2_pred), "Alleles_list.txt"
            )
        else:
            self.mix_mhc2_pred_alleles_list = None
        self.mix_mhc_pred = self._check_and_load_binary(
            neofox.NEOFOX_MIXMHCPRED_ENV, optional=True
        )
        if self.mix_mhc_pred is not None:
            self.mix_mhc_pred_alleles_list = os.path.join(
                os.path.dirname(self.mix_mhc_pred), "lib", "allele_list.txt"
            )
        else:
            self.mix_mhc_pred_alleles_list = None
        self.rscript = self._check_and_load_binary(neofox.NEOFOX_RSCRIPT_ENV)
        self.net_mhc2_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHC2PAN_ENV)
        self.net_mhc_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHCPAN_ENV)


class DependenciesConfigurationForInstaller(AbstractDependenciesConfiguration):
    def __init__(self):
        self.net_mhc2_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHC2PAN_ENV)
        self.net_mhc_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHCPAN_ENV)
        self.make_blastdb = self._check_and_load_binary(neofox.NEOFOX_MAKEBLASTDB_ENV)
        self.rscript = self._check_and_load_binary(neofox.NEOFOX_RSCRIPT_ENV)


class ReferenceFolder(object):
    def __init__(self):
        self.reference_genome_folder = self._check_reference_genome_folder()
        # sets the right file names for the resources
        self.available_mhc_ii = self._get_reference_file_name(
            NETMHC2PAN_AVAILABLE_ALLELES_FILE
        )
        self.available_mhc_i = self._get_reference_file_name(
            NETMHCPAN_AVAILABLE_ALLELES_FILE
        )
        self.iedb = self._get_reference_file_name(IEDB_FOLDER)
        self.proteome_db = self._get_reference_file_name(PROTEOME_DB_FOLDER)
        self.uniprot = self._get_reference_file_name(
            os.path.join(PROTEOME_DB_FOLDER, HOMO_SAPIENS_FASTA)
        )

        self.resources = [
            self.available_mhc_ii,
            self.available_mhc_i,
            self.iedb,
            self.proteome_db,
            self.uniprot,
            os.path.join(self.iedb, IEDB_FASTA),
            os.path.join(self.proteome_db, HOMO_SAPIENS_FASTA),
        ]
        self._check_resources()
        self._log_configuration()
        self.__available_alleles = None

    def get_available_alleles(self):
        # this enforces lazy initialisation (useful for testing)
        if not self.__available_alleles:
            self.__available_alleles = AvailableAlleles(self)
        return self.__available_alleles

    def _check_reference_genome_folder(self):
        reference_genome_folder = os.environ.get(neofox.REFERENCE_FOLDER_ENV)
        if reference_genome_folder is None:
            raise NeofoxConfigurationException(
                "Please, set the environment variable ${} pointing to the reference genome folder!".format(
                    neofox.REFERENCE_FOLDER_ENV
                )
            )
        if not os.path.exists(reference_genome_folder):
            raise NeofoxConfigurationException(
                "The provided reference genome '{}' in ${} does not exist!".format(
                    reference_genome_folder, neofox.REFERENCE_FOLDER_ENV
                )
            )
        return reference_genome_folder

    def _check_resources(self):
        # check existence of all resources explicitly defined
        missing_resources = []
        for r in self.resources:
            if not os.path.exists(r):
                missing_resources.append(r)
        if len(missing_resources) > 0:
            raise NeofoxConfigurationException(
                "Missing resources in the reference folder: {}".format(
                    str(missing_resources)
                )
            )

    def _log_configuration(self):
        logger.info("Reference genome folder: {}".format(self.reference_genome_folder))
        logger.info("Resources")
        for r in self.resources:
            logger.info(r)

    def _get_reference_file_name(self, file_name_suffix):
        return os.path.join(self.reference_genome_folder, file_name_suffix)


class AvailableAlleles(object):
    def __init__(self, references):
        self.available_mhc_i = self._load_available_hla_alleles(
            mhc=neofox.MHC_I, references=references
        )
        self.available_mhc_ii = self._load_available_hla_alleles(
            mhc=neofox.MHC_II, references=references
        )

    def _load_available_hla_alleles(
        self, references: ReferenceFolder, mhc=neofox.MHC_I
    ):
        """
        loads file with available hla alllels for netmhcpan4/netmhcIIpan prediction, returns set
        :type references: neofox.references.ReferenceFolder
        :type mhc: str
        :rtype list:
        """
        if mhc == neofox.MHC_II:
            fileMHC = references.available_mhc_ii
        else:
            fileMHC = references.available_mhc_i
        set_available_mhc = set()
        with open(fileMHC) as f:
            for line in f:
                set_available_mhc.add(line.strip())
        return set_available_mhc

    def get_available_mhc_i(self):
        return self.available_mhc_i

    def get_available_mhc_ii(self):
        return self.available_mhc_ii
