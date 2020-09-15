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


class DependenciesConfiguration(object):

    def __init__(self):
        self.blastp = self._check_and_load_binary(neofox.NEOFOX_BLASTP_ENV)
        self.mix_mhc2_pred = self._check_and_load_binary(neofox.NEOFOX_MIXMHC2PRED_ENV)
        self.mix_mhc2_pred_alleles_list = os.path.join(os.path.dirname(self.mix_mhc2_pred), 'Alleles_list.txt')
        self.mix_mhc_pred = self._check_and_load_binary(neofox.NEOFOX_MIXMHCPRED_ENV)
        self.rscript = self._check_and_load_binary(neofox.NEOFOX_RSCRIPT_ENV)
        self.net_mhc2_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHC2PAN_ENV)
        self.net_mhc_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHCPAN_ENV)

    @staticmethod
    def _check_and_load_binary(variable_name):
        variable_value = os.environ.get(variable_name, "")
        if not variable_value:
            raise NeofoxConfigurationException(
                "Please, set the environment variable ${} pointing to the right binary!".format(
                    variable_name))
        if not os.path.exists(variable_value):
            raise NeofoxConfigurationException("The provided binary '{}' in ${} does not exist!".format(
                variable_value, variable_name))
        return variable_value


class ReferenceFolder(object):

    def __init__(self):
        self.reference_genome_folder = self._check_reference_genome_folder()
        # sets the right file names for the resources
        self.available_mhc_ii = self._get_reference_file_name('avail_mhcII.txt')
        self.available_mhc_i = self._get_reference_file_name('MHC_available.csv')
        self.iedb = self._get_reference_file_name('iedb')
        self.proteome_db = self._get_reference_file_name('proteome_db')
        self.uniprot = self._get_reference_file_name('uniprot_human_with_isoforms.fasta')

        self.resources = [
            self.available_mhc_ii, self.available_mhc_i, self.iedb, self.proteome_db, self.uniprot]
        self._check_resources(self.resources)
        self._log_configuration()

    @staticmethod
    def _check_reference_genome_folder():
        reference_genome_folder = os.environ.get(neofox.REFERENCE_FOLDER_ENV, "")
        if not reference_genome_folder:
            raise NeofoxConfigurationException(
                "Please, set the environment variable ${} pointing to the reference genome folder!".format(
                    neofox.REFERENCE_FOLDER_ENV))
        if not os.path.exists(reference_genome_folder):
            raise NeofoxConfigurationException("The provided reference genome '{}' in ${} does not exist!".format(
                reference_genome_folder, neofox.REFERENCE_FOLDER_ENV))
        return reference_genome_folder

    @staticmethod
    def _check_resources(resources):
        missing_resources = []
        for r in resources:
            if not os.path.exists(r):
                missing_resources.append(r)
        if len(missing_resources) > 0:
            raise NeofoxConfigurationException(
                "Missing resources in the reference folder: {}".format(str(missing_resources)))

    def _log_configuration(self):
        logger.info("Reference genome folder: {}".format(self.reference_genome_folder))
        logger.info("Resources")
        for r in self.resources:
            logger.info(r)

    def _get_reference_file_name(self, file_name_suffix):
        return os.path.join(self.reference_genome_folder, file_name_suffix)
