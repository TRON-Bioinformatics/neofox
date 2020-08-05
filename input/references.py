import os

from logzero import logger
import pandas as pd

import input
from input.exceptions import INPuTConfigurationException


class DependenciesConfiguration(object):

    def __init__(self):
        self.blastp = self._check_and_load_binary(input.INPUT_BLASTP_ENV)
        self.mix_mhc2_pred = self._check_and_load_binary(input.INPUT_MIXMHC2PRED_ENV)
        self.mix_mhc2_pred_alleles_list = os.path.join(os.path.dirname(self.mix_mhc2_pred), 'Alleles_list.txt')
        self.mix_mhc_pred = self._check_and_load_binary(input.INPUT_MIXMHCPRED_ENV)
        self.rscript = self._check_and_load_binary(input.INPUT_RSCRIPT_ENV)
        self.net_mhc2_pan = self._check_and_load_binary(input.INPUT_NETMHC2PAN_ENV)
        self.net_mhc_pan = self._check_and_load_binary(input.INPUT_NETMHCPAN_ENV)

    @staticmethod
    def _check_and_load_binary(variable_name):
        variable_value = os.environ.get(variable_name, "")
        if not variable_value:
            raise INPuTConfigurationException(
                "Please, set the environment variable ${} pointing to the right binary!".format(
                    variable_name))
        if not os.path.exists(variable_value):
            raise INPuTConfigurationException("The provided binary '{}' in ${} does not exist!".format(
                variable_value, variable_name))
        return variable_value


class ReferenceFolder(object):

    def __init__(self):
        self.reference_genome_folder = self._check_reference_genome_folder()
        # sets the right file names for the resources
        self.available_mhc_ii = self._get_reference_file_name('avail_mhcII.txt')
        self.available_mhc_i = self._get_reference_file_name('MHC_available.csv')
        self.prov_scores_mapped3 = self._get_reference_file_name('PROV_scores_mapped3.csv')
        self.iedb = self._get_reference_file_name('iedb')
        self.proteome_db = self._get_reference_file_name('proteome_db')
        self.uniprot = self._get_reference_file_name('uniprot_human_with_isoforms.fasta')

        self.resources = [
            self.available_mhc_ii, self.available_mhc_i, self.prov_scores_mapped3, self.iedb, self.proteome_db,
            self.uniprot]
        self._check_resources(self.resources)
        self._log_configuration()

    @staticmethod
    def _check_reference_genome_folder():
        reference_genome_folder = os.environ.get(input.REFERENCE_FOLDER_ENV, "")
        if not reference_genome_folder:
            raise INPuTConfigurationException(
                "Please, set the environment variable ${} pointing to the reference genome folder!".format(
                    input.REFERENCE_FOLDER_ENV))
        if not os.path.exists(reference_genome_folder):
            raise INPuTConfigurationException("The provided reference genome '{}' in ${} does not exist!".format(
                reference_genome_folder, input.REFERENCE_FOLDER_ENV))
        return reference_genome_folder

    @staticmethod
    def _check_resources(resources):
        missing_resources = []
        for r in resources:
            if not os.path.exists(r):
                missing_resources.append(r)
        if len(missing_resources) > 0:
            raise INPuTConfigurationException(
                "Missing resources in the reference folder: {}".format(str(missing_resources)))

    def _log_configuration(self):
        logger.info("Reference genome folder: {}".format(self.reference_genome_folder))
        logger.info("Resources")
        for r in self.resources:
            logger.info(r)

    def _get_reference_file_name(self, file_name_suffix):
        return os.path.join(self.reference_genome_folder, file_name_suffix)
