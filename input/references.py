import os
import input
import logging
from input.exceptions import INPuTConfigurationException


class ReferenceFolder(object):

    def __init__(self):
        self.reference_genome_folder = self._check_reference_genome_folder()
        # sets the right file names for the resources
        self.alleles_list_pred = self._get_reference_file_name('Alleles_list_pred2.txt')
        self.available_mhc_ii = self._get_reference_file_name('avail_mhcII.txt')
        self.available_mhc_i = self._get_reference_file_name('MHC_available.csv')
        self.aa_freq_prot = self._get_reference_file_name('20181108_AA_freq_prot.csv')
        self.four_mer_freq = self._get_reference_file_name('20181108_4mer_freq.csv')
        self.aaindex1 = self._get_reference_file_name('aaindex1')
        self.aaindex2 = self._get_reference_file_name('aaindex2')
        self.prov_scores_mapped3 = self._get_reference_file_name('PROV_scores_mapped3.csv')
        self.iedb = self._get_reference_file_name('iedb')
        self.proteome_db = self._get_reference_file_name('proteome_db')
        self.blosum62 = self._get_reference_file_name('BLOSUM62-2.matrix.txt')
        self.tcell_predictor_sir_data = self._get_reference_file_name('SIRdata.mat')
        self.tcell_predictor_gene_expression = self._get_reference_file_name('genes-expression.pickle')
        self.tcell_predictor_aa_features = self._get_reference_file_name('amino-acids-features.pickle')

        # TODO: add this files self.alleles_list_pred, self.avail_mhc_ii
        self.resources = [self.alleles_list_pred, self.available_mhc_ii, self.available_mhc_i, self.aa_freq_prot,
                          self.four_mer_freq, self.aaindex1, self.aaindex2, self.prov_scores_mapped3, self.iedb,
                          self.proteome_db, self.blosum62, self.tcell_predictor_aa_features,
                          self.tcell_predictor_gene_expression, self.tcell_predictor_sir_data]
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
        logging.info("Reference genome folder: {}".format(self.reference_genome_folder))
        logging.info("Resources")
        for r in self.resources:
            logging.info(r)

    def _get_reference_file_name(self, file_name_suffix):
        return os.path.join(self.reference_genome_folder, file_name_suffix)
