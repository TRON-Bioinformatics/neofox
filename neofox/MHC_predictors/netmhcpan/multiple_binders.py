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
import numpy as np
import scipy.stats as stats
from logzero import logger
from neofox import MHC_I, MHC_II


class MultipleBinding:

    def transform_mhc_prediction_output(self, prediction_out, mhc=MHC_I):
        """
        Takes netmhcpan4 output or netmhcpanII output as neofox (parsed with Netmhc[II]panBestPrediction().filter_binding_predictions) and
        returns tuple of mhc binding rank scores, epitope and HLA allele for all predicted epitopes as list
        """
        pred_data = prediction_out[1]
        list_of_tuples = []
        for ii, i in enumerate(pred_data):
            if mhc == MHC_II:
                # rank, affinity, epitope sequence, allele
                list_of_tuples.append((float(i[9]), float(i[8]), i[2], i[1]))
            else:
                # rank, affinity, epitope sequence, allele
                list_of_tuples.append((float(i[13]), float(i[12]), i[2], i[1]))
        list_of_tuples.sort(key=lambda x: x[0])     # sort by rank
        return list_of_tuples


    def check_for_homozygosity(self, patient_alleles):
        """
        returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles. Otherwise retunrs empty list
        """
        return [allele for allele in patient_alleles if patient_alleles.count(allele) > 1]

    def extract_best_epi_per_alelle(self, tuple_epis, alleles):
        """
        this function returns the predicted epitope with the lowest binding score for each patient allele, considering homozyogosity
        """
        homo_alleles = self.check_for_homozygosity(alleles)
        dict_allels = {}
        for allele in alleles:
            for epi in tuple_epis:
                if allele == epi[-1]:
                    if allele not in dict_allels:
                        dict_allels[allele] = [epi]
                    else:
                        dict_allels[allele].append(epi)
        best_epis_per_allele = []
        for allele in dict_allels:
            dict_allels[allele].sort(key=lambda x: float(x[0]))
            best_epis_per_allele.append(dict_allels[allele][0])
            if allele in homo_alleles:
                # append homozygous allleles two times
                homo_numbers = homo_alleles.count(allele)
                if homo_numbers == 1:
                    best_epis_per_allele.append(dict_allels[allele][0])
                else:
                    homo_best_epi = dict_allels[allele][0]
                    homo_best_epi_all = []
                    # allele already one time represented in list --> add n-t times
                    [homo_best_epi_all.append(tuple(homo_best_epi)) for i in range(homo_numbers - 1)]
                    best_epis_per_allele.extend(tuple(homo_best_epi_all))
        return best_epis_per_allele

    def scores_to_list(self, tuple_epis):
        """
        Takes list of epitope tuple as neofox and returns a list of mhc rank scores of these tuples
        """
        return [epi[0] for epi in tuple_epis]

    def affinities_to_list(self, tuple_epis):
        """
        Takes list of epitope tuple as neofox and returns a list of mhc rank scores of these tuples
        """
        return [epi[1] for epi in tuple_epis]

    def determine_number_of_binders(self, list_scores, threshold=2):
        """
        Determines the number of HLA binders per mutation based on a threshold. Default is set to 2, which is threshold for weak binding using netmhcpan4.
        """
        number_binders = 0
        for score in list_scores:
            if float(score) < threshold:
                number_binders += 1
        number_binders = number_binders if not len(list_scores) == 0 else None
        return number_binders


