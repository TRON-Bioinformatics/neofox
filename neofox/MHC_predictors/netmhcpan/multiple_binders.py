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
from collections import defaultdict
from typing import List, Tuple
from logzero import logger
from neofox.model.neoantigen import MhcAllele, Mhc1, Zygosity, Mhc2, Mhc2Isoform
from neofox.model.conversion import ModelValidator


class MultipleBinding:

    @staticmethod
    def transform_mhc_prediction_output(prediction_out):
        """
        Takes netmhcpan4 output as neofox (parsed with Netmhc[II]panBestPrediction().filter_binding_predictions) and
        returns tuple of mhc binding rank scores, epitope and HLA allele for all predicted epitopes as list
        """
        # TODO: this method is parsing netmhcpan and netmhc2pan output format in a data structure. Move this out to
        #  the NetMhcPan and NetMhc2Pan objects
        list_of_tuples = [(float(i[13]), float(i[12]), i[2], 
                           ModelValidator.validate_mhc_allele_representation(MhcAllele(name=i[1]))) 
                          for i in prediction_out[1]]
        list_of_tuples.sort(key=lambda x: x[0])  # sort by rank
        return list_of_tuples
    
    @staticmethod
    def transform_mhc2_prediction_output(prediction_out):
        """
        Takes netmhcpanII output as neofox (parsed with Netmhc[II]panBestPrediction().filter_binding_predictions) and
        returns tuple of mhc binding rank scores, epitope and HLA allele for all predicted epitopes as list
        """
        # TODO: this method is parsing netmhcpan and netmhc2pan output format in a data structure. Move this out to
        #  the NetMhcPan and NetMhc2Pan objects
        # rank, affinity, epitope sequence, allele combination
        list_of_tuples = [(float(i[9]), float(i[8]), i[2],
                           ModelValidator.validate_mhc2_isoform_representation(Mhc2Isoform(name=i[1])))
                          for i in prediction_out[1]]
        list_of_tuples.sort(key=lambda x: x[0])     # sort by rank
        return list_of_tuples

    @staticmethod
    def _get_homozygous_mhc1_alleles(mhc_isoforms: List[Mhc1]) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [a.name for m in mhc_isoforms for a in m.alleles if m.zygosity == Zygosity.HOMOZYGOUS]

    @staticmethod
    def _get_heterozygous_or_hemizygous_mhc1_alleles(mhc_isoforms: List[Mhc1]) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [a.name for m in mhc_isoforms for a in m.alleles
                if m.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]]

    @staticmethod
    def _get_homozygous_mhc2_alleles(mhc_isoforms: List[Mhc2]) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [a.name for m in mhc_isoforms for g in m.genes for a in g.alleles if g.zygosity == Zygosity.HOMOZYGOUS]

    @staticmethod
    def _get_heterozygous_or_hemizygous_mhc2_alleles(mhc_isoforms: List[Mhc2]) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [a.name for m in mhc_isoforms for g in m.genes for a in g.alleles
                if g.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]]

    @staticmethod
    def extract_best_epitope_per_alelle(tuple_epitopes: List[Tuple], mhc_isoforms: List[Mhc1]):
        """
        This function returns the predicted epitope with the lowest binding score for each patient allele,
        considering homozyogosity
        """
        homozygous_alleles = MultipleBinding._get_homozygous_mhc1_alleles(mhc_isoforms)
        hetero_hemizygous_alleles = MultipleBinding._get_heterozygous_or_hemizygous_mhc1_alleles(mhc_isoforms)
        return MultipleBinding._get_sorted_epitopes(hetero_hemizygous_alleles, homozygous_alleles, tuple_epitopes)

    @staticmethod
    def _get_sorted_epitopes(hetero_hemizygous_alleles, homozygous_alleles, tuple_epitopes):

        # groups epitopes by allele
        epitopes_by_allele = {}
        for epitope in tuple_epitopes:
            allele = epitope[-1].name
            epitopes_by_allele.setdefault(allele, []).append(epitope)
        # chooses the best epitope per allele while considering zygosity
        best_epis_per_allele = []
        for allele, epitopes in epitopes_by_allele.items():
            epitopes.sort(key=lambda x: float(x[0]))    # sort by rank to choose the best epitope
            best_epitope = epitopes[0]
            if best_epitope[-1].name in hetero_hemizygous_alleles:
                best_epis_per_allele.append(best_epitope)  # adds the epitope once
            if best_epitope[-1].name in homozygous_alleles:
                best_epis_per_allele.append(best_epitope)
                best_epis_per_allele.append(best_epitope)  # adds the epitope twice
        return best_epis_per_allele
    
    @staticmethod
    def _get_sorted_epitopes_mhc2(hetero_hemizygous_alleles, homozygous_alleles, tuple_epitopes):

        # groups epitopes by allele
        epitopes_by_allele = {}
        for epitope in tuple_epitopes:
            allele = epitope[-1].name
            epitopes_by_allele.setdefault(allele, []).append(epitope)

        # chooses the best epitope per allele anc considers zygosity
        best_epitopes_per_allele = []
        for allele, epitopes in epitopes_by_allele.items():
            epitopes.sort(key=lambda x: float(x[0]))    # sort by rank to choose the best epitope
            best_epitope = epitopes[0]
            num_repetitions = 0
            if best_epitope[-1].alpha_chain.name in hetero_hemizygous_alleles or \
                    best_epitope[-1].beta_chain.name in hetero_hemizygous_alleles:
                # adds the epitope once if alleles heterozygous
                num_repetitions = 1
            if best_epitope[-1].alpha_chain.name in homozygous_alleles or \
                    best_epitope[-1].beta_chain.name in homozygous_alleles:
                # adds the epitope twice if one allele is homozygous
                num_repetitions = 2
            if best_epitope[-1].alpha_chain.name in homozygous_alleles and \
                    best_epitope[-1].beta_chain.name in homozygous_alleles:
                # adds the epitope four times if both alleles are homozygous
                num_repetitions = 4
            best_epitopes_per_allele.extend([best_epitope for _ in range(num_repetitions)])
        return best_epitopes_per_allele

    @staticmethod
    def extract_best_epitope_per_mhc2_alelle(tuple_epitopes: List[Tuple], mhc_isoforms: List[Mhc2]):
        """
        This function returns the predicted epitope with the lowest binding score for each patient allele,
        considering homozyogosity
        """
        homozygous_alleles = MultipleBinding._get_homozygous_mhc2_alleles(mhc_isoforms)
        hetero_hemizygous_alleles = MultipleBinding._get_heterozygous_or_hemizygous_mhc2_alleles(mhc_isoforms)
        return MultipleBinding._get_sorted_epitopes_mhc2(
            hetero_hemizygous_alleles, homozygous_alleles, tuple_epitopes)

    @staticmethod
    def determine_number_of_binders(list_scores, threshold=2):
        """
        Determines the number of HLA binders per mutation based on a threshold. Default is set to 2, which is threshold for weak binding using netmhcpan4.
        """
        number_binders = 0
        for score in list_scores:
            if float(score) < threshold:
                number_binders += 1
        number_binders = number_binders if not len(list_scores) == 0 else None
        return number_binders


