#!/usr/bin/env python
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

import tempfile
from typing import List, Set

from logzero import logger

from neofox.helpers import data_import
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import AbstractNetMhcPanPredictor
from neofox.model.neoantigen import HlaAllele
from neofox.model.validation import DPB1, DPA1, DRB1, DQA1, DQB1


class NetMhcIIPanPredictor(EpitopeHelper, AbstractNetMhcPanPredictor):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def generate_mhc_ii_alelle_combinations(self, hla_alleles: List[HlaAllele]) -> List[str]:
        """ given list of HLA II alleles, returns list of HLA-DRB1 (2x), all possible HLA-DPA1/HLA-DPB1 (4x)
        and HLA-DQA1/HLA-DPQ1 (4x)
        """
        drb1_alleles = list(map(
            lambda x: self._represent_drb1_allele(x),
            filter(lambda x: x.gene == DRB1, hla_alleles)))
        dpa1_alleles = list(filter(lambda x: x.gene == DPA1, hla_alleles))
        dpb1_alleles = list(filter(lambda x: x.gene == DPB1, hla_alleles))
        dqa1_alleles = list(filter(lambda x: x.gene == DQA1, hla_alleles))
        dqb1_alleles = list(filter(lambda x: x.gene == DQB1, hla_alleles))
        dp_alleles = [self._represent_dp_and_dq_allele(a, b) for a in dpa1_alleles for b in dpb1_alleles]
        dq_alleles = [self._represent_dp_and_dq_allele(a, b) for a in dqa1_alleles for b in dqb1_alleles]

        return drb1_alleles + dp_alleles + dq_alleles

    @staticmethod
    def _represent_drb1_allele(hla_allele: HlaAllele):
        return "{gene}_{group}{protein}".format(
            gene=hla_allele.gene, group=hla_allele.group, protein=hla_allele.protein)

    @staticmethod
    def _represent_dp_and_dq_allele(hla_a_allele: HlaAllele, hla_b_allele: HlaAllele):
        return "HLA-{gene_a}{group_a}{protein_a}-{gene_b}{group_b}{protein_b}".format(
            gene_a=hla_a_allele.gene, group_a=hla_a_allele.group, protein_a=hla_a_allele.protein,
            gene_b=hla_b_allele.gene, group_b=hla_b_allele.group, protein_b=hla_b_allele.protein)

    def mhcII_prediction(self, hla_alleles: List[str], tmpfasta, tmppred):
        """ Performs netmhcIIpan prediction for desired hla alleles and writes result to temporary file.
        """
        # TODO: integrate generate_mhc_ii_alelle_combinations() here to easu utilisation
        tmp_folder = tempfile.mkdtemp(prefix="tmp_netmhcIIpan_")
        lines, _ = self.runner.run_command([
            self.configuration.net_mhc2_pan,
            "-a", ",".join(hla_alleles),
            "-f", tmpfasta,
            "-tdir", tmp_folder,
            "-dirty"])
        counter = 0
        # TODO: avoid writing a file here, just return some data structure no need to go to the file system
        with open(tmppred, "w") as f:
            for line in lines.splitlines():
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("#", "-", "Number", "Temporary")):
                        continue
                    if counter == 0 and line.startswith("Seq"):
                        counter += 1
                        line = line.split()
                        line = line[0:-1] if len(line) > 12 else line
                        f.write(";".join(line) + "\n")
                        continue
                    elif counter > 0 and line.startswith("Seq"):
                        continue
                    line = line.split()
                    line = line[0:-2] if len(line) > 11 else line
                    f.write(";".join(line) + "\n")

    def filter_binding_predictions(self, position_xmer_sequence, tmppred):
        """
        filters prediction file for predicted epitopes that cover mutations
        """
        header, data = data_import.import_dat_general(tmppred)
        epitopes_covering_mutation = []
        pos_epi = header.index("Seq")
        epi = header.index("Peptide")
        for ii, i in enumerate(data):
            if self.epitope_covers_mutation(position_xmer_sequence, i[pos_epi], len(i[epi])):
                epitopes_covering_mutation.append(data[ii])
        return header, epitopes_covering_mutation

    def minimal_binding_score(self, prediction_tuple, rank=True):
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used.
        if rank = False, Aff(nM) is used
        """
        # TODO: generalize this method with netmhcpan_prediction.py + change neofox
        header = prediction_tuple[0]
        data = prediction_tuple[1]
        if rank:
            mhc_sc = header.index("%Rank")
        else:
            mhc_sc = header.index("Affinity(nM)")
        max_score = float(1000000000)
        best_predicted_epitope = []
        for ii, i in enumerate(data):
            mhc_score = float(i[mhc_sc])
            if mhc_score < max_score:
                max_score = mhc_score
                best_predicted_epitope = i
        return header, best_predicted_epitope

    def filter_for_wt_epitope_position(self, prediction_tuple, mut_seq, position_epitope_in_xmer):
        """returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        """
        header = prediction_tuple[0]
        data = prediction_tuple[1]
        seq_col = header.index("Peptide")
        pos_col = header.index("Seq")
        epitopes_wt = []
        for ii, i in enumerate(data):
            wt_seq = i[seq_col]
            wt_pos = i[pos_col]
            if (len(wt_seq) == len(mut_seq)) & (wt_pos == position_epitope_in_xmer):
                epitopes_wt.append(i)
        all_epitopes_wt = (header, epitopes_wt)
        return self.minimal_binding_score(all_epitopes_wt)

    @staticmethod
    def _get_alleles_netmhcpan_representation(hla_alleles: List[HlaAllele]) -> List[str]:
        return list(map(lambda x: "HLA-{gene}{group}:{protein}".format(gene=x.gene, group=x.group, protein=x.protein),
                        hla_alleles))

    @staticmethod
    def _get_only_available_alleles(hla_alleles: List[str], set_available_mhc: Set[str]) -> str:
        patients_available_alleles = ",".join(list(filter(
            lambda x: x in set_available_mhc,
            AbstractNetMhcPanPredictor._get_alleles_netmhcpan_representation(hla_alleles))))
        patients_not_available_alleles = list(set(hla_alleles).difference(set(set_available_mhc)))
        if len(patients_not_available_alleles) > 0:
            logger.warning("Some alleles from the patient are not available in netmhcpan: {}".format(
                ",".join(patients_not_available_alleles)))
        return patients_available_alleles
