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

from logzero import logger

from neofox.helpers import data_import
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import AbstractNetMhcPanPredictor


class NetMhcIIPanPredictor(EpitopeHelper, AbstractNetMhcPanPredictor):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def check_format_allele(self, allele):
        """
        sometimes genotyping may be too detailed. (e.g. HLA-DRB1*04:01:01 should be HLA-DRB1*04:01)
        :param allele: HLA-allele
        :return: HLA-allele in correct format
        """
        if allele.count(":") > 1:
            allele_correct = ":".join(allele.split(":")[0:2])
        else:
            allele_correct = allele
        return allele_correct

    def generate_mhcII_alelles_combination_list(self, hla_alleles, set_available_mhc):
        """ given list of HLA II alleles, returns list of HLA-DRB1 (2x), all possible HLA-DPA1/HLA-DPB1 (4x)
        and HLA-DQA1/HLA-DPQ1 (4x)
        """
        allels_for_prediction = []
        dqa_alleles = []
        dpa_alleles = []
        dqb_alleles = []
        dpb_alleles = []
        for allele in hla_alleles:
            allele = self.check_format_allele(allele)
            if allele.startswith("HLA-DRB1"):
                allele = allele.replace("HLA-", "").replace("*", "_").replace(":", "")
                if allele in set_available_mhc:
                    allels_for_prediction.append(allele)
                else:
                    logger.warn(allele + " not available")
            else:
                allele = allele.replace("*", "").replace(":", "")
                if allele.startswith("HLA-DPA"):
                    dpa_alleles.append(allele)
                elif allele.startswith("HLA-DPB"):
                    dpb_alleles.append(allele)
                elif allele.startswith("HLA-DQA"):
                    dqa_alleles.append(allele)
                elif allele.startswith("HLA-DQB"):
                    dqb_alleles.append(allele)
        dp_alleles = ["-".join([x, y.replace("HLA-", "")]) for x in dpa_alleles for y in dpb_alleles]
        dq_alleles = ["-".join([x, y.replace("HLA-", "")]) for x in dqa_alleles for y in dqb_alleles]
        dp_dq_alleles = dp_alleles + dq_alleles
        for allele in dp_dq_alleles:
            if allele in set_available_mhc:
                allels_for_prediction.append(allele)
        return allels_for_prediction

    def mhcII_prediction(self, hla_alleles, set_available_mhc, tmpfasta, tmppred):
        """ Performs netmhcIIpan prediction for desired hla alleles and writes result to temporary file.
        """
        alleles_for_prediction = self.generate_mhcII_alelles_combination_list(hla_alleles, set_available_mhc)
        hla_allele = ",".join(alleles_for_prediction)
        tmp_folder = tempfile.mkdtemp(prefix="tmp_netmhcIIpan_")
        lines, _ = self.runner.run_command([
            self.configuration.net_mhc2_pan,
            "-a", hla_allele,
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
