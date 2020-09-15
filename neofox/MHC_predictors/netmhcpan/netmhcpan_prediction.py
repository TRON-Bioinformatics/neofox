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

from logzero import logger

from neofox.helpers import data_import
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import AbstractNetMhcPanPredictor


class NetMhcPanPredictor(EpitopeHelper, AbstractNetMhcPanPredictor):

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
        # TODO: was added to netMHCIIpan too --> combine
        if allele.count(":") > 1:
            allele_correct = ":".join(allele.split(":")[0:2])
        else:
            allele_correct = allele
        return allele_correct


    def mhc_prediction(self, hla_alleles, set_available_mhc, tmpfasta, tmppred):
        """ Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file.
        """
        alleles_for_prediction = []
        for allele in hla_alleles:
            allele = self.check_format_allele(allele)
            allele = allele.replace("*", "")
            if allele in set_available_mhc:
                alleles_for_prediction.append(allele)
            else:
                logger.info(allele + "not available")
        hla_allele = ",".join(alleles_for_prediction)
        cmd = [
            self.configuration.net_mhc_pan,
            "-a", hla_allele,
            "-f", tmpfasta,
            "-BA"]
        lines, _ = self.runner.run_command(cmd)
        counter = 0
        # TODO: avoid writing a file here, just return some data structure no need to go to the file system
        with open(tmppred, "w") as f:
            for line in lines.splitlines():
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("#", "-", "HLA", "Prot")):
                        continue
                    if counter == 0 and line.startswith("Pos"):
                        counter += 1
                        line = line.split()
                        line = line[0:-1] if len(line) > 14 else line
                        f.write(";".join(line) + "\n")
                        continue
                    elif counter > 0 and line.startswith("Pos"):
                        continue
                    line = line.split()
                    line = line[0:-2] if len(line) > 14 else line
                    line = ";".join(line)
                    f.write(line + "\n")

    def filter_binding_predictions(self, position_of_mutation, tmppred):
        """filters prediction file for predicted epitopes that cover mutations
        """
        dat_prediction = data_import.import_dat_general(tmppred)
        data_mhc_prediction = dat_prediction[1]
        header = dat_prediction[0]
        data_mhc_prediction_filtered = []
        pos_epi = header.index("Pos")
        epi = header.index("Peptide")
        for ii, i in enumerate(data_mhc_prediction):
            if self.epitope_covers_mutation(position_of_mutation, i[pos_epi], len(i[epi])):
                data_mhc_prediction_filtered.append(data_mhc_prediction[ii])
        return header, data_mhc_prediction_filtered

    def minimal_binding_score(self, prediction_tuple, rank=True):
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used. if rank = False, Aff(nM) is used
        """
        # TODO: generalize this method with netmhcIIpan_prediction.py + change neofox
        header = prediction_tuple[0]
        epitope_data = prediction_tuple[1]
        if rank:
            mhc_score_column = header.index("%Rank")
        else:
            mhc_score_column = header.index("Aff(nM)")
        max_score = float(1000000000000)
        best_predicted_epitope = []
        for ii, i in enumerate(epitope_data):
            mhc_score = float(i[mhc_score_column])
            if mhc_score < max_score:
                max_score = mhc_score
                best_predicted_epitope = i
        return header, best_predicted_epitope


    def filter_for_9mers(self, prediction_tuple):
        """returns only predicted 9mers
        """
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        dat_9mers = []
        for ii, i in enumerate(dat):
            seq = i[seq_col]
            if len(seq) == 9:
                dat_9mers.append(i)
        return dat_head, dat_9mers

    def filter_for_WT_epitope(self, prediction_tuple, mut_seq, mut_allele, number_snv):
        """returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        """
        header = prediction_tuple[0]
        data = prediction_tuple[1]
        seq_col = header.index("Peptide")
        epitopes_wt = []
        for ii, i in enumerate(data):
            wt_seq = i[seq_col]
            if len(wt_seq) == len(mut_seq):
                numb_mismatch = self.hamming_check_0_or_1(mut_seq, wt_seq)
                if numb_mismatch <= number_snv:
                    epitopes_wt.append(i)
        all_epitopes_wt = (header, epitopes_wt)
        self.minimal_binding_score(all_epitopes_wt)
        return self.minimal_binding_score(all_epitopes_wt)

    def filter_for_WT_epitope_position(self, prediction_tuple, sequence_mut, position_epitope):
        """returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        """
        header = prediction_tuple[0]
        data = prediction_tuple[1]
        seq_col = header.index("Peptide")
        pos_col = header.index("Pos")
        epitopes_wt = []
        for ii, i in enumerate(data):
            wt_seq = i[seq_col]
            wt_pos = i[pos_col]
            if (len(wt_seq) == len(sequence_mut)) & (wt_pos == position_epitope):
                epitopes_wt.append(i)
        all_epitopes_wt = (header, epitopes_wt)
        return self.minimal_binding_score(all_epitopes_wt)
