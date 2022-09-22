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
from typing import List
import os
from neofox.helpers import intermediate_files
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Mhc2, Mhc2Name, Mhc2Isoform, PredictedEpitope, Neoantigen
from neofox.references.references import DependenciesConfiguration


class NetMhcIIPanPredictor:

    def __init__(
            self, runner: Runner, configuration: DependenciesConfiguration,
            blastp_runner: BlastpRunner, mhc_parser: MhcParser):

        self.runner = runner
        self.configuration = configuration
        self.mhc_parser = mhc_parser
        self.blastp_runner = blastp_runner

    @staticmethod
    def generate_mhc2_alelle_combinations(mhc_alleles: List[Mhc2]) -> List[Mhc2Isoform]:
        """given list of HLA II alleles, returns list of HLA-DRB1 (2x), all possible HLA-DPA1/HLA-DPB1 (4x)
        and HLA-DQA1/HLA-DPQ1 (4x)
        """
        dp_dq_isoforms = [
            m
            for mhc in mhc_alleles
            if mhc.name in [Mhc2Name.DQ, Mhc2Name.DP]
            for m in mhc.isoforms
        ]
        dr_isoforms = [
            m
            for mhc in mhc_alleles
            if mhc.name == Mhc2Name.DR
            for m in mhc.isoforms
        ]
        mice_isoforms = [
            m
            for mhc in mhc_alleles
            if mhc.name in [Mhc2Name.H2A_molecule, Mhc2Name.H2E_molecule]
            for m in mhc.isoforms
        ]
        return dp_dq_isoforms + dr_isoforms + mice_isoforms

    def represent_mhc2_isoforms(self, isoforms: List[Mhc2Isoform]) -> List[str]:
        return [self.mhc_parser.get_netmhc2pan_representation(i) for i in isoforms]

    def mhc2_prediction(self, mhc_alleles: List[str], sequence) -> List[PredictedEpitope]:
        """ Performs netmhcIIpan prediction for desired hla alleles and writes result to temporary file."""
        # TODO: integrate generate_mhc_ii_alelle_combinations() here to easu utilisation
        tmp_fasta = intermediate_files.create_temp_fasta(
            [sequence], prefix="tmp_singleseq_"
        )
        lines, _ = self.runner.run_command(
            [
                self.configuration.net_mhc2_pan,
                "-BA",
                "-a",
                ",".join(mhc_alleles),
                "-f",
                tmp_fasta
            ]
        )
        os.remove(tmp_fasta)
        return self._parse_netmhcpan_output(lines)
    
    def mhc2_prediction_peptide(
        self, mhc2_isoform: Mhc2Isoform, sequence ) -> PredictedEpitope:
        """ Performs netmhcIIpan prediction for desired hla allele and writes result to temporary file."""
        result = None
        tmp_peptide = intermediate_files.create_temp_peptide([sequence], prefix="tmp_singleseq_")
        lines, _ = self.runner.run_command(
            cmd=[
                self.configuration.net_mhc2_pan,
                "-BA",
                "-a",
                self.represent_mhc2_isoforms([mhc2_isoform])[0],
                "-inptype",
                "1",
                "-f",
                tmp_peptide,
            ],
            print_log=False
        )
        predicted_epitopes = self._parse_netmhcpan_output(lines)
        if predicted_epitopes:
            result = predicted_epitopes[0]
        os.remove(tmp_peptide)
        return result

    def _parse_netmhcpan_output(self, lines: str) -> List[PredictedEpitope]:
        results = []
        for line in lines.splitlines():
            line = line.rstrip().lstrip()
            if line:
                if line.startswith(("#", "-", "Number", "Temporary", "Seq", "ERROR", "Pos")):
                    continue
                line = line.split()
                line = line[0:-1] if len(line) > 12 else line
                results.append(
                    PredictedEpitope(
                        position=int(line[0]),
                        isoform_mhc_i_i=self.mhc_parser.parse_mhc2_isoform(line[1]),
                        mutated_peptide=line[2],
                        affinity_mutated=float(line[11]),
                        rank_mutated=float(line[8]),
                    )
                )
        return results

    def set_wt_netmhcpan_scores(self, predictions) -> List[PredictedEpitope]:
        for p in predictions:
            if p.wild_type_peptide is not None:
                wt_prediction = self.mhc2_prediction_peptide(
                    mhc2_isoform=p.isoform_mhc_i_i,
                    sequence=p.wild_type_peptide)
                if wt_prediction is not None:
                    # NOTE: netmhcpan in peptide mode should return only one epitope
                    p.rank_wild_type = wt_prediction.rank_mutated
                    p.affinity_wild_type = wt_prediction.affinity_mutated
        return predictions

    def get_wt_predictions(self, neoantigen: Neoantigen, patient_mhc2_isoforms):
        predictions = self.mhc2_prediction(patient_mhc2_isoforms, neoantigen.wild_type_xmer)
        predictions = EpitopeHelper.filter_peptides_covering_snv(neoantigen.position, predictions)
        return predictions

    def get_predictions(self, neoantigen: Neoantigen, patient_mhc2_isoforms, uniprot):

        predictions = self.mhc2_prediction(patient_mhc2_isoforms, neoantigen.mutated_xmer)
        if neoantigen.wild_type_xmer:
            # make sure that predicted epitopes cover mutation in case of SNVs
            predictions = EpitopeHelper.filter_peptides_covering_snv(
                position_of_mutation=neoantigen.position, predictions=predictions
            )
        # make sure that predicted neoepitopes are not part of the WT proteome
        filtered_predictions = EpitopeHelper.remove_peptides_in_proteome(predictions, uniprot)
        return filtered_predictions
