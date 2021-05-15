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
from neofox.helpers import intermediate_files
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import (
    AbstractNetMhcPanPredictor,
    PredictedEpitope,
)
from neofox.model.neoantigen import Mhc2, MhcAllele, Mhc2Name, Mhc2Isoform, Mhc2GeneName


class NetMhcIIPanPredictor(AbstractNetMhcPanPredictor):

    @staticmethod
    def generate_mhc2_alelle_combinations(mhc_alleles: List[Mhc2]) -> List[Mhc2Isoform]:
        """given list of HLA II alleles, returns list of HLA-DRB1 (2x), all possible HLA-DPA1/HLA-DPB1 (4x)
        and HLA-DQA1/HLA-DPQ1 (4x)
        """
        dp_dq_isoforms = [
            m
            for mhc in mhc_alleles
            if mhc.name != Mhc2Name.DR
            for m in mhc.isoforms
        ]
        dr_isoforms = [
            m
            for mhc in mhc_alleles
            if mhc.name == Mhc2Name.DR
            for m in mhc.isoforms
        ]
        return dp_dq_isoforms + dr_isoforms

    @staticmethod
    def _represent_drb1_allele(mhc_allele: MhcAllele):
        return "{gene}_{group}{protein}".format(
            gene=mhc_allele.gene, group=mhc_allele.group, protein=mhc_allele.protein
        )

    @staticmethod
    def _represent_dp_and_dq_allele(mhc_a_allele: MhcAllele, mhc_b_allele: MhcAllele):
        return "HLA-{gene_a}{group_a}{protein_a}-{gene_b}{group_b}{protein_b}".format(
            gene_a=mhc_a_allele.gene,
            group_a=mhc_a_allele.group,
            protein_a=mhc_a_allele.protein,
            gene_b=mhc_b_allele.gene,
            group_b=mhc_b_allele.group,
            protein_b=mhc_b_allele.protein,
        )

    @staticmethod
    def represent_mhc2_isoforms(isoforms: List[Mhc2Isoform]) -> List[str]:
        return [
            NetMhcIIPanPredictor._represent_drb1_allele(i.beta_chain)
            if i.beta_chain.gene == Mhc2GeneName.DRB1.name
            else NetMhcIIPanPredictor._represent_dp_and_dq_allele(i.alpha_chain, i.beta_chain)
            for i in isoforms
        ]

    def mhc2_prediction(
        self, mhc_alleles: List[str], sequence
    ) -> List[PredictedEpitope]:
        """ Performs netmhcIIpan prediction for desired hla alleles and writes result to temporary file."""
        # TODO: integrate generate_mhc_ii_alelle_combinations() here to easu utilisation
        tmp_fasta = intermediate_files.create_temp_fasta(
            [sequence], prefix="tmp_singleseq_"
        )
        tmp_folder = tempfile.mkdtemp(prefix="tmp_netmhcIIpan_")
        lines, _ = self.runner.run_command(
            [
                self.configuration.net_mhc2_pan,
                "-BA",
                "-a",
                ",".join(mhc_alleles),
                "-f",
                tmp_fasta,
                "-tdir",
                tmp_folder,
                "-dirty",
            ]
        )
        return self._parse_netmhcpan_output(lines)
    
    def mhc2_prediction_peptide(
        self, mhc2_isoform: Mhc2Isoform, sequence ) -> List[PredictedEpitope]:
        """ Performs netmhcIIpan prediction for desired hla allele and writes result to temporary file."""
        tmp_peptide = intermediate_files.create_temp_peptide(
            [sequence], prefix="tmp_singleseq_"
        )
        tmp_folder = tempfile.mkdtemp(prefix="tmp_netmhcIIpan_")
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
                "-tdir",
                tmp_folder,
                "-dirty",
            ],
            print_log=False
        )
        return self._parse_netmhcpan_output(lines)

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
                        pos=int(line[0]),
                        hla=self.mhc_parser.parse_mhc2_isoform(line[1]),
                        peptide=line[2],
                        affinity_score=float(line[11]),
                        rank=float(line[8]),
                    )
                )
        return results
