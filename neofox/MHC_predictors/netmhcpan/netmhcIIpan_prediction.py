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
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import Mhc2, MhcAllele, Mhc2Name


class NetMhcIIPanPredictor(AbstractNetMhcPanPredictor):
    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def generate_mhc2_alelle_combinations(self, mhc_alleles: List[Mhc2]) -> List[str]:
        """given list of HLA II alleles, returns list of HLA-DRB1 (2x), all possible HLA-DPA1/HLA-DPB1 (4x)
        and HLA-DQA1/HLA-DPQ1 (4x)
        """
        dp_dq_isoforms = [
            self._represent_dp_and_dq_allele(m.alpha_chain, m.beta_chain)
            for mhc in mhc_alleles
            if mhc.name != Mhc2Name.DR
            for m in mhc.isoforms
        ]
        dr_isoforms = [
            self._represent_drb1_allele(m.beta_chain)
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

    def mhcII_prediction(
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

    def _parse_netmhcpan_output(self, lines: str) -> List[PredictedEpitope]:
        results = []
        for line in lines.splitlines():
            line = line.rstrip().lstrip()
            if line:
                if line.startswith(("#", "-", "Number", "Temporary", "Seq", "ERROR")):
                    continue
                line = line.split()
                line = line[0:-1] if len(line) > 12 else line
                results.append(
                    PredictedEpitope(
                        pos=int(line[0]),
                        hla=line[1],
                        peptide=line[2],
                        affinity_score=float(line[8]),
                        rank=float(line[9]),
                    )
                )
        return results
