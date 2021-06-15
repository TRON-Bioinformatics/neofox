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
from typing import List, Set
from logzero import logger
from neofox.helpers import intermediate_files
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import (
    AbstractNetMhcPanPredictor,
    PredictedEpitope,
)
from neofox.model.neoantigen import Mhc1


class NetMhcPanPredictor(AbstractNetMhcPanPredictor):

    def mhc_prediction(
            self, mhc_alleles: List[Mhc1], set_available_mhc: Set, sequence

    ) -> List[PredictedEpitope]:
        """Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file."""
        input_fasta = intermediate_files.create_temp_fasta(
            sequences=[sequence], prefix="tmp_singleseq_"
        )
        cmd = [
            self.configuration.net_mhc_pan,
            "-a",
            self._get_only_available_alleles(mhc_alleles, set_available_mhc),
            "-f",
            input_fasta,
            "-BA",
        ]
        lines, _ = self.runner.run_command(cmd)
        return self._parse_netmhcpan_output(lines)

    def mhc_prediction_peptide(self, mhc_alleles: List[Mhc1], set_available_mhc: Set, sequence
    ) -> List[PredictedEpitope]:
        """Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file."""
        input_peptide = intermediate_files.create_temp_peptide(
            sequences=[sequence], prefix="tmp_singleseq_"
        )
        cmd = [
            self.configuration.net_mhc_pan,
            "-a",
            self._get_only_available_alleles(mhc_alleles, set_available_mhc),
            "-p",
            input_peptide,
            "-BA",
        ]
        lines, _ = self.runner.run_command(cmd, print_log=False)
        return self._parse_netmhcpan_output(lines)

    def _parse_netmhcpan_output(self, lines: str) -> List[PredictedEpitope]:
        results = []
        for line in lines.splitlines():
            line = line.rstrip().lstrip()
            if line:
                if line.startswith(("#", "-", "HLA", "Prot", "Pos", "No")):
                    continue
                line = line.split()
                line = line[0:-2] if len(line) > 16 else line
                results.append(
                    PredictedEpitope(
                        pos=int(line[0]),
                        hla=self.mhc_parser.parse_mhc_allele(line[1]),
                        peptide=line[2],
                        affinity_score=float(line[15]),
                        rank=float(line[12]),
                    )
                )
        return results

    @staticmethod
    def get_alleles_netmhcpan_representation(mhc_isoforms: List[Mhc1]) -> List[str]:
        return list(
            map(
                lambda x: "HLA-{gene}{group}:{protein}".format(
                    gene=x.gene, group=x.group, protein=x.protein
                ),
                [a for m in mhc_isoforms for a in m.alleles],
            )
        )

    @staticmethod
    def _get_only_available_alleles(
        mhc_isoforms: List[Mhc1], set_available_mhc: Set[str]
    ) -> str:
        hla_alleles_names = NetMhcPanPredictor.get_alleles_netmhcpan_representation(
            mhc_isoforms
        )
        patients_available_alleles = ",".join(
            list(filter(lambda x: x in set_available_mhc, hla_alleles_names))
        )
        patients_not_available_alleles = list(
            set(hla_alleles_names).difference(set(set_available_mhc))
        )
        if len(patients_not_available_alleles) > 0:
            logger.warning(
                "MHC I alleles {} are not supported by NetMHCpan and no binding or derived features will "
                "include it".format(",".join(patients_not_available_alleles))
            )
        return patients_available_alleles
