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
import os
from neofox.exceptions import NeofoxCommandException
from neofox.helpers import intermediate_files
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Mhc1, PredictedEpitope, Zygosity, Neoantigen
from neofox.references.references import DependenciesConfiguration


PEPTIDE_LENGTHS = ["8", "9", "10", "11", "12", "13", "14"]


class NetMhcPanPredictor:

    def __init__(
            self, runner: Runner, configuration: DependenciesConfiguration,
            blastp_runner: BlastpRunner, mhc_parser: MhcParser):

        self.runner = runner
        self.configuration = configuration
        self.mhc_parser = mhc_parser
        self.blastp_runner = blastp_runner

    def mhc_prediction(self, available_alleles, sequence) -> List[PredictedEpitope]:
        """Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file."""

        input_file = intermediate_files.create_temp_fasta(sequences=[sequence], prefix="tmp_singleseq_")

        if available_alleles is None or available_alleles == "":
            raise NeofoxCommandException("None of the provided MHC I alleles are supported: {}".format(available_alleles))
        cmd = [
            self.configuration.net_mhc_pan,
            "-a",
            available_alleles,
            "-f",
            input_file,
            "-BA",
            "-l {}".format(",".join(PEPTIDE_LENGTHS))
        ]

        lines, _ = self.runner.run_command(cmd)
        os.remove(input_file)
        return self._parse_netmhcpan_output(lines)

    def mhc_prediction_peptide(self, alleles, sequence) -> PredictedEpitope:
        """
        Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file.
        peptide mode cannot use FASTA format and does not provide peptide lengths
        """

        result = None
        input_file = intermediate_files.create_temp_peptide(sequences=[sequence], prefix="tmp_singleseq_")

        if alleles is None or alleles == "":
            raise NeofoxCommandException("None of the provided MHC I alleles are supported: {}".format(alleles))
        cmd = [
            self.configuration.net_mhc_pan,
            "-a",
            alleles,
            "-p",
            input_file,
            "-BA"
        ]

        lines, _ = self.runner.run_command(cmd)
        predicted_epitopes = self._parse_netmhcpan_output(lines)
        if predicted_epitopes:
            result = predicted_epitopes[0]
        os.remove(input_file)
        return result

    def _parse_netmhcpan_output(self, lines: str) -> List[PredictedEpitope]:
        results = []
        for line in lines.splitlines():
            line = line.rstrip().lstrip()
            if line:
                if line.startswith(("#", "-", "HLA", "Prot", "Pos", "No")):
                    continue
                if "Distance to training dat" in line:
                    continue
                if line.startswith("Error"):
                    raise NeofoxCommandException("netmhcpan threw an error: {}".format(line))
                line = line.split()
                line = line[0:-2] if len(line) > 16 else line
                results.append(
                    PredictedEpitope(
                        position=int(line[0]),
                        allele_mhc_i=self.mhc_parser.parse_mhc_allele(line[1]),
                        mutated_peptide=line[2],
                        affinity_mutated=float(line[15]),
                        rank_mutated=float(line[12]),
                    )
                )
        return results

    def get_alleles_netmhcpan_representation(self, mhc: List[Mhc1]) -> List[str]:
        return list(
            map(
                self.mhc_parser.get_netmhcpan_representation, [a for m in mhc for a in m.alleles],
            )
        )

    def get_only_available_alleles(self, mhc_alleles: List[Mhc1], set_available_mhc: Set[str]) -> str:
        hla_alleles_names = self.get_alleles_netmhcpan_representation(mhc_alleles)
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

    def set_wt_netmhcpan_scores(self, predictions) -> List[PredictedEpitope]:
        for p in predictions:
            if p.wild_type_peptide is not None:
                wt_prediction = self.mhc_prediction_peptide(
                    alleles=self.mhc_parser.get_netmhcpan_representation(p.allele_mhc_i), sequence=p.wild_type_peptide)
                if wt_prediction is not None:
                    # NOTE: netmhcpan in peptide mode should return only one epitope
                    p.rank_wild_type = wt_prediction.rank_mutated
                    p.affinity_wild_type = wt_prediction.affinity_mutated
        return predictions

    def get_predictions(self, available_alleles, neoantigen: Neoantigen, uniprot) -> List[PredictedEpitope]:
        predictions = self.mhc_prediction(available_alleles, neoantigen.mutated_xmer)
        if neoantigen.wild_type_xmer:
            # make sure that predicted epitopes cover mutation in case of SNVs
            predictions = EpitopeHelper.filter_peptides_covering_snv(
                position_of_mutation=neoantigen.position, predictions=predictions
            )
        # make sure that predicted neoepitopes are not part of the WT proteome
        filtered_predictions = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions, uniprot=uniprot
        )
        return filtered_predictions

    def get_wt_predictions(self, available_alleles, neoantigen) -> List[PredictedEpitope]:
        predictions = self.mhc_prediction(available_alleles, neoantigen.wild_type_xmer)
        # make sure that predicted epitopes cover mutation in case of SNVs
        predictions = EpitopeHelper.filter_peptides_covering_snv(
            position_of_mutation=neoantigen.position, predictions=predictions
        )
        return predictions
