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
from typing import List

from neofox.helpers.epitope_helper import EpitopeHelper

from neofox.model.neoantigen import Annotation, Mhc1, MhcAllele
from neofox.model.wrappers import AnnotationFactory
from neofox.helpers import intermediate_files
import pandas as pd
import os
from logzero import logger

ALLELE = "BestAllele"
RANK = "%Rank_bestAllele"
PEPTIDE = "Peptide"
SCORE = "Score_bestAllele"


class MixMHCpred:

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    @staticmethod
    def _get_mixmhc_allele_representation(mhc_alleles: List[MhcAllele]):
        return list(map(
            lambda x: "{gene}{group}{protein}".format(gene=x.gene, group=x.group, protein=x.protein), mhc_alleles))

    def _mixmhcprediction(self, mhc_isoforms: List[Mhc1], potential_ligand_sequences) -> pd.DataFrame:
        """
        Performs MixMHCpred prediction for desired hla allele and writes result to temporary file.
        """
        outtmp = intermediate_files.create_temp_file(prefix="mixmhcpred", suffix=".txt")
        tmpfasta = intermediate_files.create_temp_fasta(potential_ligand_sequences, prefix="tmp_sequence_")
        self.runner.run_command(cmd=[
            self.configuration.mix_mhc_pred,
            "-a", ",".join(self._get_mixmhc_allele_representation([a for m in mhc_isoforms for a in m.alleles])),
            "-i", tmpfasta,
            "-o", outtmp])
        results = pd.read_csv(outtmp, sep="\t", comment="#")
        os.remove(outtmp)
        return results

    def run(self, sequence_wt, sequence_mut, mhc: List[Mhc1]):
        """Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        """
        best_peptide = None
        best_rank = None
        best_allele = None
        best_score = None
        potential_ligand_sequences = EpitopeHelper.generate_nmers(
            xmer_wt=sequence_wt, xmer_mut=sequence_mut, lengths=[8, 9, 10, 11])
        if len(potential_ligand_sequences) > 0:
            results = self._mixmhcprediction(mhc, potential_ligand_sequences)
            # get best result by maximum score
            best_result = results[results[SCORE] == results[SCORE].max()]
            try:
                best_peptide = best_result[PEPTIDE].iat[0]
                best_rank = best_result[RANK].iat[0]
                best_allele = best_result[ALLELE].iat[0]
                best_score = best_result[SCORE].iat[0]
            except IndexError:
                logger.info("MixMHCpred returned no best result")
        return best_peptide, best_rank, best_allele, best_score

    def get_annotations(self, sequence_wt, sequence_mut, mhc: List[Mhc1]) -> List[Annotation]:
        best_peptide, best_rank, best_allele, best_score = self.run(
            mhc=mhc, sequence_wt=sequence_wt, sequence_mut=sequence_mut)
        return [
            AnnotationFactory.build_annotation(value=best_peptide, name="MixMHCpred_best_peptide"),
            AnnotationFactory.build_annotation(value=best_score, name="MixMHCpred_best_score"),
            AnnotationFactory.build_annotation(value=best_rank, name="MixMHCpred_best_rank"),
            AnnotationFactory.build_annotation(value=best_allele, name="MixMHCpred_best_allele"),
            ]
