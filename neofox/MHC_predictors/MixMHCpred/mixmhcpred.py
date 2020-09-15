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
from logzero import logger
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.MixMHCpred.abstract_mixmhcpred import AbstractMixMHCpred
from neofox.helpers import intermediate_files


class MixMHCpred(AbstractMixMHCpred):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self._initialise()

    def _initialise(self):
        self.best_peptide = None
        self.best_score = None
        self.best_rank = None
        self.best_allele = None


    def _mixmhcprediction(self, hla_alleles, tmpfasta, outtmp):
        """ Performs MixMHCpred prediction for desired hla allele and writes result to temporary file.
        """
        allels_for_prediction = []
        for allele in hla_alleles:
            allele = allele.replace("*", "")
            allele = allele.replace("HLA-", "")
            allels_for_prediction.append(allele)
        hla_allele = ",".join(allels_for_prediction)
        self.runner.run_command(cmd=[
            self.configuration.mix_mhc_pred,
            "-a", hla_allele,
            "-i", tmpfasta,
            "-o", outtmp])

    def _extract_best_peptide_per_mutation(self, ligand_data_tuple):
        """extract best predicted ligand per mutation
        """
        head = ligand_data_tuple[0]
        predicted_ligands = ligand_data_tuple[1]
        index_peptide = head.index("Peptide")
        index_score = head.index("Score_bestAllele")
        index_allele = head.index("BestAllele")
        index_rank = head.index("%Rank_bestAllele")
        min_value = -1000000000000000000
        for ii, i in enumerate(predicted_ligands):
            ligand_information = [str(i[index_peptide]), str(i[index_score]), str(i[index_rank]), str(i[index_allele])]
            # best ligand per mutation
            if float(i[index_score]) > float(min_value):
                min_value = i[index_score]
                best_ligand = ligand_information
        head_new = ["Peptide", "Score_bestAllele", "%Rank_bestAllele", "BestAllele"]
        return head_new, best_ligand

    def run(self, sequence_wt, sequence_mut, alleles):
        """Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        """
        self._initialise()
        tmp_prediction = intermediate_files.create_temp_file(prefix="mixmhcpred", suffix=".txt")
        seqs = self.generate_nmers(xmer_wt=sequence_wt, xmer_mut=sequence_mut, lengths=[8, 9, 10, 11])
        tmp_fasta = intermediate_files.create_temp_fasta(seqs, prefix="tmp_sequence_")
        self._mixmhcprediction(alleles, tmp_fasta, tmp_prediction)
        all_predicted_ligands = self.read_mixmhcpred(tmp_prediction)
        if len(all_predicted_ligands[1]) > 0:
            best_predicted_ligand = self._extract_best_peptide_per_mutation(all_predicted_ligands)
            self.best_peptide = self.add_best_epitope_info(best_predicted_ligand, "Peptide")
            self.best_score = float(self.add_best_epitope_info(best_predicted_ligand, "Score_bestAllele"))
            self.best_rank = self.add_best_epitope_info(best_predicted_ligand, "%Rank_bestAllele")
            self.best_allele = self.add_best_epitope_info(best_predicted_ligand, "BestAllele")


    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(value=self.best_peptide, name="MixMHCpred_best_peptide"),
            AnnotationFactory.build_annotation(value=self.best_score, name="MixMHCpred_best_score"),
            AnnotationFactory.build_annotation(value=self.best_rank, name="MixMHCpred_best_rank"),
            AnnotationFactory.build_annotation(value=self.best_allele, name="MixMHCpred_best_allele"),
            ]
