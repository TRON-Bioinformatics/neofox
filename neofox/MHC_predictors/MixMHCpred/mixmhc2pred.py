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

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.MixMHCpred.abstract_mixmhcpred import AbstractMixMHCpred
from neofox.helpers import intermediate_files


class MixMhc2Pred(AbstractMixMHCpred):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.available_alleles = self.load_available_alleles()
        self._initialise()

    def _initialise(self):
        self.best_peptide = None
        self.best_rank = None
        self.best_allele = None


    def load_available_alleles(self):
        """
        loads file with available HLA II alllels for MixMHC2pred prediction, returns set
        :return:
        """
        path_to_hlaii_file = self.configuration.mix_mhc2_pred_alleles_list
        available_alleles = []
        with open(path_to_hlaii_file) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("L", "A")):
                        continue
                    line1 = line.split()[0]
                    if line1 is not None:
                        available_alleles.append(line1)
        return available_alleles

    def combine_dq_dp_alleles(self, list_alleles):
        """ returns patient HLA-DQ/HLA-DP allele combination that are relevant for MixMHC2pred
        """
        alleles_pairs = ["__".join([allele_1, allele_2]) for allele_1 in list_alleles for allele_2 in list_alleles
                              if allele_1 != allele_2]
        alleles_triplets = ["__".join([allele_1, allele_2, allele_3]) for allele_1 in list_alleles
                                 for allele_2 in list_alleles for allele_3 in list_alleles
                                 if allele_1 != allele_2 and allele_1 != allele_3 and allele_2 != allele_3]
        all_possible_alleles = alleles_pairs + alleles_triplets
        alleles_for_prediction = [allele for allele in all_possible_alleles if allele in self.available_alleles]
        return alleles_for_prediction

    def transform_hla_ii_alleles_for_prediction(self, hla_ii_alleles):
        """ prepares list of HLA II alleles for prediction in required format
        """
        alleles_for_prediction = []
        alleles_dq = []
        alleles_dp = []
        for allele in hla_ii_alleles:
            allele = allele.replace("*", "_").replace(":", "_").replace("HLA-", "")
            if allele.startswith("DR"):
                if allele in self.available_alleles:
                    alleles_for_prediction.append(allele)
            elif allele.startswith("DP"):
                alleles_dp.append(allele)
            elif allele.startswith("DQ"):
                alleles_dq.append(allele)
        dp_alleles_for_prediction = self.combine_dq_dp_alleles(alleles_dp)
        dq_alleles_for_prediction = self.combine_dq_dp_alleles(alleles_dq)
        alleles_for_prediction = alleles_for_prediction + dq_alleles_for_prediction + dp_alleles_for_prediction
        hla_ii_alleles_transformed = " ".join(alleles_for_prediction)
        return hla_ii_alleles_transformed

    def mixmhc2prediction(self, hla_ii_alleles, tmpfasta, outtmp):
        """ Performs MixMHC2pred prediction for desired hla allele and writes result to temporary file.
        """
        hla_ii_alleles_for_prediction = self.transform_hla_ii_alleles_for_prediction(hla_ii_alleles)
        cmd = [
            self.configuration.mix_mhc2_pred,
            "-a", hla_ii_alleles_for_prediction,
            "-i", tmpfasta,
            "-o", outtmp]
        self.runner.run_command(cmd)

    def extract_best_peptide_per_mutation(self, ligand_data_tuple):
        """extract best predicted ligand per mutation
        """
        head = ligand_data_tuple[0]
        predicted_ligands = ligand_data_tuple[1]
        index_peptide = head.index("Peptide")
        index_allele = head.index("BestAllele")
        index_rank = head.index("%Rank")
        min_value = 1000000000000000000
        for ii, i in enumerate(predicted_ligands):
            ligand_information = [str(i[index_peptide]), str(i[index_rank]), str(i[index_allele])]
            # best ligand per mutation
            if float(i[index_rank]) < float(min_value):
                min_value = i[index_rank]
                best_ligand = ligand_information
        head_new = ["Peptide", "%Rank", "BestAllele"]
        return head_new, best_ligand

    def run(self, alleles, sequence_wt, sequence_mut):
        """Runs MixMHC2pred:
        prediction for peptides of length 13 to 18 based on Suppl Fig. 6 a in Racle, J., et al., Nat. Biotech. (2019).
        Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes.
        """
        self._initialise()
        tmp_prediction = intermediate_files.create_temp_file(prefix="mixmhc2pred", suffix=".txt")
        potential_ligand_sequences = self.generate_nmers(xmer_wt=sequence_wt, xmer_mut=sequence_mut,
                                                         lengths=[13, 14, 15, 16, 17, 18])
        tmp_fasta = intermediate_files.create_temp_fasta(potential_ligand_sequences, prefix="tmp_sequence_")
        # try except statement to prevent stop of neofox for mps shorter < 13aa
        # TODO: this needs to be fixed, we could filter the list of nmers by length
        if len(potential_ligand_sequences) > 0:
            try:
                self.mixmhc2prediction(alleles, tmp_fasta, tmp_prediction)
            except:
                pass
            # TODO: also all of this try-catch needs to be fixed, in general the risk here is that they hide errors
            try:
                all_predicted_ligands = self.read_mixmhcpred(tmp_prediction)
            except:
                pass
            best_predicted_ligand = self.extract_best_peptide_per_mutation(all_predicted_ligands)
            self.best_peptide = self.add_best_epitope_info(best_predicted_ligand, "Peptide")
            # TODO: improve how data is fetched so types are maintained
            self.best_rank = float(self.add_best_epitope_info(best_predicted_ligand, "%Rank"))
            self.best_allele = self.add_best_epitope_info(best_predicted_ligand, "BestAllele")

    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(value=self.best_peptide, name="MixMHC2pred_best_peptide"),
            AnnotationFactory.build_annotation(value=self.best_rank, name="MixMHC2pred_best_rank"),
            AnnotationFactory.build_annotation(value=self.best_allele, name="MixMHC2pred_best_allele"),
        ]
