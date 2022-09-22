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

import os
from neofox.helpers import intermediate_files
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation, PredictedEpitope, Neoantigen
from neofox.model.factories import AnnotationFactory


class NeoagCalculator(object):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def _apply_gbm(self, tmp_in):
        """this function calls NeoAg tool. this tool applys a gradient boosting machine based on biochemical features
        to epitopes (predicted seqs)
        """
        my_path = os.path.abspath(os.path.dirname(__file__))
        model_path = os.path.join(my_path, "neoag-master")
        tool_path = os.path.join(
            my_path, "neoag-master/NeoAg_immunogenicity_predicition_GBM.R"
        )
        cmd = [self.configuration.rscript, tool_path, model_path, tmp_in]
        output, _ = self.runner.run_command(cmd)
        return output

    def _prepare_tmp_for_neoag(
        self,
        sample_id,
        mut_peptide,
        score_mut,
        ref_peptide,
        peptide_variant_position,
        tmp_file_name,
    ):
        """writes necessary epitope information into temporary file for neoag tool; only for epitopes with
        affinity < 500 nM
        """
        header = ["Sample_ID", "mut_peptide", "Reference", "peptide_variant_position"]
        try:
            epi_row = "\t".join(
                [sample_id, mut_peptide, ref_peptide, str(peptide_variant_position)]
            )
        except ValueError:
            epi_row = "\t".join(["NA", "NA", "NA", "NA"])
        with open(tmp_file_name, "w") as f:
            f.write("\t".join(header) + "\n")
            f.write(epi_row + "\n")

    def calculate_neoag_score(self, epitope: PredictedEpitope):
        tmp_file_name = intermediate_files.create_temp_file(
            prefix="tmp_neoag_", suffix=".txt"
        )
        self._prepare_tmp_for_neoag(
            "*****",
            epitope.mutated_peptide,
            epitope.affinity_mutated,
            epitope.wild_type_peptide,
            EpitopeHelper.position_of_mutation_epitope(epitope=epitope),
            tmp_file_name,
        )
        neoag_score = self._apply_gbm(tmp_file_name)
        os.remove(tmp_file_name)
        return neoag_score

    def get_annotation(self, epitope_mhci: PredictedEpitope, neoantigen: Neoantigen) -> Annotation:
        """wrapper function to determine neoag immunogenicity score for a mutated peptide sequence"""

        neoag_score = None
        if neoantigen.wild_type_xmer and epitope_mhci.mutated_peptide and epitope_mhci.wild_type_peptide:
            neoag_score = self.calculate_neoag_score(epitope=epitope_mhci)

        annotation = AnnotationFactory.build_annotation(value=neoag_score, name="Neoag_immunogenicity")
        return annotation
