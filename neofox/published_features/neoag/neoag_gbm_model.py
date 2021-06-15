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
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)
from neofox import AFFINITY_THRESHOLD_DEFAULT


class NeoagCalculator(object):

    def __init__(self, runner, configuration, affinity_threshold=AFFINITY_THRESHOLD_DEFAULT):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.affinity_threshold = affinity_threshold
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
            if score_mut < self.affinity_threshold:
                epi_row = "\t".join(
                    [sample_id, mut_peptide, ref_peptide, str(peptide_variant_position)]
                )
            else:
                epi_row = "\t".join(["NA", "NA", "NA", "NA"])
        except ValueError:
            epi_row = "\t".join(["NA", "NA", "NA", "NA"])
        with open(tmp_file_name, "w") as f:
            f.write("\t".join(header) + "\n")
            f.write(epi_row + "\n")

    def get_annotation(
        self, sample_id, netmhcpan: BestAndMultipleBinder, peptide_variant_position, mutation
    ) -> Annotation:
        """wrapper function to determine neoag immunogenicity score for a mutated peptide sequence"""

        neoag_score = None
        if mutation.wild_type_xmer and netmhcpan.best_epitope_by_affinity.peptide and netmhcpan.best_wt_epitope_by_affinity.peptide:
            # TODO: move this tmp file creation inside the method
            tmp_file_name = intermediate_files.create_temp_file(
                prefix="tmp_neoag_", suffix=".txt"
            )
            self._prepare_tmp_for_neoag(
                sample_id,
                netmhcpan.best_epitope_by_affinity.peptide,
                netmhcpan.best_epitope_by_affinity.affinity_score,
                netmhcpan.best_wt_epitope_by_affinity.peptide,
                peptide_variant_position,
                tmp_file_name,
            )
            neoag_score = self._apply_gbm(tmp_file_name)
        annotation = AnnotationFactory.build_annotation(
            value=neoag_score, name="Neoag_immunogenicity"
        )
        return annotation
