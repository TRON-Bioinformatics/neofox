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
import os

from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import PredictedEpitope
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.references.references import ReferenceFolder, IEDB_FASTA


class Hex(object):

    def __init__(self, references: ReferenceFolder, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.iedb_fasta = os.path.join(references.iedb, IEDB_FASTA)

    def apply_hex(self, mut_peptide):
        """this function calls hex tool. this tool analyses the neoepitope candidate sequence for molecular mimicry to viral epitopes
        """
        my_path = os.path.abspath(os.path.dirname(__file__))
        tool_path = os.path.join(my_path, "hex.R")
        cmd = [self.configuration.rscript, tool_path, mut_peptide, self.iedb_fasta, my_path]
        output, _ = self.runner.run_command(cmd)
        if output == "":
            output = None
        return output

    def get_annotation(
            self, mutated_peptide_mhci: PredictedEpitope, mutated_peptide_mhcii: PredictedEpitope) -> List[Annotation]:
        """
        wrapper function for HEX (Homology evaluation of Xenopeptides) (Chiaro et al., 2021)
        """
        # TODO: add annotation of b score when re-implemented in python. The annotation is too slow for bigger datasets
        hex_aln_score_mhci = None
        hex_aln_score_mhcii = None
        # hex_b_score = None
        if mutated_peptide_mhci and mutated_peptide_mhci.peptide:
            # hex_aln_score, hex_b_score = self.apply_hex(netmhcpan.best_epitope_by_affinity.peptide).split(" ")
            hex_aln_score_mhci = self.apply_hex(mutated_peptide_mhci.peptide)
        if mutated_peptide_mhcii and mutated_peptide_mhcii.peptide:
            hex_aln_score_mhcii = self.apply_hex(mutated_peptide_mhcii.peptide)
        annotations = [
            AnnotationFactory.build_annotation(
                value=hex_aln_score_mhci, name="Hex_alignment_score_MHCI"),
            AnnotationFactory.build_annotation(
                value=hex_aln_score_mhcii, name="Hex_alignment_score_MHCII")
            # AnnotationFactory.build_annotation(
             #   value=hex_b_score, name="hex_B_score"
            #)
        ]
        return annotations
