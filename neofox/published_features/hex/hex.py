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
from neofox.model.neoantigen import Annotation, PredictedEpitope
from neofox.model.factories import AnnotationFactory
from neofox.published_features.hex.pyhex import PyHex
from neofox.references.references import ReferenceFolder


class Hex(object):

    def __init__(self, references: ReferenceFolder):
        self.iedb_fasta = references.get_iedb_fasta()
        self.pyhex = PyHex(self.iedb_fasta)

    def apply_hex(self, mut_peptide):
        """this function calls hex tool. this tool analyses the neoepitope candidate sequence for molecular mimicry to viral epitopes
        """
        return self.pyhex.run(mut_peptide)

    def get_annotation(
            self, mutated_peptide_mhci: PredictedEpitope, mutated_peptide_mhcii: PredictedEpitope) -> List[Annotation]:
        """
        wrapper function for HEX (Homology evaluation of Xenopeptides) (Chiaro et al., 2021)
        """
        # TODO: add annotation of b score when re-implemented in python. The annotation is too slow for bigger datasets
        hex_aln_score_mhci = None
        hex_aln_score_mhcii = None
        # hex_b_score = None
        if mutated_peptide_mhci and mutated_peptide_mhci.mutated_peptide:
            # hex_aln_score, hex_b_score = self.apply_hex(netmhcpan.best_epitope_by_affinity.peptide).split(" ")
            hex_aln_score_mhci = self.apply_hex(mutated_peptide_mhci.mutated_peptide)
        if mutated_peptide_mhcii and mutated_peptide_mhcii.mutated_peptide:
            hex_aln_score_mhcii = self.apply_hex(mutated_peptide_mhcii.mutated_peptide)
        annotations = [
            AnnotationFactory.build_annotation(
                value=hex_aln_score_mhci, name="HexAlignmentScore_MHCI"),
            AnnotationFactory.build_annotation(
                value=hex_aln_score_mhcii, name="HexAlignmentScore_MHCII")
            # AnnotationFactory.build_annotation(
             #   value=hex_b_score, name="hex_B_score"
            #)
        ]
        return annotations

    def get_annotations_epitope(self, epitope: PredictedEpitope) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.apply_hex(mut_peptide=epitope.mutated_peptide),
                name='hex_alignment_score')
        ]
