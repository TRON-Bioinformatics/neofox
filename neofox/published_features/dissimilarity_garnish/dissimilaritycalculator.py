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
import os.path
from typing import List

from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)


class DissimilarityCalculator(BlastpRunner):
    def __init__(self, runner, configuration, proteome_db):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        super().__init__(runner=runner, configuration=configuration)
        self.proteome_db = proteome_db

    def calculate_dissimilarity(self, mhc_mutation, mhc_affinity, filter_binder=False):
        """
        wrapper for dissimilarity calculation
        """
        dissimilarity = None
        if mhc_mutation != "-" and (not filter_binder or not mhc_affinity >= 500):
            similarity = self.run_blastp(
                peptide=mhc_mutation,
                database=os.path.join(self.proteome_db, "homo_sapiens"),
                a=32,
            )
            if similarity is not None:
                dissimilarity = 1 - similarity
        return dissimilarity

    def get_annotations(self, netmhcpan: BestAndMultipleBinder) -> List[Annotation]:
        """
        returns dissimilarity for MHC I (affinity) MHC II (affinity)
        """
        annotations = []
        if netmhcpan.best_epitope_by_affinity:
            dissimilarity = self.calculate_dissimilarity(mhc_mutation=netmhcpan.best_epitope_by_affinity.peptide,
                                                         mhc_affinity=netmhcpan.best_epitope_by_affinity.affinity_score,
                                                         filter_binder=True, )
            if dissimilarity is not None:
                annotations = [
                    AnnotationFactory.build_annotation(
                        value=dissimilarity,
                        name="Dissimilarity_MHCI_cutoff500nM",
                    ),
                ]
        return annotations
