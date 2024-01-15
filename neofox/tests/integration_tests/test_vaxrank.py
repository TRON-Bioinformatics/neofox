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
from unittest import TestCase

import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.model.factories import AnnotationFactory

from neofox.helpers.runner import Runner
from neofox.model.neoantigen import PredictedEpitope
from neofox.published_features.vaxrank.vaxrank import (
    VaxRank,
)
from neofox.tests.fake_classes import FakePredictedEpitope

class TestVaxRank(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        self.vaxrank = VaxRank()
        self.AnnotationFactory = AnnotationFactory()
        self.hla_database = self.references.get_mhc_database()

        # Generate epitopes for test case
        # NOTE: these should remain the same, to check if something in the 
        #       vaxrank calculation has changed
        epitopes = [
            "VIFSGEGSDEF", "VIFSGEGSDEF", "VIFSGEGSDEF", "IFSGEGSDEFT", "IFSGEGSDEFT", "IFSGEGSDEFT"
        ]
        alleles = [
            "HLA-A*24:02", "HLA-A*02:01", "HLA-B*15:01", "HLA-A*24:02", "HLA-A*02:01", "HLA-B*15:01"
        ]
        affinities = [1.74, 7.67, 802.78, 2057.1, 770.35, 5002.44]
        self.pred_epitopes = [FakePredictedEpitope(e, m, a) for e, m, a in zip(epitopes, alleles, affinities)]

        # expected binding value
        self.TRUE_BINDING_SCORE = 2.1082247484257532

    def test_vaxrank_binding(self):
        """
        Test if the VaxRank binding score is calculated correctly.

        The binding score of one list of neoepitopes is calculated and compared
        with the manually calculated value.
        """

        binding_score = self.vaxrank.total_binding(self.pred_epitopes)

        self.assertAlmostEqual(binding_score, self.TRUE_BINDING_SCORE)

    def test_vaxrank(self):
        """
        Test whether the VaxRank calculation works correct. 

        Check if the total score and the total imputed score are calculated correctly.
        """
        EXPRESSION_SCORE = 10.4
        IMPUTED_SCORE = 7.6

        binding_score = self.vaxrank.total_binding(self.pred_epitopes)

        # expected annotation
        expected_total_score = AnnotationFactory.build_annotation(
                value=EXPRESSION_SCORE * binding_score, 
                name="Vaxrank_totalScore"
            )
        expected_total_score_imputed = AnnotationFactory.build_annotation(
                value=IMPUTED_SCORE * binding_score, 
                name="Vaxrank_totalScore_imputed"
            )

        # calculated annotation
        annot = self.vaxrank.get_annotations(
            epitope_predictions = self.pred_epitopes,
            expression_score = EXPRESSION_SCORE, 
            imputed_score = IMPUTED_SCORE
            )

        # access the values of the Vaxrank annotation even if the order changes
        annot_dict = {a.name: a.value for a in annot}

        # since strings with the rounded value are returned by the get_annotations
        # function, we can compare here on equality
        assert annot_dict["Vaxrank_totalScore"] == expected_total_score.value
        assert annot_dict["Vaxrank_totalScore_imputed"] == expected_total_score_imputed.value
