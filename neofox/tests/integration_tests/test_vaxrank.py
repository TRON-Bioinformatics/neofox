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

from neofox.helpers.runner import Runner
from neofox.model.neoantigen import PredictedEpitope
from neofox.published_features.vaxrank.vaxrank import (
    VaxRank,
)


class TestVaxRank(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        self.vaxrank = VaxRank()
        self.hla_database = self.references.get_mhc_database()

    def test_vaxrank(self):
        """
        Test whether the VaxRank calculation works correct. 

        The binding score of one list of neoepitopes is calculated and compared
        with the true value. As floating point calculations can differ between 
        architectures and the vaxrank annotations are rounded, we only compare on 
        almost equal.
        Additionaly it is checked if the total score and the total imputed score
        is correctly calculated.
        """
        EXPRESSION_SCORE = 10.4
        IMPUTED_SCORE = 7.6

        neoepitopes = integration_test_tools.get_test_epitopes(self.hla_database)

        annot = VaxRank().get_annotations(
            epitope_predictions = neoepitopes,
            expression_score = EXPRESSION_SCORE, 
            imputed_score = IMPUTED_SCORE
            )

        # to access the values of the Vaxrank annotation even if the order changes
        annot_dict = {a.name: float(a.value) for a in annot}

        # the binding score calculated by hand
        assert abs(annot_dict["Vaxrank_bindingScore"] - 2.108220805) < 0.001
        assert abs(annot_dict["Vaxrank_bindingScore"] * EXPRESSION_SCORE - annot_dict["Vaxrank_totalScore"]) < 0.001
        assert abs(annot_dict["Vaxrank_bindingScore"] * IMPUTED_SCORE - annot_dict["Vaxrank_totalScore_imputed"]) < 0.001
