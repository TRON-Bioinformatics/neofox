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

from unittest import TestCase

from neofox.model.wrappers import NOT_AVAILABLE_VALUE
from neofox.published_features.expression import Expression


class TestExpression(TestCase):
    def test_calculate_expression_mutation(self):
        result = Expression(transcript_expression=12.0, vaf_rna=0.2).get_annotations()[
            0
        ]
        self.assertGreater(result.value, "0.0")
        # no reads for mut
        result = Expression(transcript_expression=12.0, vaf_rna=0.0).get_annotations()[
            0
        ]
        self.assertEqual(result.value, "0")
        # no reads for mut/wt
        result = Expression(transcript_expression=12.0, vaf_rna=-1).get_annotations()[0]
        self.assertEqual(result.value, NOT_AVAILABLE_VALUE)
        result = Expression(transcript_expression=None, vaf_rna=-1).get_annotations()[0]
        self.assertEqual(result.value, NOT_AVAILABLE_VALUE)
