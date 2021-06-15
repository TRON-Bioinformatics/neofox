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

from neofox.model.wrappers import AnnotationFactory, NOT_AVAILABLE_VALUE


class TestAnnotationFactory(TestCase):
    def test_float_values_precision_limit(self):
        self.assertEqual(
            "0.12346", AnnotationFactory.build_annotation("test", 0.123456789).value
        )
        self.assertEqual(
            "0.12346", AnnotationFactory.build_annotation("test", 0.1234567891234).value
        )
        self.assertEqual(
            "0.12345", AnnotationFactory.build_annotation("test", 0.12345).value
        )
        self.assertEqual(
            "0.123", AnnotationFactory.build_annotation("test", 0.123).value
        )

    def test_none(self):
        self.assertEqual(NOT_AVAILABLE_VALUE, AnnotationFactory.build_annotation("test", None).value)

    def test_booleans(self):
        self.assertEqual("1", AnnotationFactory.build_annotation("test", True).value)
        self.assertEqual("0", AnnotationFactory.build_annotation("test", False).value)

    def test_integers(self):
        self.assertEqual("123", AnnotationFactory.build_annotation("test", 123).value)
        self.assertEqual("-123", AnnotationFactory.build_annotation("test", -123).value)

    def test_strings(self):
        # TODO: I am wondering if we should transform this to NA. Franziska what do you think?
        self.assertEqual("", AnnotationFactory.build_annotation("test", "").value)
        self.assertEqual(
            "blabla", AnnotationFactory.build_annotation("test", "blabla").value
        )
