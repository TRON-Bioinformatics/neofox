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
from neofox.model.neoantigen import Annotation, Mhc2, Mhc2GeneName, MhcAllele, Mhc2Name
import re

NOT_AVAILABLE_VALUE = "NA"


class AnnotationFactory(object):
    @staticmethod
    def build_annotation(name, value):
        if isinstance(value, bool):
            # prints booleans as 0/1 strings
            value = "1" if value else "0"
        if not isinstance(value, str) and not isinstance(value, type(None)):
            if isinstance(value, float):
                value = "{0:.5g}".format(round(value, 5))
            else:
                value = str(value)
        if value is None:
            value = NOT_AVAILABLE_VALUE
        return Annotation(name=name, value=value)


def get_alleles_by_gene(
    mhc_isoforms: List[Mhc2], gene: Mhc2GeneName
) -> List[MhcAllele]:
    return [
        a for m in mhc_isoforms for g in m.genes if g.name == gene for a in g.alleles
    ]


def get_mhc2_isoform_name(a: MhcAllele, b: MhcAllele):
    # NOTE: this is needed as jus setting alpha chain to None wouldn't work with protobuf
    if a is not None and a.name:
        return "{}-{}".format(a.name, b.name.replace("HLA-", ""))
    else:
        return b.name
