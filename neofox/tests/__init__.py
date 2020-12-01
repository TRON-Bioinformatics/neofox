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
from neofox.model.conversion import ModelConverter

TEST_MHC_ONE = ModelConverter.parse_mhc1_alleles(
    [
        "HLA-A*24:02",
        "HLA-A*02:01",
        "HLA-B*15:01",
        "HLA-B*44:02",
        "HLA-C*07:02",
        "HLA-C*05:01",
    ]
)

TEST_MHC_TWO = ModelConverter.parse_mhc2_alleles(
    [
        "HLA-DRB1*04:02",
        "HLA-DRB1*08:01",
        "HLA-DQA1*03:01",
        "HLA-DQA1*04:01",
        "HLA-DQB1*03:02",
        "HLA-DQB1*04:02",
        "HLA-DPA1*01:03",
        "HLA-DPA1*02:01",
        "HLA-DPB1*13:01",
        "HLA-DPB1*04:01",
    ]
)
