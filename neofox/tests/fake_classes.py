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

import neofox
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.references.references import ReferenceFolder


class FakeReferenceFolder(ReferenceFolder):

    @staticmethod
    def _check_reference_genome_folder():
        return os.environ.get(neofox.REFERENCE_FOLDER_ENV, "")

    @staticmethod
    def _check_resources(resources):
        pass


class FakeBestAndMultipleBinder(BestAndMultipleBinder):

    def __init__(self, mutated_epitope, affinity, wild_type_epitope):
        self.best4_affinity_epitope = mutated_epitope
        self.best4_affinity = affinity
        self.best4_affinity_epitope_WT = wild_type_epitope
