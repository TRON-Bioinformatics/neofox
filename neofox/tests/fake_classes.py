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
import pkg_resources
import neofox.tests
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import (
    PredictedEpitope,
)
from neofox.references.references import (
    ReferenceFolder,
    AvailableAlleles,
    DependenciesConfiguration, HlaDatabase,
)


class FakeReferenceFolder(ReferenceFolder):

    available_mhc_i_alleles = ["A", "B"]
    available_mhc_ii_alleles = ["A", "B"]

    def _check_reference_genome_folder(self):
        return os.environ.get(neofox.REFERENCE_FOLDER_ENV, "")

    def _check_resources(self):
        pass

    def get_available_alleles(self):
        return FakeAvailableAlleles(
            available_mch_i=self.available_mhc_i_alleles,
            available_mch_ii=self.available_mhc_ii_alleles,
        )


class FakeDependenciesConfiguration(DependenciesConfiguration):
    def _check_and_load_binary(self, variable_name, optional=False):
        return os.environ.get(variable_name, "some_non_empty_fake_value")


class FakeBestAndMultipleBinder(BestAndMultipleBinder):
    def __init__(self, mutated_epitope, affinity, wild_type_epitope):
        self.best_epitope_by_affinity = PredictedEpitope(
            peptide=mutated_epitope, affinity_score=affinity, pos=1, hla="sdf", rank=1
        )
        self.best_wt_epitope_by_affinity = PredictedEpitope(
            peptide=wild_type_epitope, affinity_score=30, pos=1, hla="sdf", rank=1
        )


class FakeAvailableAlleles(AvailableAlleles):

    def __init__(self, available_mch_i=[], available_mch_ii=[]):
        self.available_mhc_i = available_mch_i
        self.available_mhc_ii = available_mch_ii

    def get_available_mhc_i(self):
        return self.available_mhc_i

    def get_available_mhc_ii(self):
        return self.available_mhc_ii


class FakeHlaDatabase(HlaDatabase):

    def __init__(self):
        super().__init__(hla_database_filename=pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/hla_database.txt"))
