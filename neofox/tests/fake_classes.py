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
