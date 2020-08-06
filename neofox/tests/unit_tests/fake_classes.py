import os

import neofox
from neofox.references.references import ReferenceFolder


class FakeReferenceFolder(ReferenceFolder):

    @staticmethod
    def _check_reference_genome_folder():
        return os.environ.get(neofox.REFERENCE_FOLDER_ENV, "")

    @staticmethod
    def _check_resources(resources):
        pass
