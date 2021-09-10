from datetime import datetime
from unittest import TestCase
import dotenv
import pkg_resources
import neofox.tests
from neofox.references.installer import NeofoxReferenceInstaller


class TestInstaller(TestCase):

    def setUp(self):
        dotenv.load_dotenv(override=True)
        pass

    def test_installer(self):
        reference_folder = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/test_installation_{:%Y%m%d%H%M%S}".format(datetime.now()))
        NeofoxReferenceInstaller(
            reference_folder=reference_folder, install_r_dependencies=False
        ).install()
