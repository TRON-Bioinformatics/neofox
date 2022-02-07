import os
from typing import List
from unittest import TestCase
import pkg_resources
import neofox
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import Neoantigen, Patient
from neofox.tests.integration_tests import integration_test_tools
from neofox.tests.synthetic_data.data_generator import DataGenerator


class TestSyntheticData(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.data_generator = DataGenerator(self.references, self.configuration)

    def _write_data(self, neoantigens: List[Neoantigen], neoantigens_filename,
                   patients: List[Patient], patients_filename):
        ModelConverter._neoantigens2table(neoantigens).to_csv(neoantigens_filename, sep="\t", index=False)
        ModelConverter.patients2table(patients).to_csv(patients_filename, sep="\t", index=False)

    def test_1patients_10neoantigens(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=1, num_neoantigens_per_patient=10, wildtype=True)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_1patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_1patients_10neoantigens.{}.txt".format(i)))

    def test_10patients_10neoantigens(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=10, num_neoantigens_per_patient=10, wildtype=True)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_10patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_10patients_10neoantigens.{}.txt".format(i)))

    def test_100patients_10neoantigens(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=100, num_neoantigens_per_patient=10, wildtype=True)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_100patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_100patients_10neoantigens.{}.txt".format(i)))

    def test_1000patients_10neoantigens(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=1000, num_neoantigens_per_patient=10, wildtype=True)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_1000patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_1000patients_10neoantigens.{}.txt".format(i)))

    def test_1patients_10neoantigens_no_wt(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=1, num_neoantigens_per_patient=10, wildtype=False)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_no_wt_1patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_no_wt_1patients_10neoantigens.{}.txt".format(i)))

    def test_10patients_10neoantigens_no_wt(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=10, num_neoantigens_per_patient=10, wildtype=False)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_no_wt_10patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_no_wt_10patients_10neoantigens.{}.txt".format(i)))

    def test_100patients_10neoantigens_no_wt(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=100, num_neoantigens_per_patient=10, wildtype=False)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_no_wt_100patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_no_wt_100patients_10neoantigens.{}.txt".format(i)))

    def test_1000patients_10neoantigens_no_wt(self):
        for i in range(5):
            patients, neoantigens = self.data_generator.generate_data(
                num_patients=1000, num_neoantigens_per_patient=10, wildtype=False)
            self._write_data(
                neoantigens=neoantigens,
                neoantigens_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/neoantigens_no_wt_1000patients_10neoantigens.{}.txt".format(i)),
                 patients=patients,
                patients_filename=pkg_resources.resource_filename(
                    neofox.tests.__name__,
                    "resources/synthetic_data/patients_no_wt_1000patients_10neoantigens.{}.txt".format(i)))
