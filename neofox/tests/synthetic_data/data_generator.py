import os
from typing import List, Tuple
from faker import Faker
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMHC2pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Patient, Neoantigen
from neofox.references.references import ReferenceFolder, HOMO_SAPIENS_FASTA, DependenciesConfiguration
from neofox.tests.synthetic_data.factories import PatientProvider, NeoantigenProvider


class DataGenerator:

    def __init__(self, reference_folder: ReferenceFolder, configuration: DependenciesConfiguration):

        self.hla_database = reference_folder.get_mhc_database()

        faker = Faker()
        mixmhcpred_alleles = set(self.load_mhc1_alleles(
            MixMHCpred(None, configuration=configuration, mhc_parser=None).available_alleles))
        netmhcpan_alleles = set(self.load_mhc1_alleles(
            reference_folder.get_available_alleles().get_available_mhc_i()))
        mhc1_alleles = mixmhcpred_alleles.union(netmhcpan_alleles)

        mixmhc2pred_alleles = set(self.load_mhc2_alleles(
            MixMHC2pred(runner=None, configuration=configuration, mhc_parser=None).available_alleles))
        netmhc2pan_alleles = set(self.load_mhc2_alleles(
            reference_folder.get_available_alleles().get_available_mhc_ii()))
        mhc2_isoforms = mixmhc2pred_alleles.union(netmhc2pan_alleles)

        self.patient_provider = PatientProvider(faker, mhc1_alleles, mhc2_isoforms, self.hla_database)
        self.neoantigen_provider = NeoantigenProvider(
            faker, proteome_fasta=os.path.join(reference_folder.proteome_db, HOMO_SAPIENS_FASTA))

    def load_mhc1_alleles(self, available_alleles: List[str]):
        mhc_alleles = []
        for a in available_alleles:
            try:
                parsed_allele = MhcParser.get_mhc_parser(self.hla_database).parse_mhc_allele(a)
            except AssertionError:
                continue
            mhc_alleles.append(parsed_allele.name)
        return mhc_alleles

    def load_mhc2_alleles(self, available_alleles: List[str]):
        mhc_alleles = []
        for a in available_alleles:
            try:
                parsed_allele = MhcParser.get_mhc_parser(self.hla_database).parse_mhc2_isoform(a)
            except AssertionError:
                continue
            mhc_alleles.append(parsed_allele.name)
        return mhc_alleles

    def generate_data(self, num_patients, num_neoantigens_per_patient, wildtype) -> Tuple[List[Patient], List[Neoantigen]]:
        patients = []
        neoantigens = []
        for _ in range(num_patients):
            patient = self.patient_provider.patient()
            patients.append(patient)
            for _ in range(num_neoantigens_per_patient):
                neoantigen = self.neoantigen_provider.neoantigen(
                    patient_identifier=patient.identifier, wildtype=wildtype)
                neoantigens.append(neoantigen)
        return patients, neoantigens
