import os
from typing import List, Tuple
from faker import Faker

from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import Patient, Neoantigen, Mhc2Isoform, MhcAllele
from neofox.model.wrappers import get_mhc2_isoform_name
from neofox.references.references import ReferenceFolder, HOMO_SAPIENS_FASTA, DependenciesConfiguration
from neofox.tests.synthetic_data.factories import PatientProvider, NeoantigenProvider


class DataGenerator:

    def __init__(self, reference_folder: ReferenceFolder, configuration: DependenciesConfiguration):
        faker = Faker()
        mixmhcpred_alleles = set(self.load_mhc1_alleles(
            MixMHCpred(None, configuration=configuration).available_alleles))
        netmhcpan_alleles = set(self.load_mhc1_alleles(
            reference_folder.get_available_alleles().get_available_mhc_i()))
        mhc1_alleles = mixmhcpred_alleles.union(netmhcpan_alleles)

        mixmhc2pred_alleles = set(self.load_mhc2_alleles(
            MixMhc2Pred(None, configuration=configuration).available_alleles, fix=False))
        netmhc2pan_alleles = set(self.load_mhc2_alleles(
            reference_folder.get_available_alleles().get_available_mhc_ii(), fix=True))
        mhc2_isoforms = mixmhc2pred_alleles.union(netmhc2pan_alleles)

        self.patient_provider = PatientProvider(faker, mhc1_alleles, mhc2_isoforms)
        self.neoantigen_provider = NeoantigenProvider(
            faker, proteome_fasta=os.path.join(reference_folder.proteome_db, HOMO_SAPIENS_FASTA))

    def load_mhc1_alleles(self, available_alleles: List[str]):
        mhc_alleles = []
        for a in available_alleles:
            try:
                parsed_allele = ModelConverter.parse_mhc_allele(self.dirty_fix_mhc_representation(a))
            except AssertionError:
                continue
            mhc_alleles.append(parsed_allele.name)
        return mhc_alleles

    def load_mhc2_alleles(self, available_alleles: List[str], fix=False):
        mhc_alleles = []
        for a in available_alleles:
            try:
                parsed_allele = self.parse_mhc2_isoform(a, fix=fix)
            except AssertionError:
                continue
            mhc_alleles.append(parsed_allele.name)
        return mhc_alleles

    # TODO: temporary fix valid only for netmhc2pan alleles until the homonimous method in conversion.py is fixed
    def parse_mhc2_isoform(self, isoform: str, fix=False) -> Mhc2Isoform:
        # TODO: this method currently fails for netmhc2pan alleles which are like 'HLA-DQA10509-DQB10630'
        # infers gene, group and protein from the name
        isoform = isoform.strip("HLA-")
        if "DQA" in isoform or "DPA" in isoform:
            alpha_chain = ModelConverter.parse_mhc_allele(
                self.dirty_fix_mhc_representation(isoform.split("-")[0]) if fix else isoform.split("__")[0])
            beta_chain = ModelConverter.parse_mhc_allele(
                self.dirty_fix_mhc_representation(isoform.split("-")[1]) if fix else isoform.split("__")[1])
        else:
            alpha_chain = MhcAllele()
            beta_chain = ModelConverter.parse_mhc_allele(
                self.dirty_fix_mhc_representation(isoform) if fix else isoform)
        # builds the final allele representation and validates it just in case
        name = get_mhc2_isoform_name(alpha_chain, beta_chain)
        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)

    def dirty_fix_mhc_representation(self, allele: str):
        # not intended to use in production this is a fix for data simulation while we integrate a better support
        # of ambiguous HLA alleles
        return allele[0:-2] + ":" + allele[-2:]

    def generate_data(self, num_patients, num_neoantigens_per_patient) -> Tuple[List[Patient], List[Neoantigen]]:
        patients = []
        neoantigens = []
        for _ in range(num_patients):
            patient = self.patient_provider.patient()
            patients.append(patient)
            for _ in range(num_neoantigens_per_patient):
                neoantigen = self.neoantigen_provider.neoantigen(patient_identifier=patient.identifier)
                neoantigens.append(neoantigen)
        return patients, neoantigens
