from typing import List

from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACData
from faker.providers.address import Provider

from neofox.exceptions import NeofoxDataValidationException
from neofox.expression_imputation.expression_imputation import ExpressionAnnotator
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.conversion import ModelConverter, ModelValidator
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Patient, Mhc1Name, Neoantigen, Mutation, Mhc2Name, Mhc2Isoform, \
    MhcAllele
from neofox.model.wrappers import get_mhc2_isoform_name
from neofox.references.references import HlaDatabase


class PatientProvider(Provider):

    def __init__(self, generator, mhc1_alleles, mhc2_alleles, hla_database: HlaDatabase):
        Provider.__init__(self, generator)
        self.hla_database = hla_database
        self.mhc_parser = MhcParser(hla_database)
        # gets available alleles from netmhcpan and netmhc2pan
        self.available_mhc1_alleles = self.load_mhc1_alleles(mhc1_alleles)
        self.available_mhc2_isoforms = self.load_mhc2_isoforms(mhc2_alleles)
        # gets available tumor types
        self.available_tumor_types = ExpressionAnnotator().cohort_indices.keys()


    def load_mhc1_alleles(self, available_alleles: List[str]):
        mhc_alleles = []
        for a in available_alleles:
            try:
                parsed_allele = self.mhc_parser.parse_mhc_allele(a)
            except AssertionError:
                continue
            mhc_alleles.append(parsed_allele)
        return mhc_alleles

    def load_mhc2_isoforms(self, available_alleles: List[str]) -> List[Mhc2Isoform]:
        mhc_isoforms = []
        for a in available_alleles:
            try:
                parsed_allele = self.parse_mhc2_isoform(a)
            except AssertionError:
                continue
            mhc_isoforms.append(parsed_allele)
        return mhc_isoforms

    # TODO: temporary fix valid only for netmhc2pan alleles until the homonimous method in conversion.py is fixed
    def parse_mhc2_isoform(self, isoform: str) -> Mhc2Isoform:
        # TODO: this method currently fails for netmhc2pan alleles which are like 'HLA-DQA10509-DQB10630'
        # infers gene, group and protein from the name
        isoform = isoform.strip("HLA-")
        if "DQA" in isoform or "DPA" in isoform:
            alpha_chain = self.mhc_parser.parse_mhc_allele(isoform.split("-")[0])
            beta_chain = self.mhc_parser.parse_mhc_allele(isoform.split("-")[1])
        else:
            alpha_chain = MhcAllele()
            beta_chain = self.mhc_parser.parse_mhc_allele(isoform)
        # builds the final allele representation and validates it just in case
        name = get_mhc2_isoform_name(alpha_chain, beta_chain)
        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)

    def patient(self) -> Patient:

        patient = None
        found = False

        while not found:
            dr_isoforms = self.random_elements(self.get_hla_ii_alleles_by_gene(Mhc2Name.DR), unique=True, length=2)
            dp_isoforms = self.random_elements(self.get_hla_ii_alleles_by_gene(Mhc2Name.DP), unique=True, length=2)
            dq_isoforms = self.random_elements(self.get_hla_ii_alleles_by_gene(Mhc2Name.DQ), unique=True, length=2)

            # NOTE: for some reason some DP alleles are malformed and cause a validation error, most do not.
            # thus I retry until I get a valid combination of HLA alleles, will clarify in another reincarnation
            try:
                patient = Patient(
                    identifier=self.generator.unique.uuid4(),
                    is_rna_available=True,
                    tumor_type=self.random_elements(self.available_tumor_types, length=1)[0],
                    # by setting unique=True we enforce that all patients are heterozygous
                    mhc1=ModelConverter.parse_mhc1_alleles(
                        self.random_elements(self.get_hla_i_alleles_by_gene(Mhc1Name.A), unique=True, length=2) +
                        self.random_elements(self.get_hla_i_alleles_by_gene(Mhc1Name.B), unique=True, length=2) +
                        self.random_elements(self.get_hla_i_alleles_by_gene(Mhc1Name.C), unique=True, length=2),
                        self.hla_database
                    ),
                    mhc2=ModelConverter.parse_mhc2_alleles(
                        [i.alpha_chain.name for i in dp_isoforms] +
                        [i.beta_chain.name for i in dp_isoforms] +
                        [i.alpha_chain.name for i in dq_isoforms] +
                        [i.beta_chain.name for i in dq_isoforms] +
                        [i.beta_chain.name for i in dr_isoforms],
                        self.hla_database
                    )
                )
                patient = ModelValidator.validate_patient(patient)
            except NeofoxDataValidationException:
                continue
            found = True
        return patient

    def get_hla_i_alleles_by_gene(self, gene_name: Mhc1Name) -> List[str]:
        return [a.name for a in self.available_mhc1_alleles if a.gene == gene_name.name]

    def get_hla_ii_alleles_by_gene(self, gene_name: Mhc2Name) -> List[Mhc2Isoform]:
        return [a for a in self.available_mhc2_isoforms if gene_name.name in a.name]


class NeoantigenProvider(Provider):

    def __init__(self, generator, proteome_fasta, length_xmer=27):
        Provider.__init__(self, generator)
        # os.path.join(reference_folder.proteome_db, HOMO_SAPIENS_FASTA)
        self.length_xmer = length_xmer
        self.protein_list = self._load_proteome(proteome_fasta)
        self.aminoacids = list(IUPACData.protein_letters)

    def _load_proteome(self, proteome_fasta):
        # stores protein sequences of at least 27 AAs
        protein_list = []
        for record in SeqIO.parse(proteome_fasta, "fasta"):
            # TODO: accept protein with rare aminoacids when all cases are controlled
            if len(record.seq) >= self.length_xmer and not EpitopeHelper.contains_rare_amino_acid(record.seq):
                protein_list.append(str(record.seq))
        return protein_list

    def neoantigen(self, patient_identifier=None, wildtype=True) -> Neoantigen:

        neoantigen = None
        found = False
        while not found:
            try:
                neoantigen = Neoantigen(
                    identifier=self.generator.unique.uuid4(),
                    patient_identifier=self.generator.unique.uuid4() if patient_identifier is None else patient_identifier,
                    gene="BRCA2" if wildtype else None, # no gene if no wildtype provided
                    mutation=self.mutation(wildtype=wildtype),
                    rna_expression=float(self.random_number(digits=4, fix_len=True))/100,
                    dna_variant_allele_frequency=float(self.random_number(digits=3, fix_len=True))/1000,
                    rna_variant_allele_frequency=float(self.random_number(digits=3, fix_len=True))/1000
                )
                neoantigen = ModelValidator.validate_neoantigen(neoantigen)
            except NeofoxDataValidationException:
                continue
            found = True

        return neoantigen

    def mutation(self, wildtype) -> Mutation:
        wildtype_xmer = self._get_wild_type_xmer()
        mutation_position = int(self.length_xmer / 2)
        mutated_xmer = wildtype_xmer[0:mutation_position] + self._mutate_aminoacid(wildtype_xmer[mutation_position]) + \
                       wildtype_xmer[mutation_position + 1:]
        if wildtype:
            mutation = Mutation(mutated_xmer=mutated_xmer, wild_type_xmer=wildtype_xmer)
        else:
            mutation = Mutation(mutated_xmer=mutated_xmer)

        return mutation

    def _get_wild_type_xmer(self):
        random_protein = self.random_elements(self.protein_list, length=1)[0]
        random_index = self.random_int(0, len(random_protein) - self.length_xmer)
        return random_protein[random_index: random_index + self.length_xmer]

    def _mutate_aminoacid(self, aminoacid):
        found = False
        new_aminoacid = None
        while not found:
            new_aminoacid = self.random_elements(self.aminoacids, length=1)[0]
            if new_aminoacid != aminoacid:
                found = True
        return new_aminoacid
