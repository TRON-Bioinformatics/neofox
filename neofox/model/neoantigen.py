# Generated by the protocol buffer compiler.  DO NOT EDIT!
# sources: neoantigen.proto
# plugin: python-betterproto
from dataclasses import dataclass
from typing import List

import betterproto


class Zygosity(betterproto.Enum):
    """*The zygosity of a given gene"""

    # *Two equal copies of the gene
    HOMOZYGOUS = 0
    # *Two different copies of the gene
    HETEROZYGOUS = 1
    # *Only one copy of the gene
    HEMIZYGOUS = 2
    # *No copy of the gene
    LOSS = 3


class Mhc1Name(betterproto.Enum):
    """
    *Valid names for MHC I classic genesMus musculus gene names are preceded by
    the prefix H2 to avoid naming collisions.
    """

    # Homo sapiens
    A = 0
    B = 1
    C = 2
    # Mus musculus
    H2K = 3
    H2D = 4
    H2L = 5


class Mhc2GeneName(betterproto.Enum):
    """
    *Valid names for MHC II classic genes.DRA is not included in this list as
    it does not have much variability in the population and for our purpose
    isconsidered constant.For Mus musculus we do not represent alpha and beta
    chains as they are homozygotes at all their MHC loci.Hence, they can be
    treated as a single gene, like DR is for HLA.See http://www.imgt.org/IMGTre
    pertoireMH/Polymorphism/haplotypes/mouse/MHC/Mu_haplotypes.htmlMus musculus
    gene names are preceded by the prefix H2 to avoid naming collisions.
    """

    # Homo sapiens
    DRB1 = 0
    DPA1 = 1
    DPB1 = 2
    DQA1 = 3
    DQB1 = 4
    # Mus musculus
    H2A = 5
    H2E = 6


class Mhc2Name(betterproto.Enum):
    """*Valid names for MHC II classic molecules"""

    DR = 0
    DP = 1
    DQ = 2
    H2A_molecule = 3
    H2E_molecule = 4


@dataclass
class Annotation(betterproto.Message):
    """*This is a generic class to hold annotations from Neofox"""

    # *The name of the annotation
    name: str = betterproto.string_field(1)
    # *The value of the annotation
    value: str = betterproto.string_field(2)


@dataclass
class Resource(betterproto.Message):
    """*This is a class to track the version of an annotation resource"""

    # *The name of the resource
    name: str = betterproto.string_field(1)
    # *The version of the resource
    version: str = betterproto.string_field(2)
    # *The URL of the resource if applicable
    url: str = betterproto.string_field(3)
    # *The MD5 hash of the resource if applicable. This may be used when version
    # is not available
    hash: str = betterproto.string_field(4)
    # *The timestamp when the download happened
    download_timestamp: str = betterproto.string_field(5)


@dataclass
class Annotations(betterproto.Message):
    """*A set of annotations for a neoantigen candidate"""

    # *List of annotations
    annotations: List["Annotation"] = betterproto.message_field(1)
    # *The annotator
    annotator: str = betterproto.string_field(2)
    # *The version of the annotator
    annotator_version: str = betterproto.string_field(3)
    # *A timestamp determined when the annotation was created
    timestamp: str = betterproto.string_field(4)
    # *List of resources
    resources: List["Resource"] = betterproto.message_field(5)


@dataclass
class Patient(betterproto.Message):
    """
    *The metadata required for analysis for a given patient + its patient
    identifier
    """

    # *Patient identifier
    identifier: str = betterproto.string_field(1)
    # *Tumor entity in TCGA study abbrevation style as described here:
    # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-
    # abbreviations
    tumor_type: str = betterproto.string_field(2)
    # *MHC I classic molecules
    mhc1: List["Mhc1"] = betterproto.message_field(3)
    # *MHC II classic molecules
    mhc2: List["Mhc2"] = betterproto.message_field(4)


@dataclass
class Mhc1(betterproto.Message):
    """
    *Models MHC I alleles related to the same MHC I gene, i.e. 2 alleles/2
    isoforms per gene
    """

    # *MHC I gene name
    name: "Mhc1Name" = betterproto.enum_field(1)
    # *Zygosity of the gene
    zygosity: "Zygosity" = betterproto.enum_field(2)
    # *The alleles of the gene (0, 1 or 2)
    alleles: List["MhcAllele"] = betterproto.message_field(3)


@dataclass
class Mhc2(betterproto.Message):
    """
    *Models MHC II alleles related to the same MHC II protein, i.e. 4 isoforms
    related to 2 genes with 2 alleles each
    """

    # *MHC II molecule name
    name: "Mhc2Name" = betterproto.enum_field(1)
    # *List of MHC II genes
    genes: List["Mhc2Gene"] = betterproto.message_field(2)
    # *Different combinations of MHC II alleles building different isoforms
    isoforms: List["Mhc2Isoform"] = betterproto.message_field(3)


@dataclass
class Mhc2Isoform(betterproto.Message):
    """*MHC II isoform"""

    # *Name to refer to the MHC II isoform
    name: str = betterproto.string_field(1)
    # *The alpha chain of the isoform
    alpha_chain: "MhcAllele" = betterproto.message_field(2)
    # *The beta chain of the isoform
    beta_chain: "MhcAllele" = betterproto.message_field(3)


@dataclass
class Mhc2Gene(betterproto.Message):
    """*MHC II gene"""

    # *MHC II gene name
    name: "Mhc2GeneName" = betterproto.enum_field(1)
    # *Zygosity of the gene
    zygosity: "Zygosity" = betterproto.enum_field(2)
    # *The alleles of the gene (0, 1 or 2)
    alleles: List["MhcAllele"] = betterproto.message_field(3)


@dataclass
class MhcAllele(betterproto.Message):
    """
    *MHC allele representation. It does not include non synonymous changes to
    the sequence, changes in the non coding regionor changes in expression. See
    http://hla.alleles.org/nomenclature/naming.html for details
    """

    # *HLA full name as provided by the user (e.g.: HLA-DRB1*13:01:02:03N). This
    # will be parsed into name, gene and group.Any digit format is allowed for
    # this field (ie: 4, 6 or 8 digits), 2 digits names are not specific enough
    # for ourpurpose and thus invalid
    full_name: str = betterproto.string_field(1)
    # *A specific HLA protein (e.g. HLA-DRB1*13:01). Alleles whose numbers differ
    # in group and protein must differ in oneor more nucleotide substitutions
    # that change the amino acid sequence of the encoded protein.This name is
    # normalized to avoid different representations of the same allele. For
    # instance both HLA-DRB113:01 andHLA-DRB1*13:01:02:03N will be transformed
    # into their normalised version HLA-DRB1*13:01. This name is also truncatedto
    # 4 digits. 2 digits names are not specific enough for our purpose and thus
    # invalid
    name: str = betterproto.string_field(2)
    # *The gene from either MHC I or II (e.g. DRB1, A) (this information is
    # redundant with the Mhc1Gene.name andMhc2Gene.name but it is convenient to
    # have this at this level too, code will check for data coherence)
    gene: str = betterproto.string_field(3)
    # *A group of alleles defined by a common serotype ie: Serological antigen
    # carried by an allotype (e.g. 13 from HLA-DRB1*13)
    group: str = betterproto.string_field(4)
    # *A specific protein (e.g.: 02 from HLA-DRB1*13:02)
    protein: str = betterproto.string_field(5)


@dataclass
class PredictedEpitope(betterproto.Message):
    # *Not sure that we need this... this is in the old PredictedEpitope model
    position: int = betterproto.int32_field(1)
    # *The mutated peptide
    mutated_peptide: str = betterproto.string_field(2)
    # *Closest wild type peptide
    wild_type_peptide: str = betterproto.string_field(3)
    # *MHC I allele
    allele_mhc_i: "MhcAllele" = betterproto.message_field(4)
    # *MHC II isoform
    isoform_mhc_i_i: "Mhc2Isoform" = betterproto.message_field(5)
    # *MHC binding affinity for the mutated peptide. This value is estimated with
    # NetMHCpan in case of MHC-I peptidesand NetMHCIIpan in cas of MHC-II
    # peptides
    affinity_mutated: float = betterproto.float_field(6)
    # *MHC binding rank for the mutated peptide. This value is estimated with
    # NetMHCpan in case of MHC-I peptidesand NetMHCIIpan in cas of MHC-II
    # peptides
    rank_mutated: float = betterproto.float_field(7)
    # *MHC binding affinity for the wild type peptide. This value is estimated
    # with NetMHCpan in case of MHC-I peptidesand NetMHCIIpan in cas of MHC-II
    # peptides
    affinity_wild_type: float = betterproto.float_field(8)
    # *MHC binding rank for the wild type peptide. This value is estimated with
    # NetMHCpan in case of MHC-I peptidesand NetMHCIIpan in cas of MHC-II
    # peptides
    rank_wild_type: float = betterproto.float_field(9)
    # *The NeoFox neoantigen annotations
    neofox_annotations: "Annotations" = betterproto.message_field(10)
    # *Patient identifier
    patient_identifier: str = betterproto.string_field(11)
    # *The HGNC gene symbol or gene identifier
    gene: str = betterproto.string_field(12)
    # *Expression value of the transcript from RNA data. Range [0, +inf].
    rna_expression: float = betterproto.float_field(13)
    # *Expression value of the transcript from TCGA data. Range [0, +inf].
    imputed_gene_expression: float = betterproto.float_field(14)
    # *Variant allele frequency from the DNA. Range [0.0, 1.0]
    dna_variant_allele_frequency: float = betterproto.float_field(15)
    # *Variant allele frequency from the RNA. Range [0.0, 1.0]
    rna_variant_allele_frequency: float = betterproto.float_field(16)
    # *External annotations for neoepitope mode.
    external_annotations: List["Annotation"] = betterproto.message_field(17)

@dataclass
class Neoantigen(betterproto.Message):
    """*A neoantigen minimal definition"""

    # *Patient identifier
    patient_identifier: str = betterproto.string_field(1)
    # *The HGNC gene symbol or gene identifier
    gene: str = betterproto.string_field(2)
    # *The amino acid position within the neoantigen candidate sequence. 1-based,
    # starting in the N-terminus
    position: List[int] = betterproto.int32_field(3)
    # *Amino acid sequence of the WT corresponding to the neoantigen candidate
    # sequence (IUPAC 1 letter codes)
    wild_type_xmer: str = betterproto.string_field(4)
    # *Amino acid sequence of the neoantigen candidate (IUPAC 1 letter codes)
    mutated_xmer: str = betterproto.string_field(5)
    # *Expression value of the transcript from RNA data. Range [0, +inf].
    rna_expression: float = betterproto.float_field(6)
    # *Expression value of the transcript from TCGA data. Range [0, +inf].
    imputed_gene_expression: float = betterproto.float_field(7)
    # *Variant allele frequency from the DNA. Range [0.0, 1.0]
    dna_variant_allele_frequency: float = betterproto.float_field(8)
    # *Variant allele frequency from the RNA. Range [0.0, 1.0]
    rna_variant_allele_frequency: float = betterproto.float_field(9)
    # *The NeoFox neoantigen annotations
    neofox_annotations: "Annotations" = betterproto.message_field(10)
    # *List of external annotations
    external_annotations: List["Annotation"] = betterproto.message_field(11)
    # *List of predicted neoepitopes for MHC-I with feature annotation (optional)
    neoepitopes_mhc_i: List["PredictedEpitope"] = betterproto.message_field(12)
    # *List of predicted neoepitopes for MHC-II with feature annotation
    # (optional)
    neoepitopes_mhc_i_i: List["PredictedEpitope"] = betterproto.message_field(13)
