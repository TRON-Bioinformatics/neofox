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
    """*Valid names for MHC I classic genes"""

    A = 0
    B = 1
    C = 2


class Mhc2GeneName(betterproto.Enum):
    """
    *Valid names for MHC II classic genes.DRA is not included in this list as
    it does not have much variability in the population and for our purpose
    isconsidered constant.
    """

    DRB1 = 0
    DPA1 = 1
    DPB1 = 2
    DQA1 = 3
    DQB1 = 4


class Mhc2Name(betterproto.Enum):
    """*Valid names for MHC II classic molecules"""

    DR = 0
    DP = 1
    DQ = 2


@dataclass
class Mutation(betterproto.Message):
    # *The amino acid position within the neoantigen candidate sequence. 1-based,
    # starting in the N-terminus
    position: List[int] = betterproto.int32_field(1)
    # *Amino acid sequence of the WT corresponding to the neoantigen candidate
    # sequence (IUPAC 1 letter codes)
    wild_type_xmer: str = betterproto.string_field(2)
    # *Amino acid sequence of the neoantigen candidate (IUPAC 1 letter codes)
    mutated_xmer: str = betterproto.string_field(3)


@dataclass
class Annotation(betterproto.Message):
    """*This is a generic class to hold annotations from Neofox"""

    # *The name of the annotation
    name: str = betterproto.string_field(1)
    # *The value of the annotation
    value: str = betterproto.string_field(2)


@dataclass
class NeoantigenAnnotations(betterproto.Message):
    """*A set of annotations for a neoantigen"""

    # *List of annotations
    annotations: List["Annotation"] = betterproto.message_field(1)
    # *The annotator
    annotator: str = betterproto.string_field(2)
    # *The version of the annotator
    annotator_version: str = betterproto.string_field(3)
    # *A timestamp determined when the annotation was created
    timestamp: str = betterproto.string_field(4)
    # *Annotation resources MD5 hash
    resources_hash: str = betterproto.string_field(5)


@dataclass
class Neoantigen(betterproto.Message):
    """*A neoantigen minimal definition"""

    # *Patient identifier
    patient_identifier: str = betterproto.string_field(1)
    # *The HGNC gene symbol or gene identifier
    gene: str = betterproto.string_field(2)
    # *The mutation
    mutation: "Mutation" = betterproto.message_field(3)
    # *Expression value of the transcript from RNA data. Range [0, +inf].
    rna_expression: float = betterproto.float_field(4)
    # *Expression value of the transcript from TCGA data. Range [0, +inf].
    imputed_gene_expression: float = betterproto.float_field(5)
    # *Variant allele frequency from the DNA. Range [0.0, 1.0]
    dna_variant_allele_frequency: float = betterproto.float_field(6)
    # *Variant allele frequency from the RNA. Range [0.0, 1.0]
    rna_variant_allele_frequency: float = betterproto.float_field(7)
    # *The NeoFox neoantigen annotations
    neofox_annotations: "NeoantigenAnnotations" = betterproto.message_field(8)
    # *List of external annotations
    external_annotations: List["Annotation"] = betterproto.message_field(9)


@dataclass
class Patient(betterproto.Message):
    """
    *The metadata required for analysis for a given patient + its patient
    identifier
    """

    # *Patient identifier
    identifier: str = betterproto.string_field(1)
    # *Is RNA expression available?
    is_rna_available: bool = betterproto.bool_field(2)
    # *Tumor entity in TCGA study abbrevation style as described here:
    # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-
    # abbreviations
    tumor_type: str = betterproto.string_field(3)
    # *MHC I classic molecules
    mhc1: List["Mhc1"] = betterproto.message_field(4)
    # *MHC II classic molecules
    mhc2: List["Mhc2"] = betterproto.message_field(5)


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
