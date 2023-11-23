# Input data

## General information

NeoFox requires two input files: a candidate file with neoantigen or neoepitope candidates and a file with patient data.
Both the input files must be tab-separated. The file with neoantigen candidates can be provided either in tabular format or in JSON format and this file may contain 
additional user-specific input that will be kept during the annotation process. The patient file requires a tabular format.

Alternatively, NeoFox may annotate a set of neoepitope candidates for which it will require a file with neoepitope
candidates and optionally a file with patient data. Both files are required in tabular format.

## Candidate file 

### Tabular file format

#### Neoantigen candidates  

This is an dummy example of a table with neoantigen candidates in tabular format:  

| gene  | wildTypeXmer                | mutatedXmer                | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | dnaVariantAlleleFrequency | external_annotation_1 | external_annotation_2 |
|-------|-----------------------------|----------------------------|-------------------|---------------|---------------------------|---------------------------|-----------------------|-----------------------|
| BRCA2 | AAAAAAAAAAAAALAAAAAAAAAAAAA | AAAAAAAAAAAAAFAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAMAAAAAAAAAAAAA | AAAAAAAAAAAAARAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAGAAAAAAAAAAAAA | AAAAAAAAAAAAAKAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAACAAAAAAAAAAAAA | AAAAAAAAAAAAAEAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAKAAAAAAAAAAAAA | AAAAAAAAAAAAACAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |

where:
- `gene`: the HGNC gene symbol. (This field is not required for neoantigen candidates derived from other sources than SNVs)      
- `mutatedXmer`: the neoantigen candidate sequence, i.e. the mutated amino acid sequence. In case of SNVs, the mutation should be located in the middle. We advise that the point mutation is flanked by 13 amino acid on both sites (IUPAC 1 respecting casing, eg: A) to cover both MHC I and MHC II neopeptides
- `wildTypeXmer`: the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A). This field shall be empty, specially in the case of neoantigen candidates derived from other sources than SNVs.  
- `patientIdentifier`: the patient identifier
- `rnaExpression`: RNA expression. (**optional**) (see *NOTE*) This value can be in any format chosen by the user (e.g. TPM, RPKM) but it is recommended to be consistent for data that should be compared.
- `rnaVariantAlleleFrequency`: the variant allele frequency (VAF) calculated from the RNA (**optional**)
- `dnaVariantAlleleFrequency`: the VAF calculated from the DNA. (**optional**)

**NOTE:** 

- Neofox annotates gene expression in TCGA cohort indicated in the `tumorType` in the patient data (see below) which might be helpful if rnaExpression is unknown. Please, not that this does not work for mouse data. 


#### Neoepitope candidates 

This is an dummy example of a table with neoepitope candidates in tabular format:  

| gene  | mutatedPeptide      | wildTypePeptide             | alleleMhcI  | isoformMhcII | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | dnaVariantAlleleFrequency | 
|-------|---------------------|-----------------------------|-------------|--------------|-------------------|---------------------------|---------------------------|---------------------------|
| BRCA2 | AAAALAAAAA | AAAAFAAAAA                  | HLA-A*01:01 |              | Ptx               | 7.942         | 0.85                      | 0.34                      |
| BRCA2 | AAAAAAAAAAAAAMAAAAAAAAAAAAA | AAAAAAAAAAAAARAAAAAAAAAAAAA |          | DRB1*01:01   | Ptx               | 7.942         | 0.85                      | 0.34                      |
| BRCA2 | AAAAGAAAAA | AAAAKAAAAA                  |          |              | Ptx               | 7.942         | 0.85                      | 0.34                      |
| BRCA2 | AAAAAAAAAAAAACAAAAAAAAAAAAA | AAAAAAAAAAAAAEAAAAAAAAAAAAA |          |              | Ptx               | 7.942         | 0.85                      | 0.34                      |
| BRCA2 | AAAAAAAAAAAAAKAAAAAAAAAAAAA | AAAAAAAAAAAAACAAAAAAAAAAAAA |          |              | Ptx               | 7.942         | 0.85                      | 0.34                      |

where:
- `mutatedPeptide`: the neoepitope candidate sequence, i.e. the mutated amino acid sequence. MHC-I neoepitopes should have a length between 8 and 14 amino acids, MHC-II neoepitopes should have a length between 9 and 20000 amino acids. 
- `wildTypePeptide`: the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A). This field shall be empty, specially in the case of neoepitope candidates derived from other sources than SNVs.  
- `alleleMhcI`: the MHC-I allele to which this neoepitope is linked (**optional**)
- `isoformMhcII`: the MHC-II isoform to which this neoepitope is linked (**optional**)
- `patientIdentifier`: the patient identifier (**only required if alleleMhcI and isoformMhcII are not provided**)
- `gene`: the HGNC gene symbol. (This field is optional)      
- `rnaExpression`: RNA expression. (**optional**) (see *NOTE*) This value can be in any format chosen by the user (e.g. TPM, RPKM) but it is recommended to be consistent for data that should be compared.
- `rnaVariantAlleleFrequency`: the variant allele frequency (VAF) calculated from the RNA (**optional**)
- `dnaVariantAlleleFrequency`: the VAF calculated from the DNA. (**optional**)

**NOTE:** 

- Neoepitopes with a value for `alleleMhcI` are considered MHC-I neoepitopes, likewise neoepitopes with a value for `isoformMhcII` are considered MHC-II neoepitopes. Both fields cannot be provided for the same neoepitope.
- If none of `alleleMhcI` and `isoformMhcII` are provided then the `patientIdentifier` is required and one neoepitope sharing the same sequence will be annotated for each MHC-I allele and MHC-II isoform according to the patient HLA type.
- Neofox annotates gene expression in TCGA cohort indicated in the `tumorType` in the patient data (see below) which might be helpful if rnaExpression is unknown. Please, not that this does not work for mouse data.


### JSON file format

#### Neoantigen candidates  

Besides tabular format, neoantigen candidates can be provided as a list of neoantigen models in JSON format as shown below. To simplify, only one full neoantigen model is shown. The terminology follows the descriptions for the [tabular file format](#tabular-file-format). For a more detailed description of the models, please refer to [here](05_models.md):  

```json
[{
    "patientIdentifier": "Ptx",
    "gene": "BRCA2",
    "wildTypeXmer": "AAAAAAAAAAAAALAAAAAAAAAAAAA",
    "mutatedXmer": "AAAAAAAAAAAAAFAAAAAAAAAAAAA"
}]
``` 

#### Neoepitope candidates 

Not supported at the moment.

## File with patient data

### Human

This is an dummy example of a patient file in tabular format:  

| identifier | mhcIAlleles                                                                  | mhcIIAlleles                                                                                                                                                   | tumorType |
|------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
| Ptx        | HLA-A\*03:01,HLA-A\*29:02,HLA-B\*07:02,HLA-B\*44:03,HLA-C\*07:02,HLA-C*16:01 | HLA-DRB1\*03:01,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*05:01,HLA-DQB1\*01:01,HLA-DQB1\*04:02,HLA-DPA1\*01:03,HLA-DPA1\*03:01,HLA-DPB1\*13:01,HLA-DPB1*04:02 | HNSC      |
| Pty        | HLA-A\*02:01,HLA-A\*30:01,HLA-B\*07:34,HLA-B\*44:03,HLA-C\*07:02,HLA-C*07:02 | HLA-DRB1\*04:02,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*04:01,HLA-DQB1\*03:02,HLA-DQB1\*14:01,HLA-DPA1\*01:03,HLA-DPA1\*02:01,HLA-DPB1\*02:01,HLA-DPB1*04:01 | HNSC      |

where:
- `identifier`: the patient identifier
- `mhcIAlleles`: comma separated MHC I alleles of the patient for HLA-A, HLA-B and HLA-C. If homozygous, the allele should be added twice.
- `mhcIIAlleles`: comma separated  MHC II alleles of the patient for HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1 and HLA-DPB1. If homozygous, the allele should be added twice.  
- `tumorType`: tumour entity in TCGA study abbreviation format (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). This field is required for imputation with gene expression and at the moment the following tumor types are supported:

| Study Name                                                         | Abbreviation |
|--------------------------------------------------------------------|-------------------|
| Adrenocortical carcinoma                                           | ACC               |
| Bladder Urothelial Carcinoma                                       | BLCA              |
| Breast invasive carcinoma                                          | BRCA              |
| Cervical squamous cell carcinoma and endocervical adenocarcinoma   | CESC              |
| Cholangiocarcinoma                                                 | CHOL              |
| Colon adenocarcinoma                                               | COAD              |
| Esophageal carcinoma                                               | ESCA              |
| Glioblastoma multiforme                                            | GBM               |
| Head and Neck squamous cell carcinoma                              | HNSC              |
| Kidney Chromophobe                                                 | KICH              |
| Kidney renal papillary cell carcinoma                              | KIRP              |
| Liver hepatocellular carcinoma                                     | LIHC              |
| Lung adenocarcinoma                                                | LUAD              |
| Lung squamous cell carcinoma                                       | LUSC              |
| Ovarian serous cystadenocarcinoma                                  | OV                |
| Pancreatic adenocarcinoma                                          | PAAD              |
| Prostate adenocarcinoma                                            | PRAD              |
| Rectum adenocarcinoma                                              | READ              |
| Sarcoma                                                            | SARC              |
| Skin Cutaneous Melanoma                                            | SKCM              |
| Testicular Germ Cell Tumors                                        | TGCT              |
| Uterine Corpus Endometrial Carcinoma                               | UCEC              |


### Mouse

This is a dummy example of a "patient" file in tabular format for mouse:

| identifier | mhcIAlleles                                                                  | mhcIIAlleles                                                                                                                                                   |
|------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Ptz        | H2Db,H2Db,H2Kb,H2Kb,H2Lb,H2Lb | H2Ab,H2Ab,H2Eb,H2Eb |


**WARNING**: NeoFox requires MHC alleles in homozygosity to be provided twice, also for mouse. 
Otherwise they are considered as hemizygous. 
For instance, each gene would be interpreted as hemizygous when `H2Db,H2Kb,H2Lb` is provided.
In the case of inbred mouse strains the MHC alleles are homozygous state in most of cases.

While using NeoFox to annotate neoantigen candidates from mouse, it should be considered that the nomenclature of the major histocompatibility complex for Mus musculus, 
H-2, is not described with the same level of detail as in humans (HLA).  
The Mus musculus strains used in laboratory experimentation are inbred mouse strains with limited variability. 
NetMHCpan and netMHCIIpan, and by extension NeoFox, support a subset of the H-2 alleles found in laboratory mice which
is again a small subset of the wild type.
As a consequence, the MHC II nomenclature is highly simplified. The alpha and beta chain genes are always considered to be 
part of the same haplotype. Furthermore, only homozygosity is considered. 
Thus there is only one possible MHC II isoform for each pair of genes, as opposed to four in human.

Murine H-2 alleles are not standardized as the human HLA alleles are by the WHO Nomenclature Committee for Factors of the HLA System.  
H-2 alleles are represented as follows in NeoFox:  
- **MHC I** genes K, D and L: H2K, H2D and H2L
- **MHC II** genes A and E: H2A and H2E  

A given allele is represented by a last small case single letter (eg: d, k, p) with an optional number (eg: d1, p2)  .

These are examples of H-2 alleles: H2Kd, H2Dd, H2Lp 




