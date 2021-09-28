# Input data

## General information

NeoFox requires two input files: a file with neoantigen candidates derived from point mutations and a file with patient data. 
The file with neoantigen candidates can be provided either in tabular format or in JSON format and this file may contain 
additional user-specific input that will be kept during the annotation process. The patient file is required in tabular format.

## File with neoantigen candidates

#### Tabular file format  

This is an dummy example of a table with neoantigen candidates in tabular format:  

| gene  | mutation.wildTypeXmer       | mutation.mutatedXmer        | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | dnaVariantAlleleFrequency | external_annotation_1 | external_annotation_2 |
|-------|-----------------------------|-----------------------------|-------------------|---------------|---------------------------|---------------------------|-----------------------|-----------------------|
| BRCA2 | AAAAAAAAAAAAALAAAAAAAAAAAAA | AAAAAAAAAAAAAFAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAMAAAAAAAAAAAAA | AAAAAAAAAAAAARAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAGAAAAAAAAAAAAA | AAAAAAAAAAAAAKAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAACAAAAAAAAAAAAA | AAAAAAAAAAAAAEAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAKAAAAAAAAAAAAA | AAAAAAAAAAAAACAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |

where:
- `gene`: the HGNC gene symbol. (This field is not required for neoantigen candidates derived from other sources than SNVs)      
- `mutation.mutatedXmer`: the neoantigen candidate sequence, i.e. the mutated amino acid sequence. In case of SNVs, the mutation should be located in the middle. We advise that the point mutation is flanked by 13 amino acid on both sites (IUPAC 1 respecting casing, eg: A) to cover both MHC I and MHC II neoepitopes
- `mutation.wildTypeXmer`: the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A). This field shall be empty, specially in the case of neoantigen candidates derived from other sources than SNVs.  
- `patientIdentifier`: the patient identifier
- `rnaExpression`: RNA expression. (**optional**) (see *NOTE*) This value can be in any format chosen by the user (e.g. TPM, RPKM) but it is recommended to be consistent for data that should be compared.
- `rnaVariantAlleleFrequency`: the variant allele frequency calculated from the RNA (Should be empty or NA if no value available, this will be estimated using the `dnaVariantAlleleFrequency` if not available)
- `dnaVariantAlleleFrequency`: the variant allele frequency calculated from the DNA. Should be empty or NA if no value available.

**NOTE:** If rnaExpression is not provided, expression will be estimated by gene expression in TCGA cohort indicated in the `tumorType` in the patient data (see below). 

### JSON file format

Besides tabular format, neoantigen candidates can be provided as a list of neoantigen models in JSON format as shown below. To simplify, only one full neoantigen model is shown. The terminology follows the descriptions for the [tabular file format](#tabular-file-format). For a more detailed description of the models, please refer to [here](05_models.md):  

```json
[{
    "patientIdentifier": "Ptx",
    "gene": "BRCA2",
    "mutation": {
        "wildTypeXmer": "AAAAAAAAAAAAALAAAAAAAAAAAAA",
        "mutatedXmer": "AAAAAAAAAAAAAFAAAAAAAAAAAAA"
    }
}]
``` 

## File with patient information

This is an dummy example of a patient file in tabular format:  

| identifier | mhcIAlleles                                                                  | mhcIIAlleles                                                                                                                                                   | tumorType |
|------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
| Ptx        | HLA-A\*03:01,HLA-A\*29:02,HLA-B\*07:02,HLA-B\*44:03,HLA-C\*07:02,HLA-C*16:01 | HLA-DRB1\*03:01,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*05:01,HLA-DQB1\*01:01,HLA-DQB1\*04:02,HLA-DPA1\*01:03,HLA-DPA1\*03:01,HLA-DPB1\*13:01,HLA-DPB1*04:02 | HNSC      |
| Pty        | HLA-A\*02:01,HLA-A\*30:01,HLA-B\*07:34,HLA-B\*44:03,HLA-C\*07:02,HLA-C*07:02 | HLA-DRB1\*04:02,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*04:01,HLA-DQB1\*03:02,HLA-DQB1\*14:01,HLA-DPA1\*01:03,HLA-DPA1\*02:01,HLA-DPB1\*02:01,HLA-DPB1*04:01 | HNSC      |

where:
- `identifier`: the patient identifier
- `mhcIAlleles`: comma separated MHC I alleles of the patient for HLA-A, HLA-B and HLA-C. If homozygous, the allele should be added twice.
- `mhcIIAlleles`: comma separated  MHC II alleles of the patient for HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1 and HLA-DPB1. If homozygous, the allele should be added twice.  
- `tumorType`: tumour entity in TCGA study abbreviation format (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). This field is required for expression imputation and at the moment the following tumor types are supported:

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


An example for Mus musculus with H-2 alleles (note that the tumour type is not applicable):

| identifier | mhcIAlleles                                                                  | mhcIIAlleles                                                                                                                                                   |
|------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Ptz        | H2Db,H2Db,H2Kb,H2Kb,H2Lb,H2Lb | H2Ab,H2Ab,H2Eb,H2Eb |


**WARNING**: alleles in homozygous state need to be provided twice, otherwise they are considered as hemizygous. 
 For instance `["HLA-A*01:01"]` would be interpreted as hemizygous and `["HLA-A*01:01", "HLA-A*01:01"]` as homozygous.
In the case of Mus musculus strains you will want to provide genes in homozygous state in most of cases.

