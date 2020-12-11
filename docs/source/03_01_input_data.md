# Input data

## General information

NeoFox requires two input files: a file with neoantigen candidates derived from point mutations and a file with patient data. The file with neoantigen candidates can be provided either in tabular format or in JSON format and this file may contain additional user-specific input that will be kept during the annotation process. The patient file is required in tabular format.

## Tabular format

### File with neonatigen candidates

We allow two different tabular formats of the neoantigen candidate file: `model-file` or `candidate-file` format. The neoantigen candidate input file must follow one of these two. However, additional user-specific columns are allowed and will be kept during the annotation process.   

1. **model-file format**  

This is an dummy example of a table with neoantigen candidates in `model-file` format:  

| gene  | mutation.wildTypeXmer       | mutation.mutatedXmer        | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | dnaVariantAlleleFrequency | external_annotation_1 | external_annotation_2 |
|-------|-----------------------------|-----------------------------|-------------------|---------------|---------------------------|---------------------------|-----------------------|-----------------------|
| BRCA2 | AAAAAAAAAAAAALAAAAAAAAAAAAA | AAAAAAAAAAAAAFAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAMAAAAAAAAAAAAA | AAAAAAAAAAAAARAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAGAAAAAAAAAAAAA | AAAAAAAAAAAAAKAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAACAAAAAAAAAAAAA | AAAAAAAAAAAAAEAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAKAAAAAAAAAAAAA | AAAAAAAAAAAAACAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |

where:
- `gene`: the HGNC gene symbol   
- `mutation.mutatedXmer`: the neoantigen candidate sequence, i.e. the mutated amino acid sequence. The mutation should be located in the middle, flanked by 13 amino acid on both sites (IUPAC 1 respecting casing, eg: A)
- `mutation.wildTypeXmer`: the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A)
- `patientIdentifier`: the patient identifier
- `rnaExpression`: the transcript expression. Should be empty if no value available
- `rnaVariantAlleleFrequency`: the variant allele frequency calculated from the RNA (**optional**, this will be estimated using the `dnaVariantAlleleFrequency` if not available)
- `dnaVariantAlleleFrequency`: the variant allele frequency calculated from the DNA (**optional**)  



2. **candidate-file format**  

Alternatively, neoantigen candidates can be provided in `candidate-file` format. In principle the columns are the same as in the `model-file`. Of note,`candidate-file` allows for an optional patient id in the data table. This is an dummy example:  

|     patient |     gene  | substitution |     transcript_expression |     +-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) |     [WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) |     VAF_in_tumor |     VAF_in_RNA    |
|-------------|-----------|--------------|---------------------------|--------------------------------------------|-------------------------------------------------|------------------|-------------------|
|     Ptx     |     BRCA2 | I547T        |     0.51950689            |     AAAAAAAAAAAAAFAAAAAAAAAAAAA            |     AAAAAAAAAAAAALAAAAAAAAAAAAA                 |     0.294        |     0.857         |
|     Ptx     |     BRCA2 | E135S        |     0.71575659            |     AAAAAAAAAAAAAMAAAAAAAAAAAAA            |     AAAAAAAAAAAAARAAAAAAAAAAAAA                 |     0.173        |     0.556         |

where:
- `patient` is the patient id (**optional**). If this column is not provided, `--patient-id` must be given as input when starting NeoFox (see [here](/03_03_usage.md)). Of note, providing this column allows to put the neoantigen candidates of several patients into one table.
- `gene` is the HGNC gene symbol
- `substitution`  represents a single amino acid substitution with single letter amino acids (eg: I547T). This column allows the detection of INDEL sequences which are removed from the dataset and not processed.  
- `+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)` the neoantigen candidate sequence, i.e. the mutated amino acid sequence. The mutation should be located in the middle, flanked by 13 amino acid on both sites (IUPAC 1 respecting casing, eg: A)
- `[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)` the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A)
- `transcript_expression` the transcript expression. Should be empty if no value available
- `VAF_in_tumor` variant allele frequency in the DNA (**optional**)
- `VAF_in_RNA` variant allele frequency in the RNA (**optional**, this will be estimated using the `VAF_in_tumor` if not available)


### File with patient information

This is an dummy example of a patient file in tabular format:  

| identifier | mhcIAlleles                                                                  | mhcIIAlleles                                                                                                                                                   | isRnaAvailable | tumorType |
|------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|-----------|
| Ptx        | HLA-A\*03:01,HLA-A\*29:02,HLA-B\*07:02,HLA-B\*44:03,HLA-C\*07:02,HLA-C*16:01 | HLA-DRB1\*03:01,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*05:01,HLA-DQB1\*01:01,HLA-DQB1\*04:02,HLA-DPA1\*01:03,HLA-DPA1\*03:01,HLA-DPB1\*13:01,HLA-DPB1*04:02 | TRUE           | HNSC      |
| Pty        | HLA-A\*02:01,HLA-A\*30:01,HLA-B\*07:34,HLA-B\*44:03,HLA-C\*07:02,HLA-C*07:02 | HLA-DRB1\*04:02,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*04:01,HLA-DQB1\*03:02,HLA-DQB1\*14:01,HLA-DPA1\*01:03,HLA-DPA1\*02:01,HLA-DPB1T\*02:01,HLA-DPB1*04:01 | FALSE          | HNSC      |

where:
- `identifier`: the patient identifier
- `mhcIAlleles`: comma separated MHC I alleles of the patient for HLA-A, HLA-B and HLA-C. If homozygous, the allele should be added twice.
- `mhcIIAlleles`: comma separated  MHC II alleles of the patient for HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1 and HLA-DPB1. If homozygous, the allele should be added twice.
- `isRnaAvailable`: whether RNA was available for the analysis. ***If  false, then expression value will be imputed from TCGA gene expression data.*** If true, then the `VAF_in_RNA` field will be used when available, else `VAF_in_DNA` will be used.
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

## JSON format

## Neoantigen candidates in JSON format 

Besides tabular format, neoantigen candidates can be provided as a list of neoantigen models in JSON format as shown below. To simplify, only one full neoantigen model is shown. The terminology follows the descriptions for the [model file](#tabular-format). For a more detailed description of the models, please refer to [here](05_models.md):  

```json
[{
    "identifier": "odJ99FdqvJoK1znK+iCpWQ==",
    "patientIdentifier": "Pt29",
    "gene": "BRCA2",
    "mutation": {
        "wildTypeXmer": "AAAAAALAAAAA",
        "mutatedXmer": "AAAAAAFAAAAA"
    }
}]
``` 
