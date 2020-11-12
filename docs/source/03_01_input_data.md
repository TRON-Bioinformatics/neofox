# Input data

## General information
NeoFox requires two input files: a file with neoantigen candidates that derive from point mutation and a file with patient data. These files can be provided either in tabular format when using the the command lined or in JSON format when using the API.  
Of note, the neoantigen candidate file may contain additional user-specific input that will be kept during the annotation process.

## Tabular format
### File with neonatigen candidates
We allow two different tabular formats of the neoantigen candidate file: in `model-file` or `candidate-file`format. The neoantigen candidate input file must follow one of these designs. However, additional user-specific columns are allowed and will be kept during the annotation process.   

1. **model-file format**  
   This is an dummy example of a table with neoantigen candidates in `model-file` format:  
   
   | transcript.assembly | transcript.gene | transcript.identifier | mutation.leftFlankingRegion | mutation.mutatedAminoacid | mutation.position | mutation.rightFlankingRegion | mutation.wildTypeAminoacid | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | dnaVariantAlleleFrequency |
   |---------------------|-----------------|-----------------------|-----------------------------|---------------------------|-------------------|------------------------------|----------------------------|-------------------|---------------|---------------------------|---------------------------|
   | hg19                | BRCA2           | uc003kii.3            | AAAAAA                      | L                         | 935               | AAAAA                        | F                          | Ptx               | 4.512         | 0.4675                    | 0.36103                   |
   | hg19                | BRCA2           | uc003kii.3            | AAAAAA                      | M                         | 518               | AAAAA                        | R                          | Ptx               | 0.154         | 0.015404                  | 0.034404                  |
   | hg19                | BRCA2           | uc003kii.3            | AAAAAA                      | G                         | 285               | AAAAA                        | K                          | Ptx               | 8.841207      | 0.89387                   | 0.51924                   |
   
   where:
   - `transcript.assembly` - the assembly of the reference genome (only hg19 is supported)
   - `transcript.gene` - the HGMC gene symbol   
   - `transcript.identifier` - a transcript identifier
   - `mutation.leftFlankingRegion` - the amino acids flanking the mutation on the left (in IUPAC one letter symbols)
   - `mutation.mutatedAminoacid` - the mutated amino acid (IUPAC 1 or 3 letters respecting casing, eg: A and Ala)
   - `mutation.position` - the 1 based position of the mutation in the protein
   - `mutation.rightFlankingRegion` - the amino acids flanking the mutation on the right (in IUPAC one letter symbols)
   - `mutation.wildTypeAminoacid` - the wild type amino acid (IUPAC 1 or 3 letters respecting casing, eg: A and Ala)
   - `patientIdentifier` - the patient identifier
   - `rnaExpression` - the transcript expression. Should be empty if no value available
   - `rnaVariantAlleleFrequency` - the variant allele frequency calculated from the RNA (**optional**, this will be estimated using the `dnaVariantAlleleFrequency` if not available)
   - `dnaVariantAlleleFrequency` - the variant allele frequency calculated from the DNA (**optional**)  



2. **candidate-file format**  
   Alternatively, neoantigen candidates can be provided in `candidate-file` format. This is an dummy example:  
   
   |  patient | gene   | UCSC_transcript | transcript_expression | substitution | +-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) | [WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) | VAF_in_tumor | VAF_in_RNA |
   |----------|--------|-----------------|-----------------------|--------------|----------------------------------------|---------------------------------------------|--------------|------------|
   |  Ptx     | VCAN   | uc003kii.3      | 0.519506894           | I547T        | DEVLGEPSQDILVTDQTRLEATISPET            | DEVLGEPSQDILVIDQTRLEATISPET                 |  0.294       |  0.857     |
   |  Ptx     | TRIM25 | uc001zii.3      | 0.715756594           | E135S        | PQLHKNTVLCNVVSQFLQADLAREPPA            | PQLHKNTVLCNVVEQFLQADLAREPPA                 |  0.173       |  0.556     |

   where:
   - `patient` is the patient id (**optional**). If this column is not provided, `--patient-id` must be given as input when starting NeoFox (see [here](/03_03_usage.md)). Of note, providing this column allows to put the neoantigen candidates of several patients into one table.
   - `gene` is the HGNC gene symbol
   - `UCSC_trancript` is the UCSC transcript id including the version. (The user can enter a non-UCSC transcript id, if no UCSC transcript id is available)
   - `substitution` represents a single amino acid substitution with single letter amino acids (eg: I547T)
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
| Pty        | HLA-A\*02:01,HLA-A\*30:01,HLA-B\*07:34,HLA-B\*44:03,HLA-C\*07:02,HLA-C*07:02 | HLA-DRB1\*04:02,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*04:01,HLA-DQB1\*03:02,HLA-DQB1\*14:01,HLA-DPA1\*01:03,HLA-DPA1\*02:01,HLA-DPB1\*02:01,HLA-DPB1*04:01 | FALSE          | HNSC      |

where:
- `identifier` - the patient identifier
- `mhcIAlleles` - the list of MHC I alleles of the patient for HLA-A, HLA-B and HLA-C. If homozygous, an allele should be added twice.
- `mhcIIAlleles` - the list of MHC II alleles of the patient for HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1 and HLA-DPB1. If homozygous, an allele should be added twice.
- `isRnaAvailable` - whether RNA was available for the analysis. ***If  false, then expression value will be imputed from TCGA gene expression data.*** If true, then the `VAF_in_RNA` field will be used, else `VAF_in_DNA` will be used.
- `tumorType` - tumour entity in TCGA study abbreviation format (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). This field is required for expression imputation and at the moment the following tumor types are supported:

   | Study Name                                                         | Sudy Abbreviation |
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

````json
[{
  "identifier": "jETwpX0R9iEiQz2SMpHkPQ==",
  "patientIdentifier": "P123",
  "transcript": {
    "identifier": "uc003kii.3",
    "assembly": "hg19",
    "gene": "VCAN"
  },
  "mutation": {
    "position": 1007,
    "wildTypeXmer": "DEVLGEPSQDILVIDQTRLEATISPET",
    "wildTypeAminoacid": "I",
    "mutatedXmer": "DEVLGEPSQDILVTDQTRLEATISPET",
    "mutatedAminoacid": "T",
    "leftFlankingRegion": "DEVLGEPSQDILV",
    "sizeLeftFlankingRegion": 13,
    "rightFlankingRegion": "DQTRLEATISPET",
    "sizeRightFlankingRegion": 13
  },
  "rnaExpression": 0.519506894,
  "dnaVariantAlleleFrequency": 0.294573643,
  "rnaVariantAlleleFrequency": 0.857142857
}, {
    next_neoantigen_candidate
}]
```` 

## Patient file in JSON format 
This is an example of a patient file with patient. Note that the patient models contain models of the MHC I and MHC II alleles. These models are shown in more details [here](05_models.md) in JSON format. To simplify, only one full patient model is shown: 

````json
[{
       "identifier": "P123",
       "isRnaAvailable": true,
       "tumor_type": true,
       "mhc1": [
          {
             "zygosity": "HETEROZYGOUS",
             "alleles": [
                {
                   "fullName": "HLA-A*01:01:02:03N",
                   "name": "HLA-A*01:01",
                   "gene": "A",
                   "group": "01",
                   "protein": "01"
                },
                {
                   "fullName": "HLA-A*01:02:02:03N",
                   "name": "HLA-A*01:02",
                   "gene": "A",
                   "group": "01",
                   "protein": "02"
                }
             ]
          },
          {
             "name": "B",
             "alleles": [
                {
                   "fullName": "HLA-B*01:01:02:04N",
                   "name": "HLA-B*01:01",
                   "gene": "B",
                   "group": "01",
                   "protein": "01"
                }
             ]
          },
          {
             "name": "C",
             "zygosity": "HEMIZYGOUS",
             "alleles": [
                {
                   "fullName": "HLA-C*01:01",
                   "name": "HLA-C*01:01",
                   "gene": "C",
                   "group": "01",
                   "protein": "01"
                }
             ]
          }
       ],
       "mhc2": [
          {
             "genes": [
                {
                   "alleles": [
                      {
                         "fullName": "HLA-DRB1*01:01",
                         "name": "HLA-DRB1*01:01",
                         "gene": "DRB1",
                         "group": "01",
                         "protein": "01"
                      }
                   ]
                }
             ],
             "isoforms": [
                {
                   "name": "HLA-DRB1*01:01",
                   "betaChain": {
                      "fullName": "HLA-DRB1*01:01",
                      "name": "HLA-DRB1*01:01",
                      "gene": "DRB1",
                      "group": "01",
                      "protein": "01"
                   }
                }
             ]
          },
          {
             "name": "DP",
             "genes": [
                {
                   "name": "DPA1",
                   "zygosity": "HETEROZYGOUS",
                   "alleles": [
                      {
                         "fullName": "HLA-DPA1*01:01",
                         "name": "HLA-DPA1*01:01",
                         "gene": "DPA1",
                         "group": "01",
                         "protein": "01"
                      },
                      {
                         "fullName": "HLA-DPA1*01:02",
                         "name": "HLA-DPA1*01:02",
                         "gene": "DPA1",
                         "group": "01",
                         "protein": "02"
                      }
                   ]
                },
                {
                   "name": "DPB1",
                   "alleles": [
                      {
                         "fullName": "HLA-DPB1*01:01",
                         "name": "HLA-DPB1*01:01",
                         "gene": "DPB1",
                         "group": "01",
                         "protein": "01"
                      }
                   ]
                }
             ],
             "isoforms": [
                {
                   "name": "HLA-DPA1*01:01-DPB1*01:01",
                   "alphaChain": {
                      "fullName": "HLA-DPA1*01:01",
                      "name": "HLA-DPA1*01:01",
                      "gene": "DPA1",
                      "group": "01",
                      "protein": "01"
                   },
                   "betaChain": {
                      "fullName": "HLA-DPB1*01:01",
                      "name": "HLA-DPB1*01:01",
                      "gene": "DPB1",
                      "group": "01",
                      "protein": "01"
                   }
                },
                {
                   "name": "HLA-DPA1*01:02-DPB1*01:01",
                   "alphaChain": {
                      "fullName": "HLA-DPA1*01:02",
                      "name": "HLA-DPA1*01:02",
                      "gene": "DPA1",
                      "group": "01",
                      "protein": "02"
                   },
                   "betaChain": {
                      "fullName": "HLA-DPB1*01:01",
                      "name": "HLA-DPB1*01:01",
                      "gene": "DPB1",
                      "group": "01",
                      "protein": "01"
                   }
                }
             ]
          },
          {
             "name": "DQ",
             "genes": [
                {
                   "name": "DQA1",
                   "zygosity": "LOSS"
                },
                {
                   "name": "DQB1",
                   "zygosity": "LOSS"
                }
             ]
          }
       ]
    }, {
    next_patient
}]
````

