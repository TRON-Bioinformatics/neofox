# NeoFox - NEOantigen Feature tOolboX

[![DOI](https://zenodo.org/badge/294667387.svg)](https://zenodo.org/badge/latestdoi/294667387)
[![PyPI version](https://badge.fury.io/py/neofox.svg)](https://badge.fury.io/py/neofox)
[![Downloads](https://crate.io/packages/neofox)](https://crate.io/packages/neofox)


NeoFox annotates neoantigen candidate sequences with published neo-epitope descriptors. For a detailed documentation, please check out [https://neofox.readthedocs.io](https://neofox.readthedocs.io/)

| Name                                                    | Reference                                                                | DOI                                                                                       |
|---------------------------------------------------------|--------------------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| MHC I binding affinity/rank score (netMHCpan-v4.0)      | Jurtz et al, 2017, The Journal of Immunology                             | https://doi.org/10.4049/jimmunol.1700893                                                  |
| MHC II binding affinity/rank score (netMHCIIpan-v3.2)   | Jensen et al, 2018, Immunology                                           | https://doi.org/10.1111/imm.12889                                                         |
| MixMHCpred score v2.1                                   | Bassani-Sternberg et al., 2017, PLoS Comp Bio; Gfeller, 2018, J Immunol. | https://doi.org/10.1371/journal.pcbi.1005725 ,   https://doi.org/10.4049/jimmunol.1800914 |
| MixMHC2pred score v1.2                                  | Racle et al, 2019, Nat. Biotech. 2019                                    | https://doi.org/10.1038/s41587-019-0289-6                                                 |
| Differential Agretopicity Index (DAI)                   | Duan et al, 2014, JEM; Ghorani et al., 2018, Ann Oncol.                  | https://doi.org/10.1084/jem.20141308                                                      |
| Self-Similarity                                         | Bjerregaard et al, 2017, Front Immunol.                                  | https://doi.org/10.3389/fimmu.2017.01566                                                  |
| IEDB immunogenicity                                     | Calis et al, 2013, PLoS Comput Biol.                                     | https://doi.org/10.1371/journal.pcbi.1003266                                              |
| Neoantigen dissimilarity                                | Richman et al, 2019, Cell Systems                                        | https://doi.org/10.1016/j.cels.2019.08.009                                                |
| PHBR-I                                                  | Marty et al, 2017, Cell                                                  | https://doi.org/10.1016/j.cell.2017.09.050                                                |
| PHBR-II                                                 | Marty Pyke et al, 2018, Cell                                             | https://doi.org/10.1016/j.cell.2018.08.048                                                |
| Generator rate                                          | Rech et al, 2018, Cancer Immunology Research                             | https://doi.org/10.1158/2326-6066.CIR-17-0559                                             |
| Recognition potential                                   | Łuksza et al, 2017, Nature; Balachandran et al, 2017, Nature             | https://doi.org/10.1038/nature24473 , https://doi.org/10.1038/nature24462                 |
| Vaxrank                                                 | Rubinsteyn, 2017, Front Immunol                                          | https://doi.org/10.3389/fimmu.2017.01807                                                  |
| Priority score                                          | Bjerregaard et al, 2017, Cancer Immunol Immunother.                      | https://doi.org/10.1007/s00262-017-2001-3                                                 |
| Tcell predictor                                         | Besser et al, 2019, Journal for ImmunoTherapy of Cancer                  | https://doi.org/10.1186/s40425-019-0595-z                                                 |
| neoag                                                   | Smith et al, 2019, Cancer Immunology Research                            | https://doi.org/10.1158/2326-6066.CIR-19-0155                                             |

## NeoFox Requirements
 
**Required Software/Tools/Dependencies:**  
- Python 3.7
- R 3.6.0
- BLAST 2.10.1
- netMHCpan 4.0
- netMHCIIpan 3.2
- MixMHCpred 2.1
- MixMHC2pred 1.2


## Usage from the command line

````commandline
neofox --model-file/--candidate-file/--json-file neoantigens_candidates.tab/neoantigens_candidates.json --patient-data patient_data.txt --output-folder /path/to/out --output-prefix out_prefix [--with-short-wide-table] [--with-tall-skinny-table] [--with-json] [--num_cpus]
````

where:
- `--candidate-file`: tab-separated values table with neoantigen candidates represented by long mutated peptide sequences
- `--model-file`: tab-separated values table with neoantigens in NeoFox model format
- `--json-file`: JSON file neoantigens in NeoFox model format
- `--patient-id`: patient identifier (*optional*, this will be used if the patient id the column `patient` is missing the candidate input file)
- `--patient-data`: a table of tab separated values containing metadata on the patient
- `--output-folder`: path to the folder to which the output files should be written 
- `--output-prefix`: prefix for the output files (*optional*)
- `--with-short-wide-table`: output file in short-wide format (*optional*)
- `--with-tall-skinny-table`: output file in tall-skinny format (*optional*)
- `--with-json`: output file in JSON format (*optional*)
- `--num_cpus`: number of CPUs to use (*optional*)

### Input data

#### model-file format  

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

#### candidate-file format  

Alternatively, neoantigen candidates can be provided in `candidate-file` format. In principle the columns are the same as in the `model-file`. Of note, `candidate-file` allows for an optional patient id in the data table. This is an dummy example:  

|     patient |     gene  | substitution |     transcript_expression |     +-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) |     [WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) |     VAF_in_tumor |     VAF_in_RNA    |
|-------------|-----------|--------------|---------------------------|--------------------------------------------|-------------------------------------------------|------------------|-------------------|
|     Ptx     |     BRCA2 | I547T        |     0.51950689            |     AAAAAAAAAAAAAFAAAAAAAAAAAAA            |     AAAAAAAAAAAAALAAAAAAAAAAAAA                 |     0.294        |     0.857         |
|     Ptx     |     BRCA2 | E135S        |     0.71575659            |     AAAAAAAAAAAAAMAAAAAAAAAAAAA            |     AAAAAAAAAAAAARAAAAAAAAAAAAA                 |     0.173        |     0.556         |

where:
- `patient` is the patient id (**optional**). If this column is not provided, `--patient-id` must be given as input when starting NeoFox. Of note, providing this column allows to put the neoantigen candidates of several patients into one table.
- `gene` is the HGNC gene symbol
- `substitution`  represents a single amino acid substitution with single letter amino acids (eg: I547T). This column allows the detection of INDEL sequences which are removed from the dataset and not processed.  
- `+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)` the neoantigen candidate sequence, i.e. the mutated amino acid sequence. The mutation should be located in the middle, flanked by 13 amino acid on both sites (IUPAC 1 respecting casing, eg: A)
- `[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)` the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A)
- `transcript_expression` the transcript expression. Should be empty if no value available
- `VAF_in_tumor` variant allele frequency in the DNA (**optional**)
- `VAF_in_RNA` variant allele frequency in the RNA (**optional**, this will be estimated using the `VAF_in_tumor` if not available)

### JSON format 

Besides tabular format, neoantigen candidates can be provided as a list of neoantigen models in JSON format as shown below. To simplify, only one full neoantigen model is shown:  

```json
[{
    "identifier": "odJ99FdqvJoK1znK+iCpWQ==",
    "patientIdentifier": "Ptx",
    "gene": "BRCA2",
    "mutation": {
        "wildTypeXmer": "AAAAAAAAAAAAALAAAAAAAAAAAAA",
        "mutatedXmer": "AAAAAAAAAAAAAFAAAAAAAAAAAAA"
    }
}]
``` 

#### patient-file format

This is an dummy example of a patient file:  

| identifier | mhcIAlleles                                                                  | mhcIIAlleles                                                                                                                                                   | isRnaAvailable | tumorType |
|------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|-----------|
| Ptx        | HLA-A\*03:01,HLA-A\*29:02,HLA-B\*07:02,HLA-B\*44:03,HLA-C\*07:02,HLA-C*16:01 | HLA-DRB1\*03:01,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*05:01,HLA-DQB1\*01:01,HLA-DQB1\*04:02,HLA-DPA1\*01:03,HLA-DPA1\*03:01,HLA-DPB1\*13:01,HLA-DPB1*04:02 | TRUE           | HNSC      |
| Pty        | HLA-A\*02:01,HLA-A\*30:01,HLA-B\*07:34,HLA-B\*44:03,HLA-C\*07:02,HLA-C*07:02 | HLA-DRB1\*04:02,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*04:01,HLA-DQB1\*03:02,HLA-DQB1\*14:01,HLA-DPA1\*01:03,HLA-DPA1\*02:01,HLA-DPB1\*02:01,HLA-DPB1*04:01 | FALSE          | HNSC      |

where:
- `identifier`: the patient identifier
- `mhcIAlleles`: comma separated MHC I alleles of the patient for HLA-A, HLA-B and HLA-C. If homozygous, the allele should be added twice.
- `mhcIIAlleles`: comma separated  MHC II alleles of the patient for HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1 and HLA-DPB1. If homozygous, the allele should be added twice.
- `isRnaAvailable`: whether RNA was available for the analysis. ***If  false, then expression value will be imputed from TCGA gene expression data.*** If true, then the `VAF_in_RNA` field will be used when available, else `VAF_in_DNA` will be used.
- `tumorType`: tumour entity in TCGA study abbreviation format (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). This field is required for expression imputation and at the moment the following tumor types are supported:

### Output data

The output data is returned in a short wide tab separated values file (`--with-short-wide-table`). Optionally, it can be provided in a tall skinny tab separated values file (`--with-tall-skinny-table`) or in JSON (`--with-json`).  

For a more information, please check out our documentation on [https://neofox.readthedocs.io](https://neofox.readthedocs.io/)