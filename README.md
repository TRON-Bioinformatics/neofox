# <img src="docs/figures/neofox_logo_small.png" width="10%">  NeoFox - NEOantigen Feature tOolboX   

[![PyPI version](https://badge.fury.io/py/neofox.svg)](https://badge.fury.io/py/neofox)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/neofox/badges/version.svg)](https://anaconda.org/bioconda/neofox)
[![Documentation Status](https://readthedocs.org/projects/neofox/badge/?version=latest)](https://neofox.readthedocs.io/en/latest/?badge=latest)


NeoFox annotates neoantigen candidate sequences with published neoantigen features. 

**For a detailed documentation, please check out [https://neofox.readthedocs.io](https://neofox.readthedocs.io/)**

If you use NeoFox, please cite the following publication:  
Franziska Lang, Pablo Riesgo-Ferreiro, Martin Löwer, Ugur Sahin, Barbara Schrörs, **NeoFox: annotating neoantigen candidates with neoantigen features**, Bioinformatics, Volume 37, Issue 22, 15 November 2021, Pages 4246–4247, https://doi.org/10.1093/bioinformatics/btab344   

## Table of Contents

[1 Implemented neoantigen features](#1-Implemented-Neoantigen-Features)  
[2 NeoFox requirements](#2-NeoFox-Requirements)  
[3 Usage from the command line](#3-Usage-from-the-command-line)  
[4 Input data](#4-input-data)  
[5 Output data](#5-output-data)  

## 1 Implemented Neoantigen Features

NeoFox covers the following neoantigen features and prediction algorithms:  

| Name                                                    | Reference                                                                | DOI                                                                                       |
|---------------------------------------------------------|--------------------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| MHC I binding affinity/rank score (netMHCpan-v4.1)      | Reynisson et al, 2020, Nucleic Acids Research                             | https://doi.org/10.4049/jimmunol.1700893                                                  |
| MHC II binding affinity/rank score (netMHCIIpan-v4.3)   | Nilsson et al, 2023, Science Adv.                                           | https://doi.org/10.1126/sciadv.adj6367                                                         |
| MixMHCpred score v2.2  §                              | Bassani-Sternberg et al., 2017, PLoS Comp Bio; Gfeller, 2018, J Immunol. | https://doi.org/10.1371/journal.pcbi.1005725 ,   https://doi.org/10.4049/jimmunol.1800914 |
| MixMHC2pred score v2.0.2  §                              | Racle et al, 2019, Nat. Biotech. 2019                                    | https://doi.org/10.1038/s41587-019-0289-6                                                 |
| Differential Agretopicity Index (DAI)                   | Duan et al, 2014, JEM; Ghorani et al., 2018, Ann Oncol.                  | https://doi.org/10.1084/jem.20141308                                                      |
| Self-Similarity                                         | Bjerregaard et al, 2017, Front Immunol.                                  | https://doi.org/10.3389/fimmu.2017.01566                                                  |
| IEDB immunogenicity                                     | Calis et al, 2013, PLoS Comput Biol.                                     | https://doi.org/10.1371/journal.pcbi.1003266                                              |
| Neoantigen dissimilarity                                | Richman et al, 2019, Cell Systems                                        | https://doi.org/10.1016/j.cels.2019.08.009                                                |
| PHBR-I  §                                                | Marty et al, 2017, Cell                                                  | https://doi.org/10.1016/j.cell.2017.09.050                                                |
| PHBR-II  §                                               | Marty Pyke et al, 2018, Cell                                             | https://doi.org/10.1016/j.cell.2018.08.048                                                |
| Generator rate                                          | Rech et al, 2018, Cancer Immunology Research                             | https://doi.org/10.1158/2326-6066.CIR-17-0559                                             |
| Recognition potential   §                                | Łuksza et al, 2017, Nature; Balachandran et al, 2017, Nature             | https://doi.org/10.1038/nature24473 , https://doi.org/10.1038/nature24462                 |
| Vaxrank                                                 | Rubinsteyn, 2017, Front Immunol                                          | https://doi.org/10.3389/fimmu.2017.01807                                                  |
| Priority score                                          | Bjerregaard et al, 2017, Cancer Immunol Immunother.                      | https://doi.org/10.1007/s00262-017-2001-3                                                 |
| PRIME  §                                                 | Schmidt et al., 2021, Cell Reports Medicine                            | https://doi.org/10.1016/j.xcrm.2021.100194                                             |
| HEX §                                                   | Chiaro et al., 2021, Cancer Immunology Research                            | https://doi.org/10.1158/2326-6066.CIR-20-0814                                             |

*§ currently not supported for mouse*

## 2 NeoFox Requirements
 
NeoFox depends on the following tools:  

- Python 3.11
- BLAST 2.10.1
- netMHCpan 4.1
- netMHCIIpan 4.3
- MixMHCpred 2.2 (optional)
- MixMHC2pred 2.0.2 (optional)
- PRIME 2.0 (optional)

Install from PyPI:
```
pip install neofox
```

Or install from bioconda:
```
conda install bioconda::neofox
```


## 3 Usage from the command line

NeoFox can be used from the command line as shown below or programmatically (see [https://neofox.readthedocs.io](https://neofox.readthedocs.io/) for more information).

````commandline
neofox --input-file neoantigens_candidates.tsv \
    --patient-data patient_data.txt \
    --output-folder /path/to/out \
    [--output-prefix out_prefix]  \
    [--organism human|mouse]  \
    [--rank-mhci-threshold 2.0] \
    [--rank-mhcii-threshold 5.0] \
    [--num-cpus] \
    [--config] \
    [--patient-id] \
    [--with-all-neoepitopes] \
    [--verbose]
````
- `--input-file`: tab-separated values table with neoantigen candidates represented by long mutated peptide sequences 
 as described [here](03_01_input_data.md#tabular-file-format) (extensions .txt and .tsv) or JSON file neoantigens in 
 NeoFox model format as  described [here](03_01_input_data.md#json-file-format) (extension .json)
- `--patient-data`: a table of tab separated values containing metadata on the patient as  described [here](03_01_input_data.md#file-with-patient-information)
- `--output-folder`: path to the folder to which the output files should be written 
- `--output-prefix`: prefix for the output files (*optional*)
- `--with-all-neoepitopes`: output annotations for all MHC-I and MHC-II neoepitopes on all HLA alleles (*optional*)
- `--rank-mhci-threshold`: MHC-I epitopes with a netMHCpan predicted rank greater than or equal than this threshold will be filtered out (*optional*)
- `--rank-mhcii-threshold`: MHC-II epitopes with a netMHCIIpan predicted rank greater than or equal than this threshold will be filtered out (*optional*)
- `--organism`: the organism to which the data corresponds. Possible values: [human, mouse]. Default value: human
- `--num-cpus`: number of CPUs to use (*optional*)
- `--config`: a config file with the paths to dependencies as shown below  (*optional*)
- `--patient-id`: patient identifier (*optional*, this is only relevant if the column `patientIdentifier` is missing in the candidate input file)
- `--verbose`: get detailed logs
                        
                        
The optional config file with the paths to the dependencies can look like this:  
````commandline
NEOFOX_REFERENCE_FOLDER=path/to/reference/folder
NEOFOX_BLASTP=path/to/ncbi-blast-2.10.1+/bin/blastp
NEOFOX_NETMHCPAN=path/to/netMHCpan-4.1/netMHCpan
NEOFOX_NETMHC2PAN=path/to/netMHCIIpan-4.3/netMHCIIpan
NEOFOX_MIXMHCPRED=path/to/MixMHCpred-2.2/MixMHCpred
NEOFOX_MIXMHC2PRED=path/to/MixMHC2pred-2.0.1/MixMHC2pred_unix
NEOFOX_MAKEBLASTDB=path/to/ncbi-blast-2.8.1+/bin/makeblastdb
NEOFOX_PRIME=/path/to/PRIME-2.0/PRIME
````

## 4 Input data

### 4.1 Neoantigen candidates in tabular format
This is an dummy example of a table with neoantigen candidates:  

| gene  | wildTypeXmer       | mutatedXmer        | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | dnaVariantAlleleFrequency | external_annotation_1 | external_annotation_2 |
|-------|-----------------------------|-----------------------------|-------------------|---------------|---------------------------|---------------------------|-----------------------|-----------------------|
| BRCA2 | AAAAAAAAAAAAALAAAAAAAAAAAAA | AAAAAAAAAAAAAFAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAMAAAAAAAAAAAAA | AAAAAAAAAAAAARAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAGAAAAAAAAAAAAA | AAAAAAAAAAAAAKAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAACAAAAAAAAAAAAA | AAAAAAAAAAAAAEAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |
| BRCA2 | AAAAAAAAAAAAAKAAAAAAAAAAAAA | AAAAAAAAAAAAACAAAAAAAAAAAAA | Ptx               | 7.942         | 0.85                      | 0.34                      | some_value            | some_value            |

where:
- `gene`: the HGNC gene symbol   
- `mutatedXmer`: the neoantigen candidate sequence, i.e. the mutated amino acid sequence. The mutation should be located in the middle, flanked by 13 amino acid on both sites (IUPAC 1 respecting casing, eg: A)
- `wildTypeXmer`: the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A)
- `patientIdentifier`: the patient identifier
- `rnaExpression`: RNA expression. (**optional**) (see *NOTE*) This value can be in any format chosen by the user (e.g. TPM, RPKM) but it is recommended to be consistent for data that should be compared.
- `rnaVariantAlleleFrequency`: the variant allele frequency calculated from the RNA (**optional**, this will be estimated using the `dnaVariantAlleleFrequency` if not available)
- `dnaVariantAlleleFrequency`: the variant allele frequency calculated from the DNA (**optional**)  

**NOTE:** If rnaExpression is not provided, expression will be estimated by gene expression in a respective TCGA cohort and this value will be used for relevant features. The TCGA cohort to be used for imputation of gene expression needs to be indicated in the `tumorType` in the patient data (see below). If `tumorType` is not provided, expression will not be imputed.  

### 4.2 Neoantigen candidates in JSON format 

Besides tabular format, neoantigen candidates can be provided as a list of neoantigen models in JSON format as shown below. To simplify, only one full neoantigen model is shown:  

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

### 4.3 Patient-data format

This is an dummy example of a patient file:  

| identifier | mhcIAlleles                                                                  | mhcIIAlleles                                                                                                                                                   | tumorType |
|------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
| Ptx        | HLA-A\*03:01,HLA-A\*29:02,HLA-B\*07:02,HLA-B\*44:03,HLA-C\*07:02,HLA-C*16:01 | HLA-DRB1\*03:01,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*05:01,HLA-DQB1\*01:01,HLA-DQB1\*04:02,HLA-DPA1\*01:03,HLA-DPA1\*03:01,HLA-DPB1\*13:01,HLA-DPB1*04:02 | HNSC      |
| Pty        | HLA-A\*02:01,HLA-A\*30:01,HLA-B\*07:34,HLA-B\*44:03,HLA-C\*07:02,HLA-C*07:02 | HLA-DRB1\*04:02,HLA-DRB1\*08:01,HLA-DQA1\*03:01,HLA-DQA1\*04:01,HLA-DQB1\*03:02,HLA-DQB1\*14:01,HLA-DPA1\*01:03,HLA-DPA1\*02:01,HLA-DPB1\*02:01,HLA-DPB1*04:01 | HNSC      |

where:
- `identifier`: the patient identifier
- `mhcIAlleles`: comma separated MHC I alleles of the patient for HLA-A, HLA-B and HLA-C. If homozygous, the allele should be added twice.
- `mhcIIAlleles`: comma separated  MHC II alleles of the patient for HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1 and HLA-DPB1. If homozygous, the allele should be added twice.
- `tumorType`: tumour entity in TCGA study abbreviation format (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). This field is required for expression imputation. The supported tumor types are listed under "Input data" in the [documentation](https://neofox.readthedocs.io/en/latest/03_01_input_data.html).


## 5 Output data

The output data is returned by default in tsv and json format. With the command line flag `--with-all-neoepitopes`, two additional files are generated containing the epitope candidates for MHCI and MHCII with NetMHCpan predictions below the given thresholds.
