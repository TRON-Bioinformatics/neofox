# NeoFox - NEOantigen Feature tOolboX

[![DOI](https://zenodo.org/badge/294667387.svg)](https://zenodo.org/badge/latestdoi/294667387)
[![PyPI version](https://badge.fury.io/py/neofox.svg)](https://badge.fury.io/py/neofox)


NeoFox annotates neoantigen candidate sequences with published neo-epitope descriptors. For a detailed documentation, please check out [https://neofox.readthedocs.io](https://neofox.readthedocs.io/)

**Published Descriptors:**
- netMHCpan *(Jurtz et al, 2017, The Journal of Immunology )*  
- netMHCIIpan *(Jensen et al, 2018, Immunology )*  
- IEDB immunogenicity *(Calis et al, 2013, PLoS Comput Biol.)*  
- Self-similarity, Conserved vs. Improved Binding  *(Bjerregaard et al, 2017, Front Immunol.)*  
- Priority Score *(Bjerregaard et al, 2017, Cancer Immunol Immunother.)*  
- DAI *(Duan et al., 2014, JEM; Ghorani et al., 2018, Ann Oncol.)*  
- Neoantigen Fitness *(Luksza et al., 2017, Nature; Balachandran et al, 2017, Nature)*  
- Residue-centric presentation score (best_rank) & Patient harmonic Best Rank (PHBR-I/II) *(Marty et al, 2017, Cell; Marty et al, 2018, Cell)*  
- Classically vs Alternatively Defined Neopitopes & Generator Rate *(Rech et al., 2018, Cancer Immunology Research)*  
- Tcell_predictor *(Besser et al, 2019, Journal for ImmunoTherapy of Cancer)*  
- neoag *(Smith et al, 2019, Cancer Immunology Research)*
- neoantigen dissimilarity *(Richman et al, 2019, Cell Systems)*
- MixMHCpred *(Bassani-Sternberg et al., 2017, PLoS Comp Bio; Gfeller, 2018; J Immunol.)*
- MixMHC2pred *(Racle et al, 2019, Nat. Biotech. 2019)*
- Vaxrank *(Rubinsteyn, 2017, Front Immunol;Wang, 2019, Bioinformatics)*

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
neofox --model-file/--candidate-file/--json-file neoantigens_candidates.tab/neoantigens_candidates.json --patient-id Ptx --patient-data/--patient-data-json patient_data.txt/patient_data.json --output-folder /path/to/out --output-prefix out_prefix [--with-short-wide-table] [--with-tall-skinny-table] [--with-json] [--num_cpus]
````

where:
- `--candidate-file`: tab-separated values table with neoantigen candidates represented by long mutated peptide sequences as described [here](03_01_input_data.md#tabular-format)
- `--model-file`: tab-separated values table with neoantigens in NeoFox model format as described [here](03_01_input_data.md#tabular-format)
- `--json-file`: JSON file neoantigens in NeoFox model format as  described [here](03_01_input_data.md#json-format)
- `--patient-id`: patient identifier (*optional*, this will be used if the patient id the column `patient` is missing the candidate input file)
- `--patient-data`: a table of tab separated values containing metadata on the patient as  described [here](03_01_input_data.md#file-with-patient-information)
- `--patient-data-json`: a table patient models as described [here](03_01_input_data.md#patient-file-in-json-format)
- `--output-folder`: path to the folder to which the output files should be written 
- `--output-prefix`: prefix for the output files (*optional*)
- `--with-short-wide-table`: output file in [short-wide](03_02_output_data.md#short-wide-format) format (*optional*)
- `--with-tall-skinny-table`: output file in [tall-skinny](03_02_output_data.md#tall-skinny-format) format (*optional*)
- `--with-json`: output file in [JSON](03_02_output_data.md#json-format) format (*optional*)
- `--num_cpus`: number of CPUs to use (*optional*)


### Input data

1. **model-file format**  

This is an dummy example of a table with neoantigen candidates in `model-file` format:  

| transcript.assembly | transcript.gene | transcript.identifier | mutation.mutatedAminoacid | mutation.position | mutation.wildTypeAminoacid | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | dnaVariantAlleleFrequency |
|---------------------|-----------------|-----------------------|---------------------------|-------------------|----------------------------|-------------------|---------------|---------------------------|---------------------------|
| hg19                | BRCA2           | uc003kii.3            | L                         | 935               | F                          | Ptx               | 4.512         | 0.4675                    | 0.36103                   |
| hg19                | BRCA2           | uc003kii.3            | M                         | 518               | R                          | Ptx               | 0.154         | 0.015404                  | 0.034404                  |
| hg19                | BRCA2           | uc003kii.3            | G                         | 285               | K                          | Ptx               | 8.841207      | 0.89387                   | 0.51924                   |

where:
- `transcript.assembly`: the assembly of the reference genome (only hg19 is supported)
- `transcript.gene`: the HGNC gene symbol   
- `transcript.identifier`: a transcript identifier
- `mutation.mutatedXmer`: the neoantigen candidate sequence, i.e. the mutated amino acid sequence. The mutation should be located in the middle, flanked by 13 amino acid on both sites (IUPAC 1 respecting casing, eg: A)
- `mutation.wildTypeXmer`: the equivalent non-mutated amino acid sequence (IUPAC 1 respecting casing, eg: A)
- `mutation.mutatedAminoacid`: the mutated amino acid (IUPAC 1 or 3 letters respecting casing, eg: A and Ala)
- `mutation.position`: the 1 based position of the mutation in the protein
- `mutation.wildTypeAminoacid`: the wild type amino acid (IUPAC 1 or 3 letters respecting casing, eg: A and Ala)
- `patientIdentifier`: the patient identifier
- `rnaExpression`: the transcript expression. Should be empty if no value available
- `rnaVariantAlleleFrequency`: the variant allele frequency calculated from the RNA (**optional**, this will be estimated using the `dnaVariantAlleleFrequency` if not available)
- `dnaVariantAlleleFrequency`: the variant allele frequency calculated from the DNA (**optional**)  



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

### Output data

The output data is returned in a short wide tab separated values file (`--with-short-wide-table`). Optionally, it can be provided in a tall skinny tab separated values file (`--with-tall-skinny-table`) or in JSON (`--with-json`).
