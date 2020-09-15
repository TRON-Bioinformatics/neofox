# NeoFox - NEOantigen Feature tOolboX


Annotation of mutated peptide sequences (mps) with published neo-epitope descriptors

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
- BLAST 2.8.1
- netMHCpan 4.0
- netMHCIIpan 3.2
- MixMHCpred 2.0.2
- MixMHC2pred 1.1

## Usage

```
neofox --model-file/--icam-file neoantigens.txt --patient-id Ptx --patient-data patient_data.txt --output-folder /path/to/out --output-prefix out_prefix [--with-short-wide-table] [--with-tall-skinny-table] [--with-json] [--num_cpus]
```

### Configuration

Setup the following environment variables pointing to the required software
```
export NEOFOX_BLASTP=/code/ncbi-blast/2.8.1+/bin/blastp
export NEOFOX_MIXMHC2PRED=/code/net/MixMHC2pred/1.1/MixMHC2pred
export NEOFOX_MIXMHCPRED=/code/MixMHCpred/2.0.2/MixMHCpred
export NEOFOX_RSCRIPT=/code/R/3.6.0/bin/Rscript
export NEOFOX_NETMHC2PAN=/code/net/MHCIIpan/3.2/netMHCIIpan
export NEOFOX_NETMHCPAN=/code/net/MHCpan/4.0/netMHCpan
```

Build the reference data as described (here)[neofox/references/build_references.sh] and set the environment variable `NEOFOX_REFERENCE_FOLDER` pointing to the reference folder.

### Input data

- `--icam-file`: tab-separated values table with neoantigens in iCaM output format
- `--model-file`: tab-separated values table with neoantigens in Neofox model format described in [protobuf model](neofox/model/neoantigen.proto)
- `--patient-id`: patient identifier (**optional**, this will be used as the patient id for neoantigens without patient)
- `--patient data`: a table of tab separated values containing metadata on the patient

**NOTE**: provide either `--icam-file` or `--model-file`

Example of iCaM neoantigens table:
```
gene	UCSC_transcript	transcript_expression	substitution	+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)	[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)	VAF_in_tumor	VAF_in_RNA
VCAN	uc003kii.3	0.519506894	I547T	DEVLGEPSQDILVTDQTRLEATISPET	DEVLGEPSQDILVIDQTRLEATISPET 0.294573643	0.857142857
```
where:
- `gene` is the HGNC gene symbol
- `UCSC_trancript` is the UCSC transcript id including the version
- `substitution` represents a single aminoacid substitution with single letter aminoacids (eg: I547T)
- `+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)` the mutated neoantigen
- `[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)` the equivalent aminoacid sequence in the normal tissue
- `transcript_expression` the transcript expression (**optional**)
- `VAF_in_tumor` variant allele frequency in the DNA (**optional**)
- `VAF_in_RNA` variant allele frequency in the RNA (**optional**, this will be estimated using the `VAF_in_tumor` if not available)

Example of Neofox model neoantigens table:
```
gene.assembly	gene.gene	gene.transcriptIdentifier	mutation.leftFlankingRegion	mutation.mutatedAminoacid	mutation.position	mutation.rightFlankingRegion	mutation.wildTypeAminoacid	patientIdentifier   rnaExpression   rnaVariantAlleleFrequency   dnaVariantAlleleFrequency
hg19	BRCA2	uc003kii.3	AAAAAA	L	935	AAAAA	F	Pt1 4.512   0.4675  0.36103
hg19	BRCA2	uc003kii.3	AAAAAA	M	518	AAAAA	R	Pt2 0.154   0.015404    0.034404
hg19	BRCA2	uc003kii.3	AAAAAA	G	285	AAAAA	K	Pt3 8.841207    0.89387 0.51924
```
where:
- `gene.assembly` - the assembly of the reference genome (only hg19 is supported)
- `gene.gene` - the HGMC gene symbol
- `gene.transcriptIdentifier` - the UCSC transcript identifier including the version number
- `mutation.leftFlankingRegion` - the aminoacids flanking the mutation on the left (in IUPAC one letter symbols)
- `mutation.mutatedAminoacid` - the mutated aminoacid (IUPAC 1 or 3 letters respecting casing, eg: A and Ala)
- `mutation.position` - the 1 based position of the mutation in the protein
- `mutation.rightFlankingRegion` - the aminoacids flanking the mutation on the left (in IUPAC one letter symbols)
- `mutation.wildTypeAminoacid` - the wild type aminoacid (IUPAC 1 or 3 letters respecting casing, eg: A and Ala)
- `patientIdentifier` - the patient identifier
- `rnaExpression` - the transcript RNA expression (**optional**)
- `rnaVariantAlleleFrequency` - the variant allele frequency calculated from the RNA (**optional**, this will be estimated using the `dnaVariantAlleleFrequency` if not available)
- `dnaVariantAlleleFrequency` - the variant allele frequency calculated from the DNA (**optional**)

Example of patients data table:
```
identifier  mhcIAlleles mhcIIAlleles   isRnaAvailable  
Pt29    HLA-A*03:01,HLA-A*02:01,HLA-B*07:02 HLA-DRB1*11:04,HLA-DRB1*15:01  True    
```
where:
- `identifier` - the patient identifier
- `mhcIAlleles` - the list of MHC I alleles in the patient
- `mhcIIAlleles` - the list of MHC II alleles in the patient
- `isRnaAvailable` - whether RNA was available for the analysis. If true then the `VAF_in_RNA` field will be used, else `VAF_in_DNA` will be used. (**optional**)

### Output data

The output data is returned in a short wide tab separated values file (`--with-short-wide-table`). Optionally, it can be provided in a tall skinny tab separated values file (`--with-tall-skinny-table`) or in JSON (`--with-json`).

## Developer guide

### Build the package

To build the package just run:
```
python setup.py bdist_wheel
```

This will create an installable wheel file under `dist/neofox-x.y.z.whl`.

### Install the package

Install the wheel file as follows:
```
pip install dist/neofox-x.y.z.whl
```

### Run integration tests

To run the integration tests make sure you have a file `.env` that contains the environment variables described in the configuration section.

Run the integration tests as follows:
```
python -m unittest discover neofox.tests.integration_tests
```

The integration tests run over some real datasets and they take some time to run.

The integration test that runs the whole program over a relevant dataset can be run as follows:
```
python -m unittest neofox.tests.integration_tests.test_neofox
```

#### Regression tests

This last test (ie: `test_neofox`) writes its output to a file named `neofox/tests/resources/output_yyyymmddHHMMSS.txt`. If there is an existing file named `neofox/tests/resources/output_previous.txt` then it loads both files in memory and compares them. It outputs whether there are some lost or gained columns and for the common columns it evaluates if the values are the same. If they are the same the file `output_previous.txt` is overwritten by the new file, otherwise it outputs the details of the differing columns.

### Run unit tests

The unit tests do not have any dependency and they finish in seconds.

Run the unit tests as follows:
```
python -m unittest discover neofox.tests.unit_tests
```

### Logging

Logs are written to the standard error and to the output folder by default. Optionally they can be written to a file by setting the environment variable `NEOFOX_LOGFILE` pointing to the desired file.
