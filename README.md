# **NeoFox - NEOantigen Feature tOolboX**


Annotation of mutated peptide sequences (mps) with published or novel potential neo-epitope descriptors

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


**Novel Potential Descriptors:**  
- Amnino Acid Index  
- Differential Expression  
- Amino acid Frequency  
- Conservation Scores (e.g PROVEAN: Choi et al, 2012, PLoS One)  
- Multiplexed Representation  


## NeoFox Requirements
 
**Required Software/Tools/Dependencies:**  
- Python-3.7.3
- BLAST-2.8.1 
- netMHCpan-4.0 
- netMHCIIpan-3.2
- MixMHCpred 2.0.2  
- R-3.6.0
- MixMHC2pred 1.1 

## **Usage**  

```
neofox --model-file/--icam-file testseq_head.txt --patient-id Ptx --patient-data patient_data.txt --output-folder /path/to/out --output-prefix out_prefix [--with-short-wide-table] [--with-tall-skinny-table] [--with-json] [--num_cpus]
```
**Specific Input:**

--icam-file: tab-separated file in iCaM output style (**required columns**: transcript_expression,+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL),[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL),  VAF_in_RNA, substitution)  <br>
--model-file: file in model format<br>
--patient-id: patient identifier <br>
--patient data: a table of tab separated values containing metadata on the patient (**required fields**: identifier, mhcIAlleles, mhcIIAlleles; **optional fields**: estimatedTumorContent, isRnaAvailable, tissue)

Example of patient data table:
```
identifier  mhcIAlleles mhcIIAlleles    estimatedTumorContent   isRnaAvailable  tissue
Pt29    HLA-A*03:01,HLA-A*02:01,HLA-B*07:02 HLA-DRB1*11:04,HLA-DRB1*15:01   69  True    skin
```


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

To run the integration tests make sure you have a file `.env` that contains the following variables with the right values:
```
export NEOFOX_REFERENCE_FOLDER=~/addannot_references
export NEOFOX_BLASTP=/code/ncbi-blast/2.8.1+/bin/blastp
export NEOFOX_MIXMHC2PRED=/code/net/MixMHC2pred/1.1/MixMHC2pred
export NEOFOX_MIXMHCPRED=/code/MixMHCpred/2.0.2/MixMHCpred
export NEOFOX_RSCRIPT=/code/R/3.6.0/bin/Rscript
export NEOFOX_NETMHC2PAN=/code/net/MHCIIpan/3.2/netMHCIIpan
export NEOFOX_NETMHCPAN=/code/net/MHCpan/4.0/netMHCpan
```

The folder `NEOFOX_REFERENCE_FOLDER` requires to contain the resources defined above.

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

Logs are written to the standard error by default. Optionally they can be written to a file by setting the environment variable `NEOFOX_LOGFILE` pointing to the desired file.
