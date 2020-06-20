# **INPuT - Immunogenictiy Neoantigen Prediction Toolbox**


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


## Input Requirements

**Specific Input:**
- icam_output.txt --> icam output file
- patient identifier --> the patient identifier to whom all neoantigens in icam output belong
- patient data --> a table of tab separated values containing metadata on the patient
  - required fields: identifier, mhcIAlleles, mhcIIAlleles
  - optional fields: estimatedTumorContent, isRnaAvailable, tissue

**Example of patient data table**
```
identifier  mhcIAlleles mhcIIAlleles    estimatedTumorContent   isRnaAvailable  tissue
Pt29    HLA-A*03:01,HLA-A*02:01,HLA-B*07:02 HLA-DRB1*11:04,HLA-DRB1*15:01   69  True    skin
```

**Required Columns of iCaM Table:**  
-   MHC_I_epitope_.best_prediction.  
- 	MHC_I_epitope_.WT.  
-   MHC_II_epitope_.best_prediction.  
- 	MHC_II_epitope_.WT.  
- 	MHC_I_score_.best_prediction.  
- 	MHC_I_score_.WT.  
- 	MHC_II_score_.best_prediction.  
- 	MHC_II_score_.WT.  
- 	MHC_I_peptide_length_.best_prediction.
- 	MHC_I_allele_.best_prediction.  
- 	MHC_II_allele_.best_prediction.  
- 	transcript_expression  
- 	VAF_in_RNA  
- 	VAF_in_tumor  
- 	X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL.
-   substitution  
-   patient.id (e.g Pt1, Ptx)    

**Required Additional Files:**  
- RNA reference *(/projects/CM27_IND_patients/GTEX_normal_tissue_data/Skin .csv, predict_all_epitopes.py)*  
- protein database *(/projects/data/human/2018_uniprot_with_isoforms/uniprot_human_with_isoforms.fasta, predict_all_epitopes.py)*  
- amino acid frequencies *(./new_features/20181108_AA_freq_prot.csv, predict_all_epitopes.py)*  
- 4mer amino acid frequnecies *(./new_features/20181108_4mer_freq.csv, predict_all_epitopes.py)*  
- PROVEAN score matrix *(./new_features/PROV_scores_mapped3.csv, predict_all_epitopes.py)*  
- available HLA I alleles for netmhcpan4 *(./netmhcpan4/MHC_available.csv, predict_all_epitopes.py)*  
- available HLA II alleles for netmhcIIpan3.2 *(./netmhcIIpan/avail_mhcII.txt, predict_all_epitopes.py)*  
- aaindex1 *("aa_index/aaindex1", predict_all_epitopes.py)*  
- aanindex2 *("aa_index/aaindex1", predict_all_epitopes.py)*  
- available HLA II alleles for MixMHC2pred *("/projects/SUMMIT/WP1.2/input/development/MixMHCpred/Alleles_list_pred2.txt")*

**Required Software/Tools/Dependencies:**  
- python2 *(anaconda/2/2018)*
- BLAST *(/code/ncbi-blast/2.8.1+/bin/blastp, neoantigen_fitness.py)*  
- netmhcpan *(/code/netMHCpan-4.0/netMHCpan, netmhcpan_prediction.py)*  
- netmhcIIpan *(/code/net/MHCIIpan/3.2/netMHCIIpan, netmhcIIpan_prediction.py)*  
- netmhcIIpan *(/code/net/MHCIIpan/3.2/netMHCIIpan, netmhcIIpan_prediction.py)*  
- MixMHCpred *(/code/MixMHCpred/2.0.2/MixMHCpred, mixmhcpred.py)*
- Tcell_predictor: python3 + scripts/pickle/mat files of Tcell_predictor tool *(/code/Anaconda/3/2018/bin/python + tool under ./Tcell_predictor, tcellpredictor_wrapper.py )*  
- Neoag: Neoag R-module *(./neoag-master, neoag_gbm_model.py)* + R *(/code/R/3.6.0/bin/Rscript)*
- MixMHCpred *(/code/MixMHCpred/2.0.2/)*
- MixMHC2pred *(/code/net/MixMHC2pred/1.1)*

## **Usage**  

```
input --icam-file testseq_head.txt --patient-id Pt123 --patient-data patients.csv [--frameshift False]
```


## Developer guide

### Build the package

To build the package just run:
```
python setup.py bdist_wheel
```

This will create an installable wheel file under `dist/input-x.y.z.whl`.

### Install the package

Install the wheel file as follows:
```
pip install dist/input-x.y.z.whl
```

### Run integration tests

To run the integration tests make sure you have a file `.env` that contains the following variables with the right values:
```
export INPUT_REFERENCE_FOLDER=~/addannot_references
export INPUT_BLASTP=/code/ncbi-blast/2.8.1+/bin/blastp
export INPUT_MIXMHC2PRED=/code/net/MixMHC2pred/1.1/MixMHC2pred
export INPUT_MIXMHCPRED=/code/MixMHCpred/2.0.2/MixMHCpred
export INPUT_RSCRIPT=/code/R/3.6.0/bin/Rscript
export INPUT_NETMHC2PAN=/code/net/MHCIIpan/3.2/netMHCIIpan
export INPUT_NETMHCPAN=/code/net/MHCpan/4.0/netMHCpan
```

The folder `$INPUT_REFERENCE_FOLDER` requires to contain the resources defined above.

Run the integration tests as follows:
```
python -m unittest discover input.tests.integration_tests
```

The integration tests run over some real datasets and they take some time to run.

The integration test that runs the whle program over a relevant dataset can be run as follows:
```
python -m unittest input.tests.integration_tests.test_input
```

### Run unit tests

The unit tests do not have any dependency and they finish in seconds.

Run the unit tests as follows:
```
python -m unittest discover input.tests.unit_tests
```