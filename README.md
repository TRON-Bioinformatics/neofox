# **INPuT - Immunogenictiy Neoantigen Prediction Toolbox**


Annotation of mutated peptide sequences (mps) with published or novel potential neo-epitope descriptors

**Published Descriptors:** \
IEDB immunogenicity (Calis et al, 2013, PLoS Comput Biol.) \
Self-similarity, Conserved vs. Improved Binding  (Bjerregaard et al, 2017, Front Immunol.) \
Priority Score (Bjerregaard et al, 2017, Cancer Immunol Immunother.) \
DAI (Duan et al., 2014, JEM; Ghorani et al., 2018, Ann Oncol.) \
Neoantigen Fitness (Luksza et al., 2017, Nature; Balachandran et al, 2017, Nature) \
Residue-centric presentation score (best_rank) & Patient harmonic Best Rank (PHBR-I/II)* (Marty et al, 2017, Cell; Marty et al, 2018, Cell) \
Classically vs Alternatively Defined Neopitopes & Generator Rate* (Rech et al., 2018, Cancer Immunology Research) \
Tcell_predictor (Besser et al, 2019, Journal for ImmunoTherapy of Cancer) \


**Novel Potential Descriptors:** \
Amnino Acid Index \
Differential Expression \
Amino acid Frequency \
Conservation Scores (e.g PROVEAN: Choi et al, 2012, PLoS One) \
Multiplexed Representation \


## **Requirements**

**Specific Input:** \
allele.csv --> ";" separated file with mhc I and mhc II alleles for all patients of a cohort \


**Required Columns of iCaM Table:** \
MHC_I_epitope_.best_prediction. \
	MHC_I_epitope_.WT. \
	MHC_II_epitope_.best_prediction. \
	MHC_II_epitope_.WT. \
	MHC_I_score_.best_prediction. \
	MHC_I_score_.WT. \
	MHC_II_score_.best_prediction. \
	MHC_II_score_.WT. \
	MHC_I_peptide_length_.best_prediction. \
	MHC_I_allele_.best_prediction. \
	MHC_II_allele_.best_prediction. \
	transcript_expression \
	VAF_in_RNA \
	VAF_in_tumor \
	X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL. \

**Required Additional Files:** \
RNA reference \
n-mer frequencies \
PROVEAN score matrix \
available HLA I alleles for netmhcpan4 \
available HLA II alleles for netmhcIIpan3.2 \



**Required Additional Tools:** \


## **Usage**
--> parallel mps annotation of patients of a cohort

sh start_annotation_multiple_patientfiles.sh cohort_folder_with_patient_icam_folders output_folder allele_table cohort_name
