**INPuT - Immunogenictiy Neoantigen Prediction Toolbox**


Annotation of mutated peptide sequences (mps) with neo-epitope descriptors

Required columns: MHC_I_epitope_.best_prediction. \
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


descriptors:
IEDB immunogenicity
self-similarity




**USAGE**
--> parallel mps annotation of patients of a cohort

sh start_annotation_multiple_patientfiles.sh cohort_folder_with_patient_icam_folders output_folder allele_table cohort_name
