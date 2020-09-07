#############################################################
##################TCGA NeoAg analysis########################
#############################################################

#Original analysis run in R v3.5.2

library(caret) #Original analysis run in v6.0-84
library(Peptides) #Original analysis run in v2.4
library(data.table) #Original analysis run in v1.12.0
library(doParallel) #Original analysis run in v1.0.14
registerDoParallel(1) #Can change to suitable number of threads

args = commandArgs(trailingOnly=TRUE)

path_to_file <- args[1]
neo_tab_path <- args[2]
outifle <- args[3]


#Input paths
# neo_tab_path = "./TCGA_neoAg_example.txt"
# GBM_model_path = "./NeoAg_immunogenicity/Final_gbm_model.rds"
#neo_tab_path = paste(path_to_file, "TCGA_neoAg_example.txt", sep = "")
GBM_model_path = paste(path_to_file, "Final_gbm_model.rds", sep = "/")
#neo_tab_path = "/projects/SUMMIT/WP1.2/input/development/neoag/tmp_neoag_kZe31d.txt"
#neo_tab_path = "/projects/SUMMIT/WP1.2/input/development/neoag/neoag-master/TCGA_neoAg_example.txt"


#Example input
neo_tab = fread(neo_tab_path)

#Function for generating independent variables for the GBM model
model_process = function(n){  
  c(
    ifelse(substr(neo_tab$mut_peptide[n],1,1) == "V",1,0),
    
    ifelse(substr(neo_tab$mut_peptide[n],nchar(neo_tab$mut_peptide[n]),nchar(neo_tab$mut_peptide[n])) == "V",1,0),
    
    ifelse(aaComp(substr(neo_tab$mut_peptide[n],nchar(neo_tab$mut_peptide[n]),nchar(neo_tab$mut_peptide[n])))[[1]][2] == 1,1,0),
    
    ifelse(aaComp(substr(neo_tab$Reference[n],neo_tab$peptide_variant_position[n],neo_tab$peptide_variant_position[n]))[[1]][8] == 1,1,0),
    
    (aaComp(substr(neo_tab$mut_peptide[n],neo_tab$peptide_variant_position[n],neo_tab$peptide_variant_position[n]))[[1]][2] -
       aaComp(substr(neo_tab$Reference[n],neo_tab$peptide_variant_position[n],neo_tab$peptide_variant_position[n]))[[1]][2]),
    
    ifelse("K" %in% unlist(strsplit(substr(neo_tab$mut_peptide[n],1,nchar(neo_tab$mut_peptide[n])-7),"|")),1,0),
    
    ifelse("V" %in% unlist(strsplit(substr(neo_tab$mut_peptide[n],1,3),"|")),1,0)
  )
}

#Multi-thread derivation of features
model_mat = foreach(n = 1:nrow(neo_tab), .combine = rbind) %dopar% model_process(n)

# only one epitope-df 
if(is(model_mat, "numeric")){
  model_mat <- as.data.frame(t(model_mat))
}

colnames(model_mat) = c("Absolute_position_1_V", "Last_position_V", "Last_position_Small", "Reference_AA_at_mutated_position_Basic", 
                        "Mutated_position_change_of_Small_feature", "Relative_site_1_K", "First_three_AA_V" ) 
##############################################################

#Formatting, binding input matrix with feature sset
neo_tab_final = cbind(neo_tab, model_mat)

#Read in the GBM R object, run on the matrix generated above
Final_model = readRDS(GBM_model_path)

#Predicting neoantigen immunogenicity scores from above GBM model
TCGA_predict = predict(Final_model, newdata = model_mat, type =  "raw")

TCGA_predict <- ifelse(length(TCGA_predict) == 0, NA, TCGA_predict )

cat(TCGA_predict)
