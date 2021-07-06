require(Biostrings)

args = commandArgs(trailingOnly=TRUE)

epitope <- args[1]
path_virus_db <- args[2]
path_to_file <- args[3]

source(paste0(path_to_file,"/EPItOMe_modified.R"))
load(paste0(path_to_file, "/BLOSUM62.rda"))
  
aln_matrix <- "BLOSUM62"

# neoepitope candidate sequence
qAASet <- AAStringSet(epitope)
qAASet <- Biostrings::as.matrix(qAASet, use.names=T)

# virus database
#path_virus_db <- "/home/franlang/neofox_test/test_references_nets/iedb/IEDB.fasta"
rAASet <- readAAStringSet(path_virus_db, format="fasta")
# do filtering based on length of neoepitope sequence
rAASet <- rAASet[which(rAASet@ranges@width == length(qAASet))]
rms <- which(alphabetFrequency(rAASet)[,31] != 0)
# remove viral epitopes with non-canonical aa
if(length(rms) > 0){
  rAASet <- rAASet[-rms]
}
rAASet <- Biostrings::as.matrix(rAASet, use.names=T)


# alignment of viral and tumor epitopes 
# alignment score
alnRes <- align_sets(query.set=qAASet, ref.set=rAASet, aln_matrix=aln_matrix)
# B-score
#cplRes <- align_cpl_sets(query.set=qAASet, ref.set=rAASet, aln_matrix=aln_matrix)

#res <- c(unique(max(alnRes)), unique(max(cplRes)))
res <- unique(max(alnRes, na.rm = T))
if(!is.numeric(res)){res = ""}

cat(res)