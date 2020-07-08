
library(seqinr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
path_to_file <- args[1]
path_to_file <- "/home/franlang/refs_test/iedb"

d2 <- read_csv(paste0(path_to_file, "/tcell_full_v3.csv"), skip = 1)

# grep human MHC class I epitopes that were positive in context of infectious disease
d2.human <- d2 %>%
  filter(Name == "Homo sapiens",
         `Object Type` == "Linear peptide",
         `Process Type` == "Occurrence of infectious disease",
         `Qualitative Measure`== "Positive",
         Class == "I") %>%
  distinct(Description, .keep_all = T) %>%
  mutate(epitope_id = gsub("http://www.iedb.org/epitope/", "", `Epitope IRI`),
        #fasta_header = substr(paste(epitope_id, `Antigen Name`, `Parent Species`, sep = "_"), 1,50),
        fasta_header = epitope_id,
        seq = str_split_fixed(Description, " ", 2)[,1])

# write fasta file for DB
path_to_fasta = paste(path_to_file, "IEDB.fasta", sep =  "/")
write.fasta(sequences = as.list(d2.human$seq), names = d2.human$fasta_header, file.out = path_to_fasta)

