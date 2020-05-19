#!/bin/bash

# remove and adjust
module load anaconda/3/2019
INPUT_MIXMHC2PRED=/code/net/MixMHC2pred/1.1/MixMHC2pred
INPUT_MAKEBLASTDB=/code/ncbi-blast/2.8.1+/bin/makeblastdb
INPUT_RSCRIPT=/code/R/3.6.0/bin/Rscript
path_to_references="/home/franlang/refs_test/"


path_to_mixMHC2pred=`echo $INPUT_MIXMHC2PRED |sed 's/\(.*\)MixMHC2pred/\1/'`
echo $path_to_mixMHC2pred

# available MHC alleles netMHCpan
/code/net/MHCpan/4.0/netMHCpan -listMHC | grep "HLA-" > "$path_to_references"/MHC_available.csv


# available MHCII alleles netMHCIIpan
/code/net/MHCIIpan/3.2/netMHCIIpan -list  > "$path_to_references"/avail_mhcII.txt


# available MHCII alleles for MixMHC2pred
cp "$path_to_mixMHC2pred"/Alleles_list.txt $path_to_references

# build IEDB blast database
mkdir "$path_to_references"/iedb
wget "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip" -O "$path_to_references"/iedb/Iedb.zip
unzip "$path_to_references"/iedb/Iedb.zip -d "$path_to_references"/iedb/
$INPUT_RSCRIPT "$path_to_references"/build_IEDB_db.R "$path_to_references"/iedb/
$INPUT_MAKEBLASTDB -in "$path_to_references"/iedb/IEDB.fasta -dbtype prot -parse_seqids -out "$path_to_references"/iedb/iedb

# human proteome database
mkdir "$path_to_references"/proteom_db
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O "$path_to_references"/proteom_db/Homo_sapiens.fa.gz
gunzip "$path_to_references"/proteom_db/Homo_sapiens.fa.gz
$INPUT_MAKEBLASTDB -in "$path_to_references"/proteom_db/Homo_sapiens.fa -dbtype prot -parse_seqids -out "$path_to_references"/proteom_db/homo_sapiens

# use Homo_sapiens.fa as uniprot_human_with_isoforms.fasta to be consistent with human proteome database
# provide BLOSUM62-2.matrix.txt as file and refer to "https://arxiv.org/abs/1205.6031"
