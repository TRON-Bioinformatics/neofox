#!/bin/bash


path_to_mixMHC2pred=`echo $INPUT_MIXMHC2PRED |sed 's/\(.*\)MixMHC2pred/\1/'`
echo $path_to_mixMHC2pred

# available MHC alleles netMHCpan
/code/net/MHCpan/4.0/netMHCpan -listMHC | grep "HLA-" > "$INPUT_REFERENCE_FOLDER"/MHC_available.csv


# available MHCII alleles netMHCIIpan
/code/net/MHCIIpan/3.2/netMHCIIpan -list  > "$INPUT_REFERENCE_FOLDER"/avail_mhcII.txt


# available MHCII alleles for MixMHC2pred
cp "$path_to_mixMHC2pred"/Alleles_list.txt $INPUT_REFERENCE_FOLDER

# build IEDB blast database
mkdir "$INPUT_REFERENCE_FOLDER"/iedb
wget "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip" -O "$INPUT_REFERENCE_FOLDER"/iedb/Iedb.zip
unzip "$INPUT_REFERENCE_FOLDER"/iedb/Iedb.zip -d "$INPUT_REFERENCE_FOLDER"/iedb/
$INPUT_RSCRIPT "$INPUT_REFERENCE_FOLDER"/build_IEDB_db.R "$INPUT_REFERENCE_FOLDER"/iedb/
$INPUT_MAKEBLASTDB -in "$INPUT_REFERENCE_FOLDER"/iedb/IEDB.fasta -dbtype prot -parse_seqids -out "$INPUT_REFERENCE_FOLDER"/iedb/iedb

# human proteome database
mkdir "$INPUT_REFERENCE_FOLDER"/proteom_db
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O "$INPUT_REFERENCE_FOLDER"/proteom_db/Homo_sapiens.fa.gz
gunzip "$INPUT_REFERENCE_FOLDER"/proteom_db/Homo_sapiens.fa.gz
$INPUT_MAKEBLASTDB -in "$INPUT_REFERENCE_FOLDER"/proteom_db/Homo_sapiens.fa -dbtype prot -parse_seqids -out "$INPUT_REFERENCE_FOLDER"/proteom_db/homo_sapiens

# use Homo_sapiens.fa as uniprot_human_with_isoforms.fasta to be consistent with human proteome database
# provide BLOSUM62-2.matrix.txt as file and refer to "https://arxiv.org/abs/1205.6031"
