#!/bin/bash


# available MHC alleles netMHCpan
$INPUT_NETMHCPAN -listMHC | grep "HLA-" > "$INPUT_REFERENCE_FOLDER"/MHC_available.csv


# available MHCII alleles netMHCIIpan
$INPUT_NETMHC2PAN -list  > "$INPUT_REFERENCE_FOLDER"/avail_mhcII.txt

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
