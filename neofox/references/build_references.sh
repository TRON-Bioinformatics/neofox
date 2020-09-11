#!/bin/bash
##
## Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
##
## This file is part of Neofox
## (see https://github.com/tron-bioinformatics/neofox).
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.##


# available MHC alleles netMHCpan
$NEOFOX_NETMHCPAN -listMHC | grep "HLA-" > "$NEOFOX_REFERENCE_FOLDER"/MHC_available.csv

# available MHCII alleles netMHCIIpan
$NEOFOX_NETMHC2PAN -list  > "$NEOFOX_REFERENCE_FOLDER"/avail_mhcII.txt

# build IEDB blast database
mkdir "$NEOFOX_REFERENCE_FOLDER"/iedb
wget "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip" -O "$NEOFOX_REFERENCE_FOLDER"/iedb/Iedb.zip
unzip "$NEOFOX_REFERENCE_FOLDER"/iedb/Iedb.zip -d "$NEOFOX_REFERENCE_FOLDER"/iedb/
$NEOFOX_RSCRIPT "$NEOFOX_REFERENCE_FOLDER"/build_IEDB_db.R "$NEOFOX_REFERENCE_FOLDER"/iedb/
$NEOFOX_MAKEBLASTDB -in "$NEOFOX_REFERENCE_FOLDER"/iedb/IEDB.fasta -dbtype prot -parse_seqids -out "$NEOFOX_REFERENCE_FOLDER"/iedb/iedb

# human proteome database
mkdir "$NEOFOX_REFERENCE_FOLDER"/proteom_db
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O "$NEOFOX_REFERENCE_FOLDER"/proteome_db/Homo_sapiens.fa.gz
gunzip "$NEOFOX_REFERENCE_FOLDER"/proteom_db/Homo_sapiens.fa.gz
$NEOFOX_MAKEBLASTDB -in "$NEOFOX_REFERENCE_FOLDER"/proteom_db/Homo_sapiens.fa -dbtype prot -parse_seqids -out "$NEOFOX_REFERENCE_FOLDER"/proteome_db/homo_sapiens

# use Homo_sapiens.fa as uniprot_human_with_isoforms.fasta to be consistent with human proteome database


# PROVEAN scores
wget ftp://ftp.jcvi.org/pub/data/provean/precomputed_scores/PROVEAN_scores_ensembl66_human.tsv.gz
# TODO: some mapping of Ensembl protein ids to UCSC ids happening here. Do we still have this code?
# eventually the mapped resource becomes comma separated...
tr ';' '\t' < PROV_scores_mapped3.csv > PROV_scores_mapped3.tab
# those Ensembl ids without a mapping to UCSC are stored as NA we need to filter out those
grep -v NA PROV_scores_mapped3.tab > PROV_scores_mapped3.filtered.tab
# we need to sort by the UCSC id and position (keeping the header)
(head -n 1 PROV_scores_mapped3.filtered.tab && tail -n +2 PROV_scores_mapped3.filtered.tab | sort -k25,25V -k2,2n) > PROV_scores_mapped3.filtered.sorted.tab
# bgzip and tabix index on the UCSC id and the position in the protein
bgzip PROV_scores_mapped3.filtered.sorted.tab
tabix -f -s 25 -b 2 -e 2 PROV_scores_mapped3.filtered.sorted.tab.gz