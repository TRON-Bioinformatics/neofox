#!/usr/bin/python
import pandas as pd
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]

# read IEDB input file
iedb = pd.read_csv(input_file, skiprows=1)

# filter entries
filtered_iedb = iedb[
    (iedb["Name"].isin(["Homo sapiens", "Homo sapiens (human)", "Homo sapiens Caucasian", "Homo sapiens Black"])) &
    (iedb["Object Type"] == "Linear peptide") &
    (iedb["Process Type"] == "Occurrence of infectious disease") &
    (iedb["Qualitative Measure"] == "Positive") &
    (iedb["Class"] == "I")
]

# sets values for identifiers and sequences
filtered_iedb["seq"] = filtered_iedb["Description"].transform(lambda x: x.strip())
filtered_iedb["fasta_header"] = filtered_iedb["Epitope IRI"].transform(lambda x: x.replace("http://www.iedb.org/epitope/", ""))
filtered_iedb.drop_duplicates(subset="seq", keep="last", inplace=True)

# writes output FASTA file
with open(output_file, "w") as fasta:
    for index, row in filtered_iedb.iterrows():
        fasta.write(">{header}\n".format(header=str(row["fasta_header"])))
        fasta.write("{sequence}\n".format(sequence=str(row["seq"])))
