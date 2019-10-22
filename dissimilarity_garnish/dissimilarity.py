#!/usr/bin/env python

import sys
import os
import os.path
import subprocess
import sys
import tempfile


my_path = os.path.abspath(os.path.dirname(__file__))
my_path2 = "/".join(my_path.split("/")[0:-1])
sys.path.insert(0, my_path2)
sys.path.insert(0, my_path)

from neoantigen_fitness import Aligner_modified

def calc_dissimilarity(fasta_file, n):
    '''
    This function determines the dissimilarity to self-proteome of epitopes as described in Richman et al
    '''
    outfile_file = tempfile.NamedTemporaryFile(prefix ="tmp_prot_", suffix = ".xml", delete = False)
    outfile = outfile_file.name
    cmd = "/code/ncbi-blast/2.8.1+/bin/blastp -gapopen 11 -gapextend 1 -outfmt 5 -query " + fasta_file + " -out "  + outfile + " -db " + os.path.join(my_path,"proteome_db","homo_sapiens.mod") + " -evalue 100000000"
    p = subprocess.Popen(cmd.split(" "),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    p_return = p.communicate()
    a = Aligner_modified.Aligner()
    # set a to 32 for dissimilarity
    a.readAllBlastAlignments(outfile)
    a.computeR(a = 32)
    kk = int(n.split("_")[1])
    x = 0.0
    x = a.Ri.get(kk)
    x_dis = 1 - x
    os.remove(fasta_file)
    os.remove(outfile)
    return x_dis if x_dis != None else "NA"


def wrap_dissimilarity(props, fastafile, filter_binder = False):
    '''wrapper for dissimilarity calculation
    '''
    mhc_mut = props["best_affinity_epitope_netmhcpan4"]
    mhc_aff = props["best_affinity_netmhcpan4"]
    #fastafile = os.path.join(my_path,"tmp.fasta")
    with open(fastafile , "w") as f:
        id = ">M_1"
        f.write(id + "\n")
        f.write(mhc_mut + "\n")
    dissim = calc_dissimilarity(fastafile, id)
    #return str(pathsim)
    if filter_binder:
        if float(mhc_aff) < 500:
            sc = str(dissim) if dissim != "NA" else "0"
        else:
            sc = "NA"
    else:
        sc = str(dissim) if dissim != "NA" else "0"
    return sc





if __name__ == '__main__':
    import sys
    basedir = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2"
    sys.path.append(basedir)
    import data_import

    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/20170713_IS_IM_data.complete.update_Dv10.csv"
    indel = False
    dat = data_import.import_dat_icam(file, indel)

    properties = {}
    for nam,char in zip(dat[0], dat[1][641]):
        properties[nam] = char


    print wrap_pathogensimilarity(properties, "mhcI")
    print amplitude_mhc(properties, "mhcI")
    #print position_of_mutation_epitope(properties, "mhcI")
