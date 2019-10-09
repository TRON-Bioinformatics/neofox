#!/usr/bin/env python

import sys
import os
import os.path
import subprocess
import Aligner_modified
import sys
import tempfile


my_path = os.path.abspath(os.path.dirname(__file__))


def calc_pathogensimilarity(fasta_file, mhc_allele, n):
    '''
    This function determines the PATHOGENSIMILARITY of epitopes according to Balachandran et al. using a blast search against the IEDB pathogenepitope database
    '''
    outfile_file = tempfile.NamedTemporaryFile(prefix ="tmp_iedb_", suffix = ".xml", delete = False)
    outfile = outfile_file.name
    #outfile = os.path.join(my_path, "tmp_iedb.xml")
    #cmd = "/kitty/code/ncbi-blast-2.2.29+/bin/blastp -outfmt 5 -query "+ fasta_file +" -out " + outfile + " -db " + os.path.join(my_path,"scratch","iedb_blast_db") + " -evalue 100000000"
    # short sequence blast
    #cmd = "/code/ncbi-blast-2.7.1+/bin/blastp -word_size 2 -gapopen 11 -threshold 11 -comp_based_stats 0 -window_size 15 -outfmt 5 -query " + fasta_file + " -out "  + outfile + " -db " + os.path.join(my_path,"scratch","iedb_blast_db") + " -evalue 100000000000000000000000000"
    # parameters used in publication
    cmd = "/code/ncbi-blast/2.8.1+/bin/blastp -gapopen 11 -gapextend 1 -outfmt 5 -query " + fasta_file + " -out "  + outfile + " -db " + os.path.join(my_path,"iedb","iedb_blast_db") + " -evalue 100000000"
    #cmd = "/code/ncbi-blast-2.7.1+/bin/blastp -gapopen -11 -gapextend -1 -outfmt 5 -query " + fasta_file + " -out "  + outfile + " -db " + os.path.join(my_path,"iedb","iedb_blast_db") + " -evalue 100000000"
    p = subprocess.Popen(cmd.split(" "),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    p_return = p.communicate()
    #if p_return[0] != "" or p_return[1] != "":
    #  print >> sys.stderr, p_return
    a = Aligner_modified.Aligner()
    a.readAllBlastAlignments(outfile)
    a.computeR()
    kk = int(n.split("_")[1])
    x = 0.0
    x = a.Ri.get(kk)
    os.remove(fasta_file)
    os.remove(outfile)
    return x if x != None else "NA"


def wrap_pathogensimilarity(props, mhc, fastafile, affinity = False, nine_mer = False):
    if mhc == "mhcI":
        if nine_mer:
            mhc_mut = props["best_affinity_epitope_netmhcpan4_9mer"]
        elif affinity:
            mhc_mut = props["best_affinity_epitope_netmhcpan4"]
        else:
            mhc_mut = props["MHC_I_epitope_.best_prediction."]
    elif mhc == "mhcII":
        mhc_mut = props["MHC_II_epitope_.best_prediction."]
    #fastafile = os.path.join(my_path,"tmp.fasta")
    with open(fastafile , "w") as f:
        id = ">M_1"
        f.write(id + "\n")
        f.write(mhc_mut + "\n")
    pathsim = calc_pathogensimilarity(fastafile, mhc, id)
    #return str(pathsim)
    return str(pathsim) if pathsim != "NA" else "0"


def amplitude_mhc(props, mhc, multiple_binding=False, affinity = False, netmhcscore = False, nine_mer= False):
    '''
    This function calculates the amplitude between mutated and wt epitope according to Balachandran et al.
    when affinity is used, use correction from Luksza et al. *1/(1+0.0003*aff_wt)
    '''
    if mhc == "mhcI":
        if multiple_binding:
            sc_mut = props["MB_score_top10_harmonic"].replace(",",".")
            sc_wt = props["MB_score_WT_top10_harmonic"].replace(",",".")
        elif affinity:
            sc_mut = props["best_affinity_netmhcpan4"].replace(",",".")
            sc_wt = props["best_affinity_netmhcpan4_WT"].replace(",",".")
        elif netmhcscore:
            sc_mut = props["best%Rank_netmhcpan4"].replace(",",".")
            sc_wt = props["best%Rank_netmhcpan4_WT"].replace(",",".")
        elif nine_mer:
            sc_mut = props["best_affinity_netmhcpan4_9mer"].replace(",",".")
            sc_wt = props["best_affinity_netmhcpan4_9mer_WT"].replace(",",".")
        else:
            sc_mut = props["MHC_I_score_.best_prediction."].replace(",",".")
            sc_wt = props["MHC_I_score_.WT."].replace(",",".")
    elif mhc == "mhcII":
        if multiple_binding:
            sc_mut = props["MB_score_MHCII_top10_harmonic"].replace(",",".")
            sc_wt = props["MB_score_MHCII_top10_WT_harmonic"].replace(",",".")
        elif affinity:
            sc_mut = props["best_affinity_netmhcIIpan"].replace(",",".")
            sc_wt = props["best_affinity_netmhcIIpan_WT"].replace(",",".")
        elif netmhcscore:
            sc_mut = props["best%Rank_netmhcIIpan"].replace(",",".")
            sc_wt = props["best%Rank_netmhcIIpan_WT"].replace(",",".")
        else:
            sc_mut = props["MHC_II_score_.best_prediction."].replace(",",".")
            sc_wt = props["MHC_II_score_.WT."].replace(",",".")
    try:
        if nine_mer or affinity:
            return str(float(sc_wt) / float(sc_mut) * (1 / (1 + 0.0003 * float(sc_wt))))
        else:
            return str(float(sc_wt) / float(sc_mut))
    except(ZeroDivisionError, ValueError) as e:
        return "NA"


def recognition_potential(props, mhc, affinity = False,  netmhcscore = False, nine_mer= False ):
    '''
    This function calculates the recognition potential, defined by the product of amplitude and pathogensimiliarity of an epitope according to Balachandran et al.
    F_alpha = - max (A_i x R_i)

    Returns (A_i x R_i) value only for nonanchor mutation and epitopes of length 9; only considered by Balachandran
    '''
    if mhc == "mhcI":
        if affinity:
            amp = props["Amplitude_mhcI_affinity"]
            pathsim = props["Pathogensimiliarity_mhcI_affinity_nmers"]
        elif netmhcscore:
            amp = props["Amplitude_mhcI_rank_netmhcpan4"]
            pathsim = props["Pathogensimiliarity_mhcI"]
        elif nine_mer:
            amp = props["Amplitude_mhcI_affinity_9mer_netmhcpan4"]
            pathsim = props["Pathogensimiliarity_mhcI_9mer"]
        else:
            amp = props["Amplitude_mhcI"]
            pathsim = props["Pathogensimiliarity_mhcI"]
    elif mhc == "mhcII":
        amp = props["Amplitude_mhcII"]
        pathsim = props["Pathogensimiliarity_mhcII"]
    try:
        if nine_mer:
            mhc_affinity_mut = props["best_affinity_netmhcpan4_9mer"]
            recog = str(float(amp) * float(pathsim)) if (props["Mutation_in_anchor"] == "0") & (float(mhc_affinity_mut) < 500) else "NA"
        else:
            recog = str(float(amp) * float(pathsim)) if props["Mutation_in_anchor"] == "0" else "NA"
    except ValueError:
        recog = "NA"
    #print >> sys.stderr, props["Mutation_in_anchor"], str(recog)
    return str(recog)


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
