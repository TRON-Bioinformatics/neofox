#Python version with Biopython
#. /etc/profile.d/modules.sh; module load software/python/python-2.7.9

'''
This script takes as input table from iCAM pipeline and calculates Literature feature of neoantigens
'''

# import modules
import sys
import os
import math
from Bio import SeqIO

from aa_index import aa_index
from IEDB_Immunogenicity import predict_immunogenicity_simple

def calc_IEDB_immunogenicity(props, mhc, affin_filtering = False):
    '''
    This function determines the IEDB immunogenicity score
    '''
    if mhc == "mhcI":
        #mhc_mut = props["MHC_I_epitope_.best_prediction."]
        #mhc_allele = props["MHC_I_allele_.best_prediction."]
        mhc_mut = props["best_affinity_epitope_netmhcpan4"]
        mhc_allele = props["bestHLA_allele_affinity_netmhcpan4"]
        mhc_score = props["best_affinity_netmhcpan4"]
    elif mhc == "mhcII":
        mhc_mut = props["MHC_II_epitope_.best_prediction."]
        mhc_allele = props["MHC_II_allele_.best_prediction."]
    try:
        if affin_filtering:
            if float(mhc_score) < 500:
                return str(predict_immunogenicity_simple.predict_immunogenicity(mhc_mut,mhc_allele.replace("*","").replace(":","")))
            else:
                return "NA"
        else:
            return str(predict_immunogenicity_simple.predict_immunogenicity(mhc_mut,mhc_allele.replace("*","").replace(":","")))
    except ValueError:
        return "NA"

def dai(props, mhc, multiple_binding=False, affinity = False, netmhcscore = False, affin_filtering = False):
    '''Calculates DAI: Returns difference between wt and mut MHC binding score. If multiple_binding= true, harmonic means of MHC scores of top10 epitope candidates related to a mps is used '''
    if mhc == "mhcI":
        if multiple_binding:
            sc_mut = props["MB_score_top10_harmonic"]
            sc_wt = props["MB_score_WT_top10_harmonic"]
        elif affinity:
            sc_mut = props["best_affinity_netmhcpan4"]
            sc_wt = props["best_affinity_netmhcpan4_WT"]
        elif netmhcscore:
            sc_mut = props["best%Rank_netmhcpan4"]
            sc_wt = props["best%Rank_netmhcpan4_WT"]
        else:
            sc_mut = props["MHC_I_score_.best_prediction."]
            sc_wt = props["MHC_I_score_.WT."]
    elif mhc == "mhcII":
        if multiple_binding:
            sc_mut = props["MB_score_MHCII_top10_harmonic"]
            sc_wt = props["MB_score_MHCII_top10_WT_harmonic"]
        elif affinity:
            sc_mut = props["best_affinity_netmhcIIpan"]
            sc_wt = props["best_affinity_netmhcIIpan_WT"]
        elif netmhcscore:
            sc_mut = props["best%Rank_netmhcIIpan"]
            sc_wt = props["best%Rank_netmhcIIpan_WT"]
        else:
            sc_mut = props["MHC_II_score_.best_prediction."]
            sc_wt = props["MHC_II_score_.WT."]
    try:
        if affin_filtering:
            if float(sc_mut) < 500:
                return str(float(sc_wt) - float(sc_mut))
            else:
                return "NA"
        else:
            return str(float(sc_wt) - float(sc_mut))
    except ValueError:
        return "NA"

def diff_number_binders(props, mhc = "mhcI", threshold = 1):
    ''' returns absolute difference of potential candidate epitopes between mutated and wt epitope
    '''
    if mhc =="mhcII":
        num_mut = props["MB_number_pep_MHCIIscore<" + str(threshold)]
        num_wt =  props["MB_number_pep_MHCIIscore<" + str(threshold) + "_WT"]
    else:
        num_mut = props["MB_number_pep_MHCscore<" + str(threshold)]
        num_wt =  props["MB_number_pep_WT_MHCscore<" + str(threshold)]
    try:
        return str(float(num_mut) - float(num_wt))
    except ValueError:
        return "NA"

def ratio_number_binders(props, mhc = "mhcI", threshold = 1):
    ''' returns ratio of number of potential candidate epitopes between mutated and wt epitope. if no WT candidate epitopes, returns number of mutated candidate epitopes per mps
    '''
    if mhc =="mhcII":
        num_mut = props["MB_number_pep_MHCIIscore<" + str(threshold)]
        num_wt =  props["MB_number_pep_MHCIIscore<" + str(threshold) + "_WT"]
    else:
        num_mut = props["MB_number_pep_MHCscore<" + str(threshold)]
        num_wt =  props["MB_number_pep_WT_MHCscore<" + str(threshold)]
    try:
        return str(float(num_mut) / float(num_wt))
    except ZeroDivisionError:
        return str(num_mut)
    except ValueError:
        return "NA"

def rna_expression_mutation(props):
    '''
    This function calculates the product of VAF in RNA and transcript expression
    to reflect the expression of the mutated transcript
    '''
    transcript_expression = props["transcript_expression"]
    try:
        vaf_rna = props["VAF_in_RNA"]
    except KeyError:
        vaf_rna = props["VAF_in_tumor"]
    try:
        return str(float(transcript_expression) * float(vaf_rna)) if float(vaf_rna) > 0 else "NA"
    except ValueError:
        return "NA"

def expression_mutation_tc(props, tumour_content):
    '''calculated expression of mutation corrected by tumour content
    '''
    transcript_expression = props["Expression_Mutated_Transcript"]
    if "patient.id" in props:
        patid = props["patient.id"]
    else:
        patid = props["patient"]
    try:
        tumour_content =  float(tumour_content[patid])/100
    except (KeyError, ValueError) as e:
        tumour_content = "NA"
    try:
        return str(float(transcript_expression) / float(tumour_content))
    except ValueError:
        return "NA"


def number_of_mismatches(props, mhc):
    '''
    This function calculates the number of mismatches between the wt and the mutated epitope
    '''
    if mhc == "mhcI":
        mhc_epitope_mut = props["MHC_I_epitope_.best_prediction."]
        mhc_epitope_wt = props["MHC_I_epitope_.WT."]
    elif mhc == "mhcII":
        mhc_epitope_mut = props["MHC_II_epitope_.best_prediction."]
        mhc_epitope_wt = props["MHC_II_allele_.best_prediction."]
    p1 = 0
    try:
        for i,aa in enumerate(mhc_epitope_mut):
            if aa != mhc_epitope_wt[i]:
                p1 += 1
        return str(p1)
    except IndexError:
        return "NA"

def match_in_proteome(props, db):
    '''
    This function checks if the mutated epitope has an exact match in a protein database (uniprot)
    Returns 0 if mutation is present in proteome and 1 if it not present
    '''
    seq = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
    try:
        return "0" if seq in db else "1"
    except:
        return "NA"

def calc_logistic_function(mhc_score):
    '''Calculates negative logistic function given mhc score
    '''
    try:
        return float(1/(1 + math.exp(5 * (float(mhc_score) - 2))))
    except (OverflowError, ValueError) as e:
        return "NA"

def calc_priority_score(props, multiple_binding=False):
    '''
    This function calculates the Priority Score using parameters for mhc I.
    '''
    vaf_tumor = props["VAF_in_tumor"]
    try:
        vaf_rna = props["VAF_in_RNA"]
    except KeyError:
        vaf_tumor = props["VAF_in_tumor"]
    transcript_expr = props["transcript_expression"]
    no_mismatch = props["Number_of_mismatches_mhcI"]
    if multiple_binding:
        score_mut = props["MB_score_top10_harmonic"]
        score_wt = props["MB_score_WT_top10_harmonic"]
    else:
        #score_mut = props["MHC_I_score_.best_prediction."]
        #score_wt = props["MHC_I_score_.WT."]
        score_mut = props["best%Rank_netmhcpan4"]
        score_wt = props["best%Rank_netmhcpan4_WT"]
    mut_in_prot = props["mutation_found_in_proteome"]
    mut_in_prot = "0" if mut_in_prot == "True" else "1" if mut_in_prot == "False" else mut_in_prot
    L_mut = calc_logistic_function(score_mut)
    L_wt = calc_logistic_function(score_wt)
    priority_score = 0.0
    try:
        if vaf_tumor != "-1":
            priority_score = (L_mut * float(vaf_tumor) * math.tanh(float(transcript_expr))) * (float(mut_in_prot) * (1 - 2**(-float(no_mismatch)) * L_wt))
        else:
            priority_score = (L_mut * float(vaf_rna) * math.tanh(float(transcript_expr))) * (float(mut_in_prot) * (1 - 2**(-float(no_mismatch)) * L_wt))
        return str(priority_score)
    except (TypeError, ValueError) as e:
        return "NA"

def wt_mut_aa(props, mut):
    '''Returns wt and mut aa.
    '''
    substitution = props["substitution"]
    try:
        if mut == "mut":
            return substitution[-1]
        elif mut == "wt":
            return substitution[0]
    except ValueError:
        return "NA"

def add_aa_index1(props, mut, key, val):
    """Adds amino acido index to dictioniary. output from aa index 1 append function = tuple of feature name and feature value (nam_wt, nam_mut, val_wt, val_mut)
    mut indicates if mutated or wt aa
    """
    if mut == "mut":
        aa = props["MUT_AA"]
    elif mut == "wt":
        aa = props["WT_AA"]
    try:
        return "_".join([key, mut]), str(val[aa])
    except KeyError:
            return "_".join([key, mut]), "NA"

def add_aa_index2(props, key, val):
    """Adds amino acido index to dictioniary. output from aa index 1 append function = tuple of feature name and feature value (nam_wt, nam_mut, val_wt, val_mut)
    mut indicates if mutated or wt aa
    """
    mut_aa = props["MUT_AA"]
    wt_aa = props["WT_AA"]
    try:
        return key, str(val[wt_aa][mut_aa])
    except KeyError:
        return key, "NA"


def write_ouptut_to_file(epitope_data):
    '''
    This function prints output, semilicon separated --> csv file
    '''
    dat_new = epitope_data[1]
    head_new = epitope_data[0]
    print ";".join(head_new)
    for ii,i in enumerate(dat_new):
          print ";".join(i)


def classify_adn_cdn(props, mhc, category):
    '''returns if an epitope belongs to classically and alternatively defined neoepitopes (CDN vs ADN) (indicate which category to examine by category)--> Rech et al, 2018
    grouping is based on affinity and affinitiy foldchange between wt and mut
    '''
    group = "NA"
    if mhc == "mhcI":
        score_mut = props["best_affinity_netmhcpan4"]
        amplitude = props["Amplitude_mhcI_affinity"]
        bdg_cutoff_classical = 50
        bdg_cutoff_alternative = 5000
        amplitude_cutoff = 10
    elif mhc == "mhcII":
        score_mut = props["MHC_II_score_.best_prediction."]
        amplitude = props["Amplitude_mhcII"]
        bdg_cutoff_classical = 1
        bdg_cutoff_alternative = 4
        amplitude_cutoff = 4
    #print >> sys.stderr, score_mut, props["best_affinity_netmhcpan4_WT"], amplitude
    try:
        if category == "CDN":
            #print >> sys.stderr, float(score_mut), float(bdg_cutoff_classical)
            #print >> sys.stderr, float(score_mut) < float(bdg_cutoff_classical)
            if float(score_mut) < float(bdg_cutoff_classical):
                group = "True"
            elif float(score_mut) > float(bdg_cutoff_classical):
                group = "False"
        elif category == "ADN":
            if float(score_mut) < float(bdg_cutoff_alternative) and float(amplitude) > float(amplitude_cutoff):
                group = "True"
            elif float(score_mut) > float(bdg_cutoff_alternative) or float(amplitude) < float(amplitude_cutoff):
                group = "False"
    except ValueError:
        group = "NA"
    #print >> sys.stderr, category + ": "+group
    return group



if __name__ == '__main__':
    import sys
    basedir = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT"
    sys.path.append(basedir)
    from helpers import data_import
    import predict_all_epitopes
    import epitope


    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/20170713_IS_IM_data.complete.update_Dv10.csv"
    #file = "/scratch/info/projects/CM04_Core_NGS_Processing/MG148/A2DR1_GBMs_M5/scratch/scratch/A2DR1_GBMs_M5_mut_set.txt.transcript.squish.somatic.freq.annotation"
    file = "/scratch/info/projects/CM04_Core_NGS_Processing/MG148/A2DR1_GBMs_M7/scratch/scratch/A2DR1_GBMs_M7_mut_set.txt.transcript.squish.somatic.freq.annotation"
    db_uniprot = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/database/mouse_uniprot_swissprot_plus_isoforms.fasta"
    indel = False
    dat = data_import.import_dat_icam(file, indel)
    if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
            dat = data_import.change_col_names(dat)
    if "mutation_found_in_proteome" not in dat[0]:
            db = predict_all_epitopes.Bunchepitopes().build_proteome_dict(db_uniprot)

    dict_all = predict_all_epitopes.Bunchepitopes().Allepit
    ## add priority_score
    for ii,i in enumerate(dat[1]):
        #print dat[1][ii]
        # dict for each epitope
        z = epitope.Epitope()
        #print z
        z.init_properties(dat[0], dat[1][ii])
        z.add_features(number_of_mismatches(z.properties, "mhcI"), "Number_of_mismatches_mhcI")
        if "mutation_found_in_proteome" not in z.properties:
            z.add_features(match_in_proteome(z.properties, db), "mutation_found_in_proteome")
        z.add_features(calc_priority_score(z.properties), "Priority_score")
        #print z.properties

        for key in z.properties:
            if key not in dict_all:
                # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                dict_all[key] = [z.properties[key]]
            else:
                dict_all[key].append(z.properties[key])
    #print dict_all

    header = dat[0]
    header.append("Priority_score")
    #print header

    print "\t".join(header)
    for i in range(len(dict_all["mutation"])):
        z = []
        [z.append(dict_all[col][i]) for col in header]
        print "\t".join(z)

    #properties = {}
    #for nam,char in zip(dat[0], dat[1][1]):
    #    properties[nam] = char

    #print apppend_aa_index1(properties)
    #print apppend_aa_index2(properties)
