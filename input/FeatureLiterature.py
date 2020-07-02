# Python version with Biopython
# . /etc/profile.d/modules.sh; module load software/python/python-2.7.9

'''
This script takes as input table from iCAM pipeline and calculates Literature feature of neoantigens
'''

# import modules
import math

from logzero import logger

from input import MHC_I, MHC_II
from input.IEDB_Immunogenicity import predict_immunogenicity_simple
from input.helpers import properties_manager


def calc_IEDB_immunogenicity(epitope, mhc_allele, mhc_score, affin_filtering=False):
    '''
    This function determines the IEDB immunogenicity score
    '''
    try:
        if affin_filtering:
            if float(mhc_score) < 500:
                return str(predict_immunogenicity_simple.predict_immunogenicity(
                    epitope, mhc_allele.replace("*", "").replace(":",
                                                                                                                    "")))
            else:
                return "NA"
        else:
            return str(predict_immunogenicity_simple.predict_immunogenicity(
                epitope, mhc_allele.replace("*", "").replace(":",
                                                                                                                "")))
    except ValueError:
        return "NA"


def dai(score_mutation, score_wild_type, affin_filtering=False):
    """
    Calculates DAI: Returns difference between wt and mut MHC binding score.
    """
    # TODO: these conversions to float need to go away from here
    try:
        if affin_filtering:
            if float(score_mutation) < 500:
                return str(float(score_wild_type) - float(score_mutation))
            else:
                return "NA"
        else:
            return str(float(score_wild_type) - float(score_mutation))
    except ValueError:
        return "NA"


def diff_number_binders(num_mutation, num_wild_type):
    """
    returns absolute difference of potential candidate epitopes between mutated and wt epitope
    """
    # TODO: this is not the absolute difference, just the difference :S
    try:
        return str(float(num_mutation) - float(num_wild_type))
    except ValueError:
        return "NA"


def ratio_number_binders(num_mutation, num_wild_type):
    """
    returns ratio of number of potential candidate epitopes between mutated and wt epitope. if no WT candidate epitopes, returns number of mutated candidate epitopes per mps
    """
    try:
        return str(float(num_mutation) / float(num_wild_type))
    except ZeroDivisionError:
        return str(num_wild_type)
    except ValueError:
        return "NA"


def rna_expression_mutation(transcript_expression, vaf_rna):
    """
    This function calculates the product of VAF in RNA and transcript expression
    to reflect the expression of the mutated transcript
    """
    try:
        return str(float(transcript_expression) * float(vaf_rna)) if float(vaf_rna) > 0 else "NA"
    except ValueError:
        return "NA"


def expression_mutation_tc(transcript_expression, tumor_content):
    """
    calculated expression of mutation corrected by tumour content
    """
    corrected_expression = "NA"
    if tumor_content != "NA":
        if tumor_content > 0.0:
            try:
                corrected_expression = str(float(transcript_expression) / tumor_content)
            except ValueError:
                pass
    return corrected_expression


def number_of_mismatches(epitope_wild_type, epitope_mutation):
    """
    This function calculates the number of mismatches between the wt and the mutated epitope
    """
    p1 = 0
    try:
        for i, aa in enumerate(epitope_mutation):
            if aa != epitope_wild_type[i]:
                p1 += 1
        return str(p1)
    except IndexError:
        return "NA"


def match_in_proteome(sequence, db):
    """
    This function checks if the mutated epitope has an exact match in a protein database (uniprot)
    Returns 0 if mutation is present in proteome and 1 if it not present
    """
    try:
        seq_in_db = [sequence in entry for entry in db]
        return "0" if any(seq_in_db) else "1"
    except:
        return "NA"


def calc_logistic_function(mhc_score):
    '''Calculates negative logistic function given mhc score
    '''
    try:
        return float(1 / (1 + math.exp(5 * (float(mhc_score) - 2))))
    except (OverflowError, ValueError) as e:
        return "NA"


def calc_priority_score(vaf_tumor, vaf_rna, transcript_expr, no_mismatch, score_mut, score_wt, mut_not_in_prot):
    """
    This function calculates the Priority Score using parameters for mhc I.
    """
    L_mut = calc_logistic_function(score_mut)
    L_wt = calc_logistic_function(score_wt)
    priority_score = "NA"
    try:
        if vaf_tumor not in ["-1", "NA"]:
            priority_score = (L_mut * float(vaf_tumor) * math.tanh(float(transcript_expr))) * (
                    float(mut_not_in_prot) * (1 - 2 ** (-float(no_mismatch)) * L_wt))
        else:
            priority_score = (L_mut * float(vaf_rna) * math.tanh(float(transcript_expr))) * (
                    float(mut_not_in_prot) * (1 - 2 ** (-float(no_mismatch)) * L_wt))
    except (TypeError, ValueError) as e:
        pass
    return str(priority_score)


def wt_mut_aa(substitution, mut):
    '''Returns wt and mut aa.
    '''
    try:
        if mut == "mut":
            return substitution[-1]
        elif mut == "wt":
            return substitution[0]
    except ValueError:
        return "NA"


def write_ouptut_to_file(epitope_data):
    '''
    This function prints output, semilicon separated --> csv file
    '''
    dat_new = epitope_data[1]
    head_new = epitope_data[0]
    print(";".join(head_new))
    for ii, i in enumerate(dat_new):
        print(";".join(i))


def classify_adn_cdn(score_mutation, amplitude, bdg_cutoff_classical, bdg_cutoff_alternative, amplitude_cutoff, category):
    """
    returns if an epitope belongs to classically and alternatively defined neoepitopes (CDN vs ADN) (indicate which category to examine by category)--> Rech et al, 2018
    grouping is based on affinity and affinitiy foldchange between wt and mut
    """
    group = "NA"
    try:
        if category == "CDN":
            if float(score_mutation) < float(bdg_cutoff_classical):
                group = "True"
            elif float(score_mutation) > float(bdg_cutoff_classical):
                group = "False"
        elif category == "ADN":
            if float(score_mutation) < float(bdg_cutoff_alternative) and float(amplitude) > float(amplitude_cutoff):
                group = "True"
            elif float(score_mutation) > float(bdg_cutoff_alternative) or float(amplitude) < float(amplitude_cutoff):
                group = "False"
    except ValueError:
        group = "NA"
    return group


# if __name__ == '__main__':
#     import sys
#
#     basedir = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT"
#     sys.path.append(basedir)
#     from input.helpers import data_import
#     from input import predict_all_epitopes, epitope
#
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/20170713_IS_IM_data.complete.update_Dv10.csv"
#     # file = "/scratch/info/projects/CM04_Core_NGS_Processing/MG148/A2DR1_GBMs_M5/scratch/scratch/A2DR1_GBMs_M5_mut_set.txt.transcript.squish.somatic.freq.annotation"
#     file = "/scratch/info/projects/CM04_Core_NGS_Processing/MG148/A2DR1_GBMs_M7/scratch/scratch/A2DR1_GBMs_M7_mut_set.txt.transcript.squish.somatic.freq.annotation"
#     db_uniprot = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/database/mouse_uniprot_swissprot_plus_isoforms.fasta"
#     indel = False
#     dat = data_import.import_dat_icam(file, indel)
#     if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
#         dat = data_import.change_col_names(dat)
#     if "mutation_found_in_proteome" not in dat[0]:
#         db = predict_all_epitopes.Bunchepitopes().load_proteome(db_uniprot)
#
#     dict_all = predict_all_epitopes.Bunchepitopes().Allepit
#     ## add priority_score
#     for ii, i in enumerate(dat[1]):
#         # dict for each epitope
#         z = epitope.Epitope()
#         z.init_properties(dat[0], dat[1][ii])
#         z.add_features(number_of_mismatches(z.properties, MHC_I), "Number_of_mismatches_mhcI")
#         if "mutation_found_in_proteome" not in z.properties:
#             z.add_features(match_in_proteome(z.properties, db), "mutation_found_in_proteome")
#         z.add_features(calc_priority_score(z.properties), "Priority_score")
#
#         for key in z.properties:
#             if key not in dict_all:
#                 # keys are are feautres; values: list of feature values associated with mutated peptide sequence
#                 dict_all[key] = [z.properties[key]]
#             else:
#                 dict_all[key].append(z.properties[key])
#
#     header = dat[0]
#     header.append("Priority_score")
#
#     print("\t".join(header))
#     for i in range(len(dict_all["mutation"])):
#         z = []
#         [z.append(dict_all[col][i]) for col in header]
#         print("\t".join(z))
