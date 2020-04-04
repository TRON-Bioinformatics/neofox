from functools import reduce
#!/usr/bin/env python

def freq_aa(props, dict_freq):
    '''
    This function annotates data frame with frequency of mutated AA in the human proteome
    '''
    mut_aa = props["MUT_AA"]
    try:
        return str(dict_freq[mut_aa])
    except KeyError:
        return "NA"


def freq_prod_4mer(props, dict_freq):
    '''
    This function extracts 4 aa that are directed to TCR (pos 4 to pos 7 within epitope) and calculates the product of aa frequencies
    '''
    mut_aa = props["MUT_AA"]
    epitope = props["MHC_I_epitope_.best_prediction."]
    try:
        epi_4mer = list(epitope[3:7])
        epi_freqs = [float(dict_freq[aa]) for aa in epi_4mer]
        freq_prod = reduce(lambda x, y: x * y, epi_freqs)
        return str(freq_prod)
    except (TypeError, KeyError) as e:
        return "NA"


def freq_4mer(props, dict_freq):
    '''Returns the frequency of 4mer directed to TCR
    '''
    epitope = props["MHC_I_epitope_.best_prediction."]
    try:
        epi_4mer = epitope[3:7]
        return str(dict_freq[epi_4mer])
    except KeyError:
        return "NA"



def build_frequency_dict(frequency_file, freq_dict):
    """Loads file with information of frequeny of nmers
    """
    c = 0
    with open(frequency_file) as f:
        for line in f:
            c += 1
            w = line.strip("\n").split(";")
            if c==1:
                header = w
                continue
            else:
                freq_dict[w[0]] = w[1]
    return freq_dict




if __name__ == '__main__':
    import sys
    basedir = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2"
    sys.path.append(basedir)
    import data_import
    from input import FeatureLiterature

    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/20170713_IS_IM_data.complete.update_Dv10.csv"
    indel = False
    dat = data_import.import_dat_icam(file, indel)

    properties = {}
    for nam,char in zip(dat[0], dat[1][1]):
        properties[nam] = char
    # add MUT AA
    properties["MUT_AA"] = FeatureLiterature.wt_mut_aa(properties, "mut")

    freq_file1 = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/new_features/20181108_AA_freq_prot.csv"
    freq_file2 = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/new_features/20181108_4mer_freq.csv"
    freq_dict = {}
    freq_dict_4mer = {}
    freq_dict = build_frequency_dict(freq_file1, freq_dict )
    freq_dict_4mer = build_frequency_dict(freq_file2,freq_dict_4mer)

    print(freq_aa(properties, freq_dict))
    print(freq_prod_4mer(properties, freq_dict))
    print(freq_4mer(properties, freq_dict_4mer))
