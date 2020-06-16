from functools import reduce


# !/usr/bin/env python

def freq_aa(mutated_aminoacid, dict_freq):
    '''
    This function annotates data frame with frequency of mutated AA in the human proteome
    '''
    try:
        return str(dict_freq[mutated_aminoacid])
    except KeyError:
        return "NA"


def freq_prod_4mer(mutation, dict_freq):
    '''
    This function extracts 4 aa that are directed to TCR (pos 4 to pos 7 within epitope) and calculates the product of aa frequencies
    '''
    try:
        epi_4mer = list(mutation[3:7])
        epi_freqs = [float(dict_freq[aa]) for aa in epi_4mer]
        freq_prod = reduce(lambda x, y: x * y, epi_freqs)
        return str(freq_prod)
    except (TypeError, KeyError) as e:
        return "NA"


def freq_4mer(mutation, dict_freq):
    '''Returns the frequency of 4mer directed to TCR
    '''
    try:
        epi_4mer = mutation[3:7]
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
            if c == 1:
                header = w
                continue
            else:
                freq_dict[w[0]] = w[1]
    return freq_dict
