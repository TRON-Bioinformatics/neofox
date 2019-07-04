import numpy as np
import scipy.io as sio
import pandas as pd
import pickle
import sys
import subprocess

load_data = sio.loadmat('SIRdata.mat')

with open('genes-expression.pickle','rb') as handle:
    dict_expression=pickle.load(handle)


with open('amino-acids-features.pickle','rb') as handle:
    dict_data=pickle.load( handle)

def seq2bin(seq):
  
    aa = [ 'A' , 'C' , 'D' , 'E' , 'F' , 'G' , 'H' , 'I' , 'K' , 'L' , 'M' , 'N' , 'P' , 'Q' , 'R' , 'S' , 'T' , 'V' , 'W' , 'Y' ]
    dict_aa=dict((i,j) for j,i in enumerate(aa))
    arr=np.zeros((1,9*20))
    for ii,letter in enumerate(seq):
        arr[0,ii*20+dict_aa.get(letter)]=1

    return arr

def get_hydrophbicity(x,dict_):
    pair_letters=[c for c in x if c.isupper()]
    res=np.abs(dict_[pair_letters[0]]-dict_[pair_letters[1]])

    return res

def get_size(x,dict_):
    pair_letters = [c for c in x if c.isupper()]
    res = np.abs(dict_[pair_letters[0]] - dict_[pair_letters[1]])

    return res

def get_charge_change(x,dict_):
    pair_letters = [c for c in x if c.isupper()]
    if dict_[pair_letters[0]] ==dict_[pair_letters[1]]:
        return 0
    else:
        return 1

def get_charge_abs(x,dict_):
    pair_letters = [c for c in x if c.isupper()]
    if dict_[pair_letters[0]] ==dict_[pair_letters[1]]:
        return 0
    else:
        return 1


def get_polar(x,dict_):
    pair_letters = [c for c in x if c.isupper()]
    res=np.abs(dict_[pair_letters[0]]-dict_[pair_letters[1]])

    return res

def get_absolute(x,dict_):
    pair_letters=[c for c in x if c.isupper()]
    res=np.abs(dict_[pair_letters[0]]-dict_[pair_letters[1]])

    return res


def get_diffetenet(x,dict_):
    pair_letters = [c for c in x if c.isupper()]
    if dict_[pair_letters[0]] ==dict_[pair_letters[1]]:
        return 0
    else:
        return 1


def get_geneExpression(gene):
    res=dict_expression.get(gene)
    if res==None:
        res=np.nan

    return res


def get_properties(amino_substitation):
    lst_fetures=[]
    lst_fetures.append(get_diffetenet(amino_substitation,dict_data['Charge']))
    lst_fetures.append( get_absolute(amino_substitation,dict_data['Size']))
    lst_fetures.append(get_absolute(amino_substitation,dict_data['Hydro']))
    lst_fetures.append(get_absolute(amino_substitation,dict_data['Charge']))
    lst_fetures.append(get_diffetenet(amino_substitation,dict_data['Polar']))

    arr=np.asarray(lst_fetures)

    return arr


def isBind(seq):
    cutOff=load_data.get('CutOffs')

    bind_mat=load_data.get('PWMs')
    bnd_score=bind_mat.dot(seq.T).T
    b=np.zeros_like(bnd_score)
    b[bnd_score>cutOff.T]=1

    return bnd_score,b


def main(f_name):
    lst_data=[]
    with open(f_name,'r') as f:
        for row in f:
            gene_name,sequence,aa_subs=row.split()
            seq_arr = seq2bin(sequence)
            # tap score
            tap_mat = load_data.get('tap')
            tap_score = tap_mat.dot(seq_arr.T).ravel()
            #cleavge score
            clv_mat = load_data.get('clv')
            clv_mat = clv_mat[0, 20:200]
            clv_score = clv_mat.dot(seq_arr.T).ravel()
            #binding score
            bnd_score, cut_off = isBind(seq_arr)

            features_aa=get_properties(aa_subs)
            #expresion
            expression_value=get_geneExpression(gene_name)

            lst_data.append(np.hstack((expression_value,features_aa,clv_score,tap_score)))
            mat_features=np.asarray(lst_data)
    return mat_features


if __name__ == "__main__":
    main(sys.argv[1])





