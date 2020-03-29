import preprocess_path_added as preprocess
import sys
import pickle
import pandas as pd
import os

my_path = os.path.abspath(os.path.dirname(__file__))


def main(f_name):

    input_file = f_name
    mat = preprocess.main(input_file)
    with open(my_path + "/" + 'Classifier.pickle', 'rb') as f:
        clf = pickle.load(f)
    scores = clf.predict_proba(mat)
    dict_ = {}
    with open(input_file, 'r') as f:
        for row, val in zip(f, scores):
            seq = row.split()[1]
            dict_[seq] = val[-1]

    return dict_


if __name__ == "__main__":
    '''

    :params: first argument file to process, second argument name of new file
    :return: None
    '''
    dict_results = main(sys.argv[-2])
    df = pd.DataFrame.from_dict(dict_results, orient='index')
    df.to_csv(sys.argv[-1])
