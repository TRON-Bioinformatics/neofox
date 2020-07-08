import os
import pickle

import pandas as pd

from input.Tcell_predictor.preprocess import Preprocessor


def main(f_name, output_file, references):
    input_file = f_name
    mat = Preprocessor(references=references).main(input_file)
    # NOTE: we do not put the Classifier.pickle in the references.py because it is code and not data what's in there
    # thus it belongs with the code
    with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'Classifier.pickle'), 'rb') as f:
        classifier = pickle.load(f)
    scores = classifier.predict_proba(mat)
    dictionary = {}
    with open(input_file, 'r') as f:
        for row, val in zip(f, scores):
            seq = row.split()[1]
            dictionary[seq] = val[-1]
    df = pd.DataFrame.from_dict(dictionary, orient='index')
    df.to_csv(output_file)
