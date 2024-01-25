# -*- coding: utf-8 -*-
"""
build MNB model for PUS-dependent psi prediction.

reference: 
    1. https://www.kaggle.com/code/singhakash/dna-sequencing-with-machine-learning
    2. https://iq.opengenus.org/text-classification-naive-bayes/
"""

import os
import sys
import pandas as pd# pip install openpyxl
import matplotlib.pyplot as plt
import numpy as np
import pickle
from sklearn.feature_extraction.text import CountVectorizer
import argparse

# Define command line arguments
parser = argparse.ArgumentParser(description='predict PUS-targeting observation by loading psiMNB model.')
parser.add_argument('-model_file', metavar='FILE', type=str, help='model file to be loaded')
parser.add_argument('-to_predict', metavar='FILE', type=str, help='file to be predicted')
parser.add_argument('-output_dir', metavar='OUTPUT_DIR', type=str, help='output file name to be saved')

# Parse the arguments
args = parser.parse_args()

# Check that required arguments are provided
if not args.model_file:
    parser.error('the -model_file FILE argument for model file is required')
if not args.to_predict:
    parser.error('the -to_predict FILE argument for ready-to-be-predicted file is required')
if not args.output_dir:
    parser.error('the -output_dir FILE argument for output dir is required')

model_file = args.model_file
to_predict = args.to_predict
output_dir = args.output_dir

# Check if directory exists and is writable
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
elif not os.access(output_dir, os.W_OK):
    raise PermissionError("Output directory is not writable")

# Display help if no arguments provided
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

#convert our training data sequences into short overlapping k-mers of length 4. 
def Kmers_funct(seq, size=4):
    return [seq[x:x+size].lower() for x in range(len(seq) - size + 1)]

#compare observation
def compare(row):
  max_val = row.max()
  max_col = row.idxmax()
  return max_col

####preform prediction####
with open(model_file, 'rb') as f:
    cv, classifier = pickle.load(f)

to_pred = pd.read_csv(to_predict,sep="\t")
words = to_pred.apply(lambda x: Kmers_funct(x['Y_extendSeq_20nt']), axis=1)

for item in range(len(words)):
    words[item] = ' '.join(words[item])

to_pred_seq_cv = cv.transform(words)
predicted_label = classifier.predict(to_pred_seq_cv)
unique, counts = np.unique(predicted_label, return_counts=True)
print(dict(zip(unique, counts)))
predicted_proba = classifier.predict_proba(to_pred_seq_cv)

predicted_proba_df = pd.DataFrame(predicted_proba, columns=list(classifier.classes_))
predicted_proba_df['PUS-targeting'] = predicted_proba_df.apply(compare, axis=1)
to_pred = pd.concat([to_pred, predicted_proba_df], axis=1)
# to_pred.to_excel(os.path.join(output_dir, model_file+"_prediction.xlsx"),index=False)
to_pred.to_csv(os.path.join(output_dir, os.path.splitext(os.path.basename(model_file))[0] +"_prediction.txt"), sep="\t", index=False)


