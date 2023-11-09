# -*- coding: utf-8 -*-
"""
build MNB model for PUS-dependent psi prediction.

reference: https://www.kaggle.com/code/singhakash/dna-sequencing-with-machine-learning
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path="G:/ModificationToolBox/PseudoTB_manuscript/psiFinder_20221012/motif_MNB"
os.chdir(path)

human_rna = pd.read_table('human_PUS_MNB_input_k-mer_PUS1.txt')
human_rna.head()

human_rna['class'].value_counts().sort_index().plot.bar()
plt.title("Class distribution of Human PUS")


####build model####
# human_rna = human_rna[human_rna['class'] == 0]

#convert our training data sequences into short overlapping k-mers of length 4. 
def Kmers_funct(seq, size=4):
    return [seq[x:x+size].lower() for x in range(len(seq) - size + 1)]

human_rna['words'] = human_rna.apply(lambda x: Kmers_funct(x['Y_extendSeq_20nt']), axis=1)
human_rna = human_rna.drop('Y_extendSeq_20nt', axis=1)

human_rna.head()

human_texts = list(human_rna['words'])
for item in range(len(human_texts)):
    human_texts[item] = ' '.join(human_texts[item])
    
#separate labels
y_human = human_rna.iloc[:, 0].values # y_human for human_rna
y_human


from sklearn.feature_extraction.text import CountVectorizer
cv = CountVectorizer(ngram_range=(1,1),min_df=5)#The n-gram size of (x,x) is previously determined by testing
X_human = cv.fit_transform(human_texts)
print(X_human.shape)

# Splitting the human dataset into the training set and test set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X_human, 
                                                    y_human, 
                                                    test_size = 0.20, 
                                                    random_state=42)

from sklearn.naive_bayes import MultinomialNB
import pickle
classifier = MultinomialNB(alpha=0.1)
classifier.fit(X_train, y_train)
with open('PUS1_multinomialnb_model.pkl', 'wb') as f:
    pickle.dump(classifier, f)
y_pred = classifier.predict(X_test)

from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
print("Confusion matrix for predictions on human test RNA sequence\n")
print(pd.crosstab(pd.Series(y_test, name='Actual'), pd.Series(y_pred, name='Predicted')))
def get_metrics(y_test, y_predicted):
    accuracy = accuracy_score(y_test, y_predicted)
    precision = precision_score(y_test, y_predicted, average='weighted')
    recall = recall_score(y_test, y_predicted, average='weighted')
    f1 = f1_score(y_test, y_predicted, average='weighted')
    return accuracy, precision, recall, f1
accuracy, precision, recall, f1 = get_metrics(y_test, y_pred)
print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f" % (accuracy, precision, recall, f1))


####preform prediction####
with open('PUS1_multinomialnb_model.pkl', 'rb') as f:
    classifier = pickle.load(f)

to_pred = pd.read_table('Day0_common_anno_group_redundance_mix.txt')
to_pred.head()
words = to_pred.apply(lambda x: Kmers_funct(x['Y_extendSeq_20nt']), axis=1)

for item in range(len(words)):
    words[item] = ' '.join(words[item])

to_pred_seq_cv = cv.transform(words)
predicted_label = classifier.predict(to_pred_seq_cv)
unique, counts = np.unique(predicted_label, return_counts=True)
print(dict(zip(unique, counts)))
predicted_proba = classifier.predict_proba(to_pred_seq_cv)

predicted_proba_df = pd.DataFrame(predicted_proba, columns=['non-PUS1', 'PUS1'])
to_pred = pd.concat([to_pred, predicted_proba_df], axis=1)
to_pred['predicted_PUS1_label'] = pd.DataFrame(predicted_label)
mapping = {0: 'non-PUS1', 1: 'PUS1'}
to_pred['enzyme'] = to_pred['predicted_PUS1_label'].map(mapping)
to_pred.head()
to_pred.to_excel("Day0_common_anno_group_redundance_mix_PUS1_prediction.xlsx",index=False)# pip install openpyxl