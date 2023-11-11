# -*- coding: utf-8 -*-
"""
build MNB model for PUS-dependent psi prediction.

reference: 
    1. https://www.kaggle.com/code/singhakash/dna-sequencing-with-machine-learning
    2. https://iq.opengenus.org/text-classification-naive-bayes/
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

def test_kmer_alpha(size,alpha):
    path="G:/ModificationToolBox/PseudoTB_manuscript/psiFinder_20221012/motif_MNB"
    os.chdir(path)

    human_rna = pd.read_table('human_PUS_MNB_input_k-mer_overall.txt')
    human_rna.head()

    # human_rna['class'].value_counts().sort_index().plot.bar()
    # plt.title("Class distribution of Human PUS")


    ####build model####
    #convert our training data sequences into short overlapping k-mers of length 4. 
    def Kmers_funct(seq, size=size):
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
    cv = CountVectorizer(ngram_range=(1,1),min_df=5)#The n-gram size of 1 is previously determined by testing
    X_human = cv.fit_transform(human_texts)
    print(X_human.shape)

    # Splitting the human dataset into the training set and test set
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X_human, 
                                                        y_human, 
                                                        test_size = 0.20, 
                                                        random_state=42)

    from sklearn.naive_bayes import MultinomialNB
    classifier = MultinomialNB(alpha=alpha)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)

    from sklearn.naive_bayes import MultinomialNB
    import pickle
    classifier = MultinomialNB(alpha=alpha)
    classifier.fit(X_train, y_train)
    with open('PUS1_multinomialnb_model.pkl', 'wb') as f:
        pickle.dump(classifier, f)
    y_pred = classifier.predict(X_test)

    from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, precision_score, recall_score
    print("Confusion matrix for predictions on human test RNA sequence\n")
    print(pd.crosstab(pd.Series(y_test, name='Actual'), pd.Series(y_pred, name='Predicted')))
    def get_metrics(y_test, y_predicted):
        accuracy = accuracy_score(y_test, y_predicted)
        precision = precision_score(y_test, y_predicted, average='weighted')
        recall = recall_score(y_test, y_predicted, average='weighted')
        f1 = f1_score(y_test, y_predicted, average='weighted')
        mcc = matthews_corrcoef(y_test, y_predicted)
        return accuracy, precision, recall, f1, mcc
    accuracy, precision, recall, f1, mcc = get_metrics(y_test, y_pred)
    # print("accuracy = %.3f \nprecision = %.3f \nrecall = %.3f \nf1 = %.3f" % (accuracy, precision, recall, f1))
    # print("f1 = %.3f" % (f1))
    return round(f1,3)
    

kmer_df = pd.DataFrame(columns=['kmer', 'F1'])
for k in range(3, 10):
    print("\n","~~~~~~~~~~~~~~~~~~~~~~k-mers=",k)
    F1_tmp = test_kmer_alpha(size=k,alpha=0.1)
    kmer_df.loc[len(kmer_df)] ={'kmer': k, 'F1': F1_tmp}

plt.plot(kmer_df['kmer'], kmer_df['F1'])
plt.xlabel('k-mer size') 
plt.ylabel('F1 score')
kmer_df.index=kmer_df["kmer"]
kmer_df = kmer_df.drop(columns=["kmer"])
plt.text(8, 0.75,str(kmer_df), fontsize = 10)
plt.savefig('kmer_plot.png') 
plt.show()


alpha = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
alpha_df = pd.DataFrame(columns=['alpha', 'F1'])
for k in alpha:
    print("\n","~~~~~~~~~~~~~~~~~~~~~~alpha=",k)
    F1_tmp = test_kmer_alpha(size=4,alpha=k)
    alpha_df.loc[len(alpha_df)] ={'alpha': k, 'F1': F1_tmp}

plt.plot(alpha_df['alpha'], alpha_df['F1'])
plt.xlabel('alpha')
plt.ylabel('F1')
alpha_df.index=alpha_df["alpha"]
alpha_df = alpha_df.drop(columns=["alpha"])
plt.text(0.7, 0.85,str(alpha_df), fontsize = 10)
plt.savefig('alpha_plot.png') 
plt.show()