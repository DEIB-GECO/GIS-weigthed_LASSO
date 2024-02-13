import numpy as np
import pandas as pd
import os
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import  MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import make_scorer
from sklearn.utils import shuffle
from sklearn.metrics import recall_score, precision_score
import random
import json
import sys
sys.path.append("/home/mongardi/Metagene_repo/src")
from run_experiments import *
from utils.utils import * 
from get_wgis_scores import*

dire_kidney= '/home/mongardi/Metagene_repo/data/data_kidney'
dire_r = '/home/mongardi/Metagene_repo/results/gis_sensitibity_analysis/rounds'
dire_1 = '/home/mongardi/Metagene_repo/data/prior_knowledge/genes_and_ids_all_red.csv'
dire_2 = '/home/mongardi/Metagene_repo/data/prior_knowledge/gene_scores_norm_go.csv'

genes_set = load_txt('/home/mongardi/Metagene/Cancer/scripts/kidney_multicol/controlled_features.txt')

df = pd.read_csv(os.path.join(dire_kidney, 'Kidney_df_tr_coding_new.csv'))
df = shuffle(df, random_state=42)
df.set_index(df.columns[0], inplace=True)
#df.head()

X = df.drop(['is_healthy'], axis=1)
y = y = df['is_healthy']

rpm = X.div(X.sum(axis=1).values, axis=0) *1e6
rpm_log = np.log2(rpm + 1)

# training and test split
X_train, X_test,  y_train, y_test = train_test_split(rpm_log, y, random_state = 42,
                                        test_size = 0.2, stratify=y)


most_relevant_genes = ['PIK3C2G', 'MT1F', 'NCKAP1', 'PCBP1', 'METTL8']
genes_set = [x for x in genes_set if x not in most_relevant_genes]


#from joblib import parallel_backend

for j in range(len(most_relevant_genes)):

    keys = range(11)
    dct = {x: [] for x in keys}
    dct_f = {x: [] for x in [0, 1]}

    print(dct)
    for k in range(100):
        random.seed(int(k))
        print('--------------------',most_relevant_genes[j],'--------------------')
        # add noise to most relevant features
        X_rf = rpm_log[most_relevant_genes[j]].copy().to_frame()
        for i in range(10):
            name = most_relevant_genes[j] + '_' + str(i+1)
            # choose sample to add noise
            samples = X_rf.sample(n=100, random_state=i+k).index
            feature = X_rf[most_relevant_genes[j]].copy()
            np.random.seed(seed=i)#
            feature.loc[samples] +=  np.random.normal(0,0.01,len(samples))
            X_rf[name] = feature

        #print(X_rf.corr())
        #X_rf.shape
        others = most_relevant_genes.copy()
        others.remove(most_relevant_genes[j])
        print('------------------')
        rpm_log_new = pd.concat([rpm_log[genes_set],rpm_log[others], X_rf], axis=1)
        # training and test split
        X_train, X_test,  y_train, y_test = train_test_split([], rpm_log_new, y, random_state = 42,
                                                test_size = 0.2, stratify=y)

        penalties = get_GIS_scores_lasso_func([most_relevant_genes[j]], k=1, v=1, go=True, reactome=False,
                                       hpo=False, notebook=True, dire_1=dire_1, dire_2=dire_2)
        #print(penalties)
        print(np.min(penalties), np.max(penalties))
        #p_g = get_GIS_scores_lasso_func([most_relevant_genes[j]], k=1, v=1)
        #print(p_g)
        #penalties[-10:] = np.linspace(p_g[0],  np.max(penalties), 10)
        #penalties[-10:] = [0.95 for x in range(10)]
        penalties[-10:] = np.linspace(0.5,  np.max(penalties), 10)
        print(penalties[-10:])
        #print(X_train.shape, X_test.shape)
        classes , train_counts = np.unique(y_train, return_counts=True)
        train_percentages = train_counts / len(y_train)
        dict(zip(classes, train_percentages))

        store_results = {'gene_list': X_train.columns.values}

        # training and testing

        #lambdas = [2.0, 1.0, 0.5]
        lambdas = [5.0, 4.0, 2.0, 1.0]
        scoring_results = pd.DataFrame(np.zeros([13, 4]), index=X_train.columns[-13:])
        scoring_results.columns = lambdas

        print(scoring_results.head())

        for i in lambdas:

            split = 0
            scoring = {'accuracy':'accuracy',
                    'rec': make_scorer(recall_score, average= 'macro'),
                    'prec': make_scorer(precision_score, average= 'macro')}

            # rec_all': make_scorer(recall_score, average= None),
            # 'prec_all': make_scorer(recall_score, average=None)
            # 'prec': make_scorer(precision_score, pos_label='LumA')
            lasso_clf = Pipeline(steps=[('Scaler', MinMaxScaler()),
                                    ('Model', LogisticRegression(solver='saga',max_iter=10000, tol=1e-4, penalty='l1', random_state=split, penalty_factor=np.array(penalties)))])
            # penalty='l1'
            # 1.0, 5.0 (12), 2.0 (11), 4.0 (10)
            C=np.array([i])
            #C=np.array([0.5])
            #C=np.array([20.0])
            params = {'Model__C': C}
            cv = GridSearchCV(lasso_clf, params, cv=5, refit='accuracy', scoring=scoring)
            cv.fit(X_train, y_train)

            results = pd.DataFrame(cv.cv_results_)
            #results
            print("The best C paramter:",cv.best_params_)
            print("The CV accuracy of the corresponding model is:",cv.best_score_)


            print()
            preds = cv.predict(X_test)
            # X_train.columns.values
            betas = pd.DataFrame(cv.best_estimator_['Model'].coef_)
            # betas.shape
            # betas
            mask = betas!=0
            #cv.best_estimator_['Model'].intercept_
            print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")
            betas.columns = X_train.columns.values
            #betas.head()
            print(betas.iloc[0][-15:])
            scoring_results.loc[:, i] = betas.iloc[0][-13:].values
            features = betas.loc[:, (betas != 0).any(axis=0)]
            #features.to_csv('/home/mongardi/Metagene/Cancer/results/kidney/kidney_tr_'+ str(lam) +'.csv',index=False)
            print(betas.iloc[0][-11:].index)
            dct_f[0].append(sum(mask.sum(axis=0)!= 0))
            dct_f[1].append(mask.loc[0][-10:].sum().astype(np.float64))
            for f, ft in enumerate(betas.iloc[0][-11:].index):
                #print(ft)
                dct[f].append(betas.iloc[0][ft])
            # storing the mean for each gene across the different classes
            store_results['split'+ str(split)] = betas.mean(axis=0)
        #scoring_results.to_csv('/home/mongardi/Metagene/Cancer/scripts/kidney_multicol/' + most_relevant_genes[j] + '_gis_1_coeff.csv')
        # visualize resuts
    save_dict(dct, os.path.join(dire_r, 'rounds_'+ most_relevant_genes[j] + '_4'))
    save_dict(dct_f,os.path.join(dire_r,'rounds_features_'+ most_relevant_genes[j] + '_4'))