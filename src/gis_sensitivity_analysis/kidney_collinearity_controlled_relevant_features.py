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
dire = '/home/mongardi/GIS_prior_knowledge'

sys.path.append(os.path.join(dire, "src"))
from run_experiments import *
from utils.utils import * 
from get_wgis_scores import*

dire_kidney = os.path.join(dire,'/data/data_kidney')
dire_r =  os.path.join(dire,'results/gis_sensitibity_analysis')
dire_1 = os.path.join(dire,'/data/prior_knowledge/genes_and_ids_all_red.csv')
dire_2 = os.path.join(dire,'/data/prior_knowledge/gene_scores_norm_go.csv')

genes_set = load_txt(os.path.join(dire_kidney, 'controlled_features.txt'))

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


def load_txt(filename):
    content = []
    with open(filename)as f:
        for line in f:
            content.append(line.strip())
    return content


def save_dict(dictionary,filename):
    with open(filename + ".txt", "w") as fp:
        json.dump(dictionary, fp)
        print("Done writing dict into .txt file")

most_relevant_genes = ['PIK3C2G', 'MT1F', 'NCKAP1', 'PCBP1', 'METTL8']
genes_set = [x for x in genes_set if x not in most_relevant_genes]


fi_scores = fisher_scores(X_train, y_train.to_numpy())

rpm_log_new = pd.concat([rpm_log[genes_set],rpm_log[most_relevant_genes]],axis=1)
# training and test split
X_train, X_test,  y_train, y_test = train_test_split(rpm_log_new, y, random_state = 42,
                                        test_size = 0.2, stratify=y)

print(X_train.shape, X_test.shape)
classes , train_counts = np.unique(y_train, return_counts=True)
train_percentages = train_counts / len(y_train)
dict(zip(classes, train_percentages))

store_results = {'gene_list': X_train.columns.values}
# training and testing  
lambdas = [5.0, 4.0, 2.0, 1.0]
genes = list(fi_scores.iloc[0][most_relevant_genes].sort_values().index)
for i in lambdas:
    
    #storing_results = pd.DataFrame(np.zeros([50, 6]))
    #storing_results.columns = ['Features_s', 'Standard', 'Features_g', 'GIS', 'Features_gm', 'GIS modi']
    split = 0

    print('--------------- Standard LASSO-----------------')
    scoring = {'accuracy':'accuracy',
            'rec': make_scorer(recall_score, average= 'macro'),
            'prec': make_scorer(precision_score, average= 'macro')}

    # rec_all': make_scorer(recall_score, average= None),
    # 'prec_all': make_scorer(recall_score, average=None)
    # 'prec': make_scorer(precision_score, pos_label='LumA')
    lasso_clf = Pipeline(steps=[('Scaler', MinMaxScaler()),
                            ('Model', LogisticRegression(solver='saga',max_iter=10000, tol=1e-4, penalty='l1', random_state=split))])
    # penalty='l1'
    # 1.0, 5.0 (12), 2.0 (11), 4.0 (10)
    C=np.array([i])
    #C=np.array([0.5])
    #C=np.array([20.0])
    params = {'Model__C': C}
    cv = GridSearchCV(lasso_clf, params, cv=5, refit='accuracy', scoring=scoring)
    cv.fit(X_train, y_train)

    results = pd.DataFrame(cv.cv_results_)
    results
    print("The best C paramter:",cv.best_params_)
    print("The CV accuracy of the corresponding model is:",cv.best_score_)

    print()
    preds = cv.predict(X_test)
    # X_train.columns.values
    betas = pd.DataFrame(cv.best_estimator_['Model'].coef_)
    # betas.shape
    # betas
    mask = betas!=0
    cv.best_estimator_['Model'].intercept_
    print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")
    betas.columns = X_train.columns.values
    betas.head()
    #print(betas.iloc[0][-15:])
    features = betas.loc[:, (betas != 0).any(axis=0)]
 
    relevant_features = genes
    dct_features = {x: [] for x in relevant_features}
    dct_accuracy = {x: [] for x in relevant_features} 
    dct_n_features = {x: [] for x in relevant_features}

    for feat in relevant_features:
        dct_features[feat].append(betas.iloc[0][feat])
        dct_accuracy[feat].append(accuracy_score(y_test, preds))
        dct_n_features[feat].append(sum(mask.sum(axis=0)!= 0))



    for f in relevant_features:
        penalties = get_GIS_scores_lasso_func([], rpm_log_new.columns, k=1, v=1,  go=True, reactome=False,
                                       hpo=False, notebook=True, dire_1=dire_1, dire_2=dire_2)

        lasso_clf = Pipeline(steps=[('Scaler', MinMaxScaler()),
                            ('Model', LogisticRegression(solver='saga',max_iter=10000, tol=1e-4, penalty='l1', random_state=split, penalty_factor=np.array(penalties)))])
       

        params = {'Model__C': C}
        cv = GridSearchCV(lasso_clf, params, cv=5, refit='accuracy', scoring=scoring)
        cv.fit(X_train, y_train)

        results = pd.DataFrame(cv.cv_results_)
        results
        print("The best C paramter:",cv.best_params_)
        print("The CV accuracy of the corresponding model is:",cv.best_score_)
        
        print()
        preds = cv.predict(X_test)

        #X_train.columns.values
        betas = pd.DataFrame(cv.best_estimator_['Model'].coef_)
        betas
        mask = betas!=0
        cv.best_estimator_['Model'].intercept_
        print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")
        betas.columns = X_train.columns.values
        betas.head()
        features = betas.loc[:, (betas != 0).any(axis=0)]

        dct_features[f].append(betas.iloc[0][f])
        dct_accuracy[f].append(accuracy_score(y_test, preds))
        dct_n_features[f].append(sum(mask.sum(axis=0)!= 0))
   
        #values = np.linspace(np.min(penalties), 1.0, 10)
        values = np.linspace(0.5, 1.0, 50)
        for value in values:
            indexes = np.where(np.isin(X_train.columns,f))[0].tolist()
            print(indexes)
            #penalties_new = [value if k in indexes else penalties[k] for k, x in enumerate(penalties)]
            
            penalties_new  = penalties.copy()
            print(penalties_new[indexes[0]])
            penalties_new[indexes[0]] = value
            print(penalties_new[indexes[0]])
         
            lasso_clf = Pipeline(steps=[('Scaler', MinMaxScaler()),
                            ('Model', LogisticRegression(solver='saga',max_iter=10000, tol=1e-4, penalty='l1', random_state=split, penalty_factor=np.array(penalties_new)))])
            params = {'Model__C': C}
            cv = GridSearchCV(lasso_clf, params, cv=5, refit='accuracy', scoring=scoring)
            cv.fit(X_train, y_train)
            results = pd.DataFrame(cv.cv_results_)
            results
            print("The best C paramter:",cv.best_params_)
            print("The CV accuracy of the corresponding model is:",cv.best_score_)
            preds = cv.predict(X_test)
            #X_train.columns.values
            betas = pd.DataFrame(cv.best_estimator_['Model'].coef_)
            # betas
            # mask = betas!=0
            cv.best_estimator_['Model'].intercept_
            print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")
            betas.columns = X_train.columns.values
            betas.head()
            features = betas.loc[:, (betas != 0).any(axis=0)]
      
            dct_features[f].append(betas.iloc[0][f])
            dct_accuracy[f].append(accuracy_score(y_test, preds))
            dct_n_features[f].append(sum(mask.sum(axis=0)!= 0))
    save_dict(dct_features,'rounds_mrf'+ str(i))
    save_dict(dct_accuracy,'rounds_acc_mrf'+ str(i))
    save_dict(dct_n_features,'rounds_n_mrf'+ str(i))
