import numpy as np
import pandas as pd
import os
import warnings
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import make_scorer
from sklearn.utils import shuffle
from sklearn.metrics import recall_score, accuracy_score, precision_score, f1_score
from joblib import parallel_backend
from get_bio_info import * 
from utils.utils import *
from get_wgis_scores import get_GIS_scores_lasso_func


def run_experiments(args, X, y, results_dire,tol =1e-4, lambdas=0.5, n_splits = 2, train_sample = True,
                    max_iter=10000,
                    cv_lambda = False, 
                    is_GIS= True, bio_info=True, 
                    is_train= True, k=1, v=1, go=True, reactome=False, hpo=False):
    

    # create results dire
    if not os.path.exists(results_dire):
        os.makedirs(results_dire)
    
    else:
        warnings.warn('Classifier may have already been built! Classifier results will be overwritten!', category=Warning)
  
    print(results_dire)
    n_classes = len(np.unique(y))
    print(X.shape)
    if is_train:
        # training and test split
        X_train, X_test,  y_train, y_test = train_test_split(X, y, random_state = 42,
                                            test_size = 0.2, stratify=y)
        
        
        if  is_GIS:
            penalties = get_GIS_scores_lasso_func(args, X.columns, k=k, v=v, go=go, reactome=reactome, hpo=hpo)
        else: 
            penalties = np.ones(len(X.columns))

        
        #n_splits = 2
        if train_sample: 
            training_results, union_genes, intersection_genes =  train(X_train, X_test, y_train, y_test, penalties, results_dire, tol=tol, lambdas = lambdas, n_splits=n_splits, cv_samples = True, cv_lambda=False, store_cv= False, max_iter=max_iter)
        
        if cv_lambda: 
            cv_results_df = train(X_train, X_test, y_train, y_test, penalties, results_dire, tol=tol, lambdas = lambdas, cv_samples = False, cv_lambda=True, store_cv= True,  max_iter=max_iter)
            print(cv_results_df)

    genes_list = X.columns

    if bio_info:
        union_genes, intersection_genes = get_union_inter(results_dire, X, n_splits)
        selected_genes = get_selected_genes(results_dire, intersection_genes, union_genes, n_splits, n_classes)
        #ann_matrix = get_bio_stats(genes_list, selected_genes, results_dire, labels=np.unique(y))
        get_bio_stats(args, genes_list, selected_genes, results_dire, labels=np.unique(y))
    #training_results
        #return ann_matrix 

    


def train(X, X_t, y, y_t, penalties, results_dire, tol=1e-4, lambdas = 0.0, n_splits=10, cv_samples = True, cv_lambda=False, store_cv= False, max_iter=10000):

    scoring = {'accuracy':'accuracy',
            'rec': make_scorer(recall_score, average= 'macro'),
            'prec': make_scorer(precision_score, average= 'macro'), 
            'f1':make_scorer(f1_score, average="macro")}

  
  
    classes = len(np.unique(y))
    if cv_samples:
        if isinstance(n_splits, list):
            splits = n_splits
        else:
            splits = np.arange(n_splits)

        store_results = {'gene_list': X.columns.values}
        accuracy = {'Cross-val':[], 'Test':[]}
        recall = {'Cross-val':[], 'Test':[]}
        precision = {'Cross-val':[], 'Test':[]}
        f1 = {'Cross-val':[], 'Test':[]}
        #from joblib import parallel_backend
        for split in splits: 
            with parallel_backend('threading', n_jobs=5):
                print('split: ', split)
                print('-------------------Training-------------------')
                lasso_clf = Pipeline(steps=[('Scaler', MinMaxScaler()),
                                ('Model', LogisticRegression(solver='saga',max_iter=max_iter, tol=tol, penalty='l1', n_jobs=-1,  penalty_factor= np.array(penalties),random_state=split, class_weight='balanced'))])

                C=np.array([lambdas])
                params = {'Model__C': C}
                cv = GridSearchCV(lasso_clf, params, cv=5, refit='accuracy', scoring=scoring)
            
                cv.fit(X, y)
                results = pd.DataFrame(cv.cv_results_)
                #print(results)
                print("The best C paramter:",cv.best_params_)
                print("The CV accuracy of the corresponding model is:",cv.best_score_)

                accuracy['Cross-val'].append(cv.best_score_)
                recall['Cross-val'].append(results['mean_test_rec'].at[0])
                precision['Cross-val'].append(results['mean_test_prec'].at[0])
                f1['Cross-val'].append(results['mean_test_f1'].at[0])

                

                print('-------------------Testing-------------------')

                preds = cv.predict(X_t)
             
                prediction_results = show_single_class_evaluation(preds, y_t, np.unique(y), results_dire)
                if split ==0:
                    prediction_results_all = prediction_results.copy()
                else:
                    prediction_results_all = pd.concat([prediction_results_all,prediction_results], axis = 0)
               
                
          
                betas = pd.DataFrame(cv.best_estimator_['Model'].coef_)
                mask = betas!=0
                #cv.best_estimator_['Model'].intercept_
                print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")
                betas.columns = X.columns.values
                #betas.head()
                features = betas.loc[:, (betas != 0).any(axis=0)]
                features.to_csv(results_dire + '/selected_genes_'+ str(split) +'.csv',index=False)
                # storing the mean for each gene across the different ]lasses
                store_results['split'+ str(split)] = betas.mean(axis=0)
        
        print('------------------- Results-------------------')

        print('Cross val accuracy:',np.mean(accuracy['Cross-val']), 'sd: ', (np.std(accuracy['Cross-val'])))
        print('Cross val recall:', np.mean(recall['Cross-val']), 'sd: ',(np.std(recall['Cross-val'])))
        print('Cross val precision:', np.mean(precision['Cross-val']), 'sd: ',(np.std(precision['Cross-val'])))

        print('----- Test Results----')
        print(prediction_results_all.mean(axis=0))
        print(prediction_results_all.std(axis=0))
            

        # analyse results
        store_results_df = pd.DataFrame(store_results)
        store_results_df = store_results_df.set_index('gene_list')
        if store_results_df.shape[0] > 1:
            union_genes = store_results_df.loc[( store_results_df != 0).any(axis=1), :].index
            intersection_genes = store_results_df.loc[( store_results_df != 0).all(axis=1), :].index
        else: 
            intersection_genes = store_results_df.loc[0][store_results_df.loc[0] != 0].index
            union_genes = intersection_genes

        print('Intersection:', len(intersection_genes))
        print('Union:', len(union_genes))
        return store_results_df, union_genes, intersection_genes


    if cv_lambda:

        if store_cv:

            cv_results_df = pd.DataFrame(np.zeros((len(lambdas), len(scoring))), columns = list(scoring.keys()), index=lambdas)
            for l in lambdas:

                #from joblib import parallel_backend
                print('-------------------Training-------------------')
                with parallel_backend('threading', n_jobs=5):
                    split = 0
                    lasso_clf = Pipeline(steps=[('Scaler', MinMaxScaler()),
                                    ('Model', LogisticRegression(solver='saga',max_iter=10000, tol=tol, penalty='l1', n_jobs=-1, penalty_factor= np.array(penalties),random_state=split, class_weight='balanced'))])

                    C=np.array([l])
                    params = {'Model__C': C}
                    cv = GridSearchCV(lasso_clf, params, cv=5, refit='accuracy', scoring=scoring)
                    cv.fit(X, y)
                    results = pd.DataFrame(cv.cv_results_)

                    print("The best C paramter:",cv.best_params_)
                    print("The CV accuracy of the corresponding model is:",cv.best_score_)

                    betas = pd.DataFrame(cv.best_estimator_['Model'].coef_)
                    mask = betas!=0
                    print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")    
                    cv_results_df.loc[l] = [cv.best_score_, results['mean_test_rec'].at[0], results['mean_test_prec'].at[0], results['mean_test_f1'].at[0]]

                 

            return cv_results_df
                   


        else: 
            ## to do
            print('-------------------Training-------------------')
            with parallel_backend('threading', n_jobs=5):
                split = 0
                lasso_clf = Pipeline(steps=[('Scaler', MinMaxScaler()),
                                ('Model', LogisticRegression(solver='saga',max_iter=max_iter, tol=tol, penalty='l1', n_jobs=-1, penalty_factor= np.array(penalties),random_state=split))])

                C=np.array([lambdas])
                params = {'Model__C': C}
                cv = GridSearchCV(lasso_clf, params, cv=5, refit=False, scoring=scoring)
                cv.fit(X, y)
                results = pd.DataFrame(cv.cv_results_)
                #print(results)
                print("The best C paramter:",cv.best_params_)
                print("The CV accuracy of the corresponding model is:",cv.best_score_)

                accuracy['Cross-val'].append(cv.best_score_)
                recall['Cross-val'].append(results['mean_test_rec'].at[0])
                precision['Cross-val'].append(results['mean_test_prec'].at[0])
                recall_all['Cross-val'].append(None)
                prec_all['Cross-val'].append(None)

                print('-------------------Testing-------------------')

                preds = cv.predict(X_t)
                accuracy['Test'].append(accuracy_score(y_t, preds))
                recall['Test'].append(recall_score(y_t, preds, average='macro'))
                precision['Test'].append(precision_score(y_t, preds, average='macro'))
                recall_all['Test'].append(recall_score(y_t, preds, average=None))
                prec_all['Test'].append(precision_score(y_t, preds, average=None))
              

                betas = pd.DataFrame(cv.best_estimator_['Model'].coef_)
                mask = betas!=0
                cv.best_estimator_['Model'].intercept_
                print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")
                betas.columns = X.columns.values
                betas.head()
                features = betas.loc[:, (betas != 0).any(axis=0)]
                features.to_csv(results_dire + '/selected_genes_cv' +'.csv',index=False)

                # storing the mean for each gene across the different classes
                store_results['split'+ str(split)] = betas.mean(axis=0)


       
