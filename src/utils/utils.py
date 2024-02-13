import numpy as np
import pandas as pd
import os
import json

from sklearn.metrics import recall_score, accuracy_score, precision_score,  balanced_accuracy_score, f1_score
from sklearn.feature_selection import mutual_info_classif
from scipy.stats import f_oneway, ttest_rel, ttest_ind, describe,  wilcoxon
from utils.fisher_score import *
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportions_ztest

def load_dictionary(filename):
        
        data_dict = open(filename)
        data_dict = json.load(data_dict)
        return data_dict 


def save_dictionary(dictionary,filename):
    with open(filename + ".json", "w") as fp:
        json.dump(dictionary, fp)
        print("Done writing dict into .json file")

def load_dataset(filename, idx=0):

    df = pd.read_csv(filename)
    #print(df.head())
    df_ind = df.set_index(df.columns[idx])
    #print(df_ind.head())
    return df_ind


def save_list(list_file, file_name):
    with open(file_name, 'w') as f:
        for s in list_file:
            f.write(str(s) + '\n')

def show_single_class_evaluation(y_pred: int, y_test: int, labels, results_dire, verbose=False):

    if verbose:
        print("Balanced accuracy: ", round(balanced_accuracy_score(y_test, y_pred), 5)) # not possible for single class
        print("Accuracy: ", round(accuracy_score(y_test, y_pred), 5)) # not possible for single class
        print('precision ', round(precision_score(y_test, y_pred, average="macro"), 5))
        print('recall ', round(recall_score(y_test, y_pred, average="macro"), 5))
        print('f1_macro ', round(f1_score(y_test, y_pred, average="macro"),5))
        print('f1_weighted ', round(f1_score(y_test, y_pred, average="weighted"),5))
        print("Precision: ", [round(i, 5) for i in precision_score(y_test, y_pred, average=None) ])
        print("Recall: ",  [round(i, 5) for i in recall_score(y_test, y_pred, average=None) ]) 
        print("F1 Score: ", [round(i, 5) for i in f1_score(y_test, y_pred, average=None) ]) 
        print('--------------------------------------------')

    dic_result = {}
    dic_result['balanced_accuracy'] = [round(balanced_accuracy_score(y_test, y_pred), 5)]
    dic_result['accuracy'] = [round(accuracy_score(y_test, y_pred), 5)]
    dic_result['precision'] = [round(precision_score(y_test, y_pred, average="macro"), 5)]
    dic_result['recall'] = [round(recall_score(y_test, y_pred, average="macro"), 5)]
    dic_result['f1_macro'] = [round(f1_score(y_test, y_pred, average="macro"),5)]
    dic_result['f1_weighted'] = [round(f1_score(y_test, y_pred, average="weighted"),5)]
    for i in range(len(labels)):
        dic_result[labels[i]+'-precision'] =  round( precision_score(y_test, y_pred, average=None)[i], 5)
    for i in range(len(labels)):
        dic_result[labels[i]+'-recall'] =  round( recall_score(y_test, y_pred, average=None)[i], 5)
    for i in range(len(labels)):   
        dic_result[labels[i]+'-f1_score'] =  round( f1_score(y_test, y_pred, average=None)[i], 5)

    df_result = pd.DataFrame.from_dict(dic_result)
    #df_result.to_csv(results_dire + '/output_detailed_scores.csv', index=False)
    return df_result


def get_stats(v):
    res = describe(v)
    return np.concatenate([
        [
            res.minmax[0],
            res.minmax[1],
            res.mean,
            res.variance,
            res.skewness,
            res.kurtosis
        ],
        np.percentile(v, q=[10, 25, 50, 75, 90])
    ])


# Predictive Power Scores and Functions

def get_subset(df):

    genes_coeffs = abs(df[df !=0].loc[0]).sort_values(ascending=False)
    genes = genes_coeffs.index.to_list()
    return genes

def mutual_information(X, y, df=None):

    #mutual information
    mi_scores = []
    if True in np.unique(y):
        levels = [True, False]
    else:
        levels = np.unique(y)
    print(levels)
    for level in levels:
        mi_class= mutual_info_classif(X, y==level)
        mi_scores.append(mi_class)
        #print(mi_class)
    mi_scores_df = pd.DataFrame(mi_scores, columns=X.columns)
    ratios = []

    if df is not None:

        for i in range(mi_scores_df.shape[0]):

            if levels[i]== True:
                gene_subset = df[df > 0].columns
            elif levels[i]== False:
                gene_subset = df[df < 0].columns
            else:
                print('to do')
            

            norm_rank  = mi_scores_df.iloc[i].rank() / np.max(mi_scores_df.iloc[i].rank())
            ratios.append([ norm_rank[gene_subset].sum()/ len(gene_subset),  mi_scores_df.iloc[i][gene_subset].sum() / len(gene_subset)])

        return ratios
    else:   
        return mi_scores_df

def get_feature_importance(X, y, df):

    gene_subset = get_subset(df)
    f_scores = fisher_scores(X, y, gene_subset)
    mi = mutual_information(X, y, df)
    return f_scores, mi
    
def one_way_anova(X, y, gene_subset):

    for gene in X.columns:

        classes = get_classes(gene, y)
        anova = f_oneway(classes)

    return anova

def get_classes(X, y):
    s = np.argsort(y)
    return np.split(X[s], np.unique(y[s], return_index=True)[1][1:])

def fisher_scores(X, y, gene_subset=None):

    scores = fisher_score(X,y)
    scores_df = pd.DataFrame(scores.tolist(), index=X.columns).T

    if gene_subset is not None:
        norm_rank  = scores_df.iloc[0].rank(ascending=False) / np.max(scores_df.iloc[0].rank())

        return norm_rank[gene_subset].sum()/ len(gene_subset), scores_df.iloc[0][gene_subset].sum() / len(gene_subset)
    
    else:
        return scores_df
    

###--------

def get_union_inter(results_dire, X, n_splits):

    files = os.listdir(results_dire) 
    files = [x for x in files if x.startswith('selected')]
    files.sort()
    #print(files)
    if isinstance(n_splits, list):
            n_splits = len(n_splits)

    ff = pd.DataFrame(np.zeros((n_splits, len(X.columns))), columns = X.columns )
    
    for i, file in enumerate(files[:n_splits]):
        #print(file)
        features_df = pd.read_csv(os.path.join(results_dire, file))
        mask = features_df!=0
        #print(sum(mask.sum(axis=0)!= 0), "features were used in this model plus the intercept")
        
        features = features_df.loc[:, (features_df != 0).any(axis=0)]
        #features.shape
        #print(features)
        #print(features.columns.values)
        ff.iloc[i][features.columns.values] = 1.0

    union_genes = ff.loc[:, (ff != 0).any(axis=0)].columns
    print(len(union_genes))
    intersection_genes = ff.loc[:, (ff != 0).all(axis=0)].columns
    print(len(intersection_genes))

    return union_genes, intersection_genes

def get_selected_genes(results_dire, intersection_genes, union_genes, n_splits, n_classes, verbose=False):

    files = os.listdir(results_dire) 
    files = [x for x in files if x.startswith('selected')]
    files.sort()
    if n_classes <=2:
        n_classes = 1
        if verbose:
            print(files)
    ff_union = pd.DataFrame(np.zeros((n_splits,len(union_genes))), columns = union_genes)
    ff_intersec = pd.DataFrame(np.zeros((n_splits,len(intersection_genes))), columns = intersection_genes)
    df_intersec = pd.DataFrame(np.zeros((n_classes,len(intersection_genes))), columns = intersection_genes)
    
    if isinstance(n_splits, list):
            n_splits = len(n_splits)

    for k in range(n_classes):
        
        for i, file in enumerate(files[:n_splits]):
    

            features_df = pd.read_csv(os.path.join(results_dire,file))
            if verbose:
                print(len(features_df.loc[k][features_df.loc[k]!=0]))
            features = features_df.loc[k].index
            ff_union.iloc[i][features] = features_df.iloc[k]
            features = [f for f in features if f in intersection_genes]
            if verbose:
                print(len(features))
            ff_intersec.iloc[i][features] = features_df.iloc[k][features]


        #ff_intersec_mean =   abs(ff_intersec.mean(axis=0)).sort_values(ascending=False)
        #ff_intersec_mean =   abs(ff_intersec.mean(axis=0))
        ff_intersec_mean =   abs(ff_intersec).min(axis=0)
        if verbose:
            print('--------------------------')
            print(len(ff_intersec_mean[ff_intersec_mean != 0]))
            print('--------------------------')
        df_intersec.iloc[k] = ff_intersec_mean

    return df_intersec


###--------

def get_info_enriched_terms(dire1, dire2, ratio=True):

      # dire1
      files_dire1 = os.listdir(dire1)
      files_dire1 = [x for x in files_dire1 if x.startswith('significant') and 'info' in x]
      files_dire1.sort()
      print('File 1')
      results_dire_1_df = view_enrichment_results(dire1)
      # dire2
      files_dire2 = os.listdir(dire2)
      files_dire2 = [x for x in files_dire2 if x.startswith('significant') and 'info' in x]
      files_dire2.sort()
      print('File 2')
      results_dire_2_df = view_enrichment_results(dire2)
      classes = results_dire_1_df['Subtype'].values
      print(classes)

    

      for i in range(len(files_dire1)):

            #class_name =  files_dire1[i].split('_')[-1].split('.')[0]
            class_name = classes[i]
            print()
            print('---------------------', class_name, '---------------------')
            print()

            dict_dire1 = load_dictionary(os.path.join(dire1, files_dire1[i]))
            dict_dire2 = load_dictionary(os.path.join(dire2, files_dire2[i]))   

            print('T-test Proportions')
            print('----GO----')
            count = np.array([results_dire_2_df.iloc[i]['GO sig'],results_dire_1_df.iloc[i]['GO sig'] ])
            nobs = np.array([results_dire_2_df.iloc[i]['GO'],results_dire_1_df.iloc[i]['GO'] ])
            stat, pval = proportions_ztest(count, nobs, alternative='larger')
            print('{0:0.3f}'.format(pval))
            print('----KEGG----')
            count = np.array([results_dire_2_df.iloc[i]['KEGG sig'],results_dire_1_df.iloc[i]['KEGG sig'] ])
            nobs = np.array([results_dire_2_df.iloc[i]['KEGG'],results_dire_1_df.iloc[i]['KEGG'] ])
            stat, pval = proportions_ztest(count, nobs, alternative='larger')
            print('{0:0.3f}'.format(pval))
            print('----Reactome----')
            count = np.array([results_dire_2_df.iloc[i]['REACTOME sig'],results_dire_1_df.iloc[i]['REACTOME sig'] ])
            nobs = np.array([results_dire_2_df.iloc[i]['REACTOME'],results_dire_1_df.iloc[i]['REACTOME'] ])
            stat, pval = proportions_ztest(count, nobs, alternative='larger')
            print('{0:0.3f}'.format(pval))
            print('----HPO----')
            count = np.array([results_dire_2_df.iloc[i]['HPO sig'],results_dire_1_df.iloc[i]['HPO sig'] ])
            nobs = np.array([results_dire_2_df.iloc[i]['HPO'],results_dire_1_df.iloc[i]['HPO'] ])
            stat, pval = proportions_ztest(count, nobs, alternative='larger')
            print('{0:0.3f}'.format(pval))

            if ratio:
                print('----GO----')
                rat = results_dire_2_df.iloc[i]['GO'] / results_dire_1_df.iloc[i]['GO'] 
                print('{0:0.3f}'.format(rat))
                print('----KEGG----')
                rat = results_dire_2_df.iloc[i]['KEGG'] / results_dire_1_df.iloc[i]['KEGG'] 
                print('{0:0.3f}'.format(rat))
                print('----Reactome----')
                rat = results_dire_2_df.iloc[i]['REACTOME'] / results_dire_1_df.iloc[i]['REACTOME'] 
                print('{0:0.3f}'.format(rat))
                print('----HPO----')
                rat = results_dire_2_df.iloc[i]['HPO'] / results_dire_1_df.iloc[i]['HPO'] 
                print('{0:0.3f}'.format(rat))
                
            if len(dict_dire1) > 0 and len(dict_dire2)>0:
                  data = [list(dict_dire1.values()), list(dict_dire2.values())]

                  print('--------Statistics--------')
                  print('File 1')
                  stats1 =  get_stats(data[0])[:4]
                  print('Min:', stats1[0],
                        'Max:', stats1[1], 
                        'Mean:', stats1[2], 
                        'std:', stats1[3])
                  print('File 2')
                  stats2 = get_stats(data[1])[:4]
                  print('Min:', stats2[0], 
                        'Max:', stats2[1], 
                        'Mean:', stats2[2], 
                        'std:', stats2[3])
                  fig = plt.figure(figsize =(10, 7))
                  plt.boxplot(data)
                  plt.title(class_name)
                  plt.ylabel('IC struct')
                  plt.xticks([1, 2], ['LASSO', 'GIS'])
                  # show plot
                  plt.show()


def get_info_coefficients(dire1, dire2, X, n_splits, n_classes):

        union_genes, intersection_genes = get_union_inter(dire1, X, n_splits=n_splits)
        df_1= get_selected_genes(dire1, intersection_genes, union_genes, n_splits, n_classes) 

        union_genes, intersection_genes = get_union_inter(dire2, X, n_splits=n_splits)
        df_2= get_selected_genes(dire2, intersection_genes, union_genes, n_splits, n_classes) 

        not_matching_gis_all = []
        not_matching_lasso = []

        for i in range(df_1.shape[0]):

                print('---------------', str(i), '--------------------')
                genes_subset_1 = df_1.iloc[i][df_1.iloc[i]!=0].index.tolist()
                #print(df_1.iloc[i][df_1.iloc[i]!=0])
                print(len(genes_subset_1))
                genes_subset_2 = df_2.iloc[i][df_2.iloc[i]!=0].index.tolist()
                print(df_2.iloc[i][df_2.iloc[i]!=0].median())
                print(len(genes_subset_2))

                matching = [x for x in genes_subset_1 if x in genes_subset_2]
                not_matching_gis = [x for x in genes_subset_2 if x not in matching]
                print(df_2.iloc[i][not_matching_gis])
                print('Shared genes: ', len([x for x in genes_subset_1 if x in genes_subset_2]))

                # paired t-test
                # measuring the effect of an applied measure
                t_test = ttest_rel(df_1.iloc[i][matching].values, df_2.iloc[i][matching].values)
                print(t_test.pvalue)

                # Wilcoxon signed-rank test
                diff = np.around(df_1.iloc[i][matching].values- df_2.iloc[i][matching].values, decimals=10)
                t_ranks = wilcoxon(diff)
                print(t_ranks.pvalue)
                #print('T test coefficients: ', t_test.)

                df_plot =abs(df_2.iloc[i][df_2.iloc[i]!=0]).sort_values()
                df_plot_2=abs(df_1.iloc[i][df_1.iloc[i]!=0]).sort_values()
                ax = df_plot.plot(figsize=(15,7), style='o-', alpha=0.4, label='GIS')
                not_matching_gis_idx = [i for i, x in enumerate(df_plot.index) if x in not_matching_gis]
                matching_gis_idx = [i for i, x in enumerate(df_plot_2.index) if x in matching]
                matching_gis_idx_2 = [i for i, x in enumerate(df_plot.index) if x in matching]
                #print(not_matching_gis_idx)
                ax.scatter(not_matching_gis_idx, df_plot.iloc[not_matching_gis_idx], color='red')
                ax.scatter(matching_gis_idx_2, df_plot_2.iloc[matching_gis_idx], color='orange', label='LASSO')
                ax.set_title('Coefficients Distribution')
                ax.set_xlabel('Genes')
                ax.set_ylabel('Coefficients (abs)')
                ax.legend()
                plt.show()
                not_matching_gis_all.extend(not_matching_gis)
                not_matching_lasso.extend([x for x in genes_subset_1 if x not in matching])
                
        return list(set(not_matching_gis_all)), list(set(not_matching_lasso))

def view_enrichment_results(results_dire):
    
    results_dire_file = os.path.join(results_dire, 'enrichment_results.csv')
    df = pd.read_csv(results_dire_file)
    print(df)
    return df