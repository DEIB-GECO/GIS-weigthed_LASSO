import numpy as np
import pandas as pd
import tqdm
import json
from scipy.stats import fisher_exact
from itertools import chain
import time

def load_txt(filename):
    content = []
    with open(filename)as f:
        for line in f:
            content.append(line.strip())
    return content

def load_dict(filename):

    with open(filename, 'r') as fp:
        dictionary = json.load(fp)
        #print('Loaded dictionary')
    return  dictionary 

def load_annotations_file(filename):    
    with open(filename) as f:
        ann = {}
        for line in f:
            entrez, go = line.strip().split(',')
            if entrez in ann:
                ann[entrez].append(go)
            else:
                ann[entrez] = []
                ann[entrez].append(go)
    #print('N of genes: ', len(ann.keys()))
    return ann

def get_indexes(go_genes,go_terms):
    gene_index = {}
    count = 0
    for g in go_genes:
        gene_index[g] = count
        count += 1
    
    go_index = {}
    count = 0
    for g in go_terms:
        go_index[g] = count
        count += 1
    return gene_index, go_index

def load_dataset(filename, idx=0):

    df = pd.read_csv(filename)
    #print(df.head())
    df_ind = df.set_index(df.columns[idx])
    #print(df_ind.head())
    return df_ind
    
class go_matrix:

    def __init__(self, go_file, alt_file, ids_file, genes_list):
        self.go_file = go_file
        self.alt_file = alt_file
        self.ids_file = ids_file
        self.genes_list = genes_list

    def build_matrix(self, gene_go_file):

        d_go_anns = load_dict(self.go_file)
        d_go_alts = load_dict(self.alt_file)
        d_gene_anns = load_annotations_file(gene_go_file)
        go_genes = d_gene_anns.keys()
        go_terms = d_go_anns.keys()
        gene_index, go_index=  get_indexes(go_genes,go_terms)
        matrix = np.zeros((len(gene_index), len(go_index)))
       
        for g in tqdm.tqdm(d_gene_anns):
            for go in d_gene_anns[g]:

                if go in go_index:
                    indice_go = go_index[go]
                    new_go = go
                else:
                    if go in d_go_alts:
                        indice_go = go_index[d_go_alts[go]]
                        new_go = d_go_alts[go]
                matrix[gene_index[g], indice_go] = 1.
                for anc in d_go_anns[new_go]:
                    if anc in d_go_anns:
                        indice_go = go_index[anc]
                    else:
                        if anc in d_go_alts:
                            indice_go = go_index[d_go_alts[anc]]
                    matrix[gene_index[g], indice_go] = 1.
            
        cc = list(map(lambda x : x[0], sorted(list(go_index.items()), key = lambda x : x[1]) ))
        ii = list(map(lambda x : x[0], sorted(list(gene_index.items()), key = lambda x : x[1]) ))
        self.df = pd.DataFrame(matrix, columns = cc, index = ii)
        #self.df = self.df.drop_duplicates()
        #return self.df

    def build_weighted_matrix(self, gene_go_file, IC_scores_file):
        
        d_go_anns = load_dict(self.go_file)
        d_go_alts = load_dict(self.alt_file)
        d_gene_anns = load_annotations_file(gene_go_file)
        go_genes = d_gene_anns.keys()
        go_terms = d_go_anns.keys()
        gene_index, go_index=  get_indexes(go_genes,go_terms)
    
        matrix = np.zeros((len(gene_index), len(go_index)))
        # check if duplucates are here!
        
        scores = load_dataset(IC_scores_file)
        #print(scores.head())
        scores_sorted = {k:scores.loc[k] for k in d_go_anns}

        for g in tqdm.tqdm(d_gene_anns):
       
            for go in d_gene_anns[g]:

       
                if go in go_index:
                    indice_go =go_index[go]
                    new_go = go
                    score_idx = go

                else:
                    if go in d_go_alts:
                        indice_go = go_index[d_go_alts[go]]
                        new_go = d_go_alts[go]
                        score_idx = d_go_alts[go]
                #print(scores_sorted[score_idx])
             
                matrix[gene_index[g], indice_go] = scores_sorted[score_idx]
    
                for anc in d_go_anns[new_go]:
                    if anc in d_go_anns:
                        indice_go = go_index[anc]
                        score_idx= anc
                    else:
                        if anc in d_go_alts:
                            indice_go = go_index[d_go_alts[anc]]
                            score_idx= d_go_alts[anc]
                    matrix[gene_index[g], indice_go] = scores_sorted[score_idx]
            
        cc = list(map(lambda x : x[0], sorted(list(go_index.items()), key = lambda x : x[1]) ))
        ii = list(map(lambda x : x[0], sorted(list(gene_index.items()), key = lambda x : x[1]) ))
        self.df_scores = pd.DataFrame(matrix, columns = cc, index = ii)
        #self.df_scores = self.df_scores.drop_duplicates()
        #return self.df_scores

    

    def get_symbol_matrix(self, m, return_df= False):

        sys_df =  load_dataset('/home/mongardi/Metagene/GO_ann/genes_and_ids_all_red.csv', idx=1)
        matrix_new = pd.DataFrame(np.zeros((len(self.genes_list), m.shape[1])),index = self.genes_list, columns= m.columns)
        #print(sys_df.head())
        # build gis dict
        #matrix = m.drop_duplicates()
        matrix = m
        print(matrix.shape)
        sum = 0
        for gene in tqdm.tqdm(self.genes_list):

            if gene not in sys_df.index:
                #if gene.split('.1')[0] in merged.index:
                    
                    #d_gis[gene] = merged.loc[gene.split('.1')[0]]['0']
                
                #else:
                #print(gene)
                sum = sum +1
                matrix_new.loc[gene] = np.zeros(matrix.shape[1])

            else: 

                df_g = sys_df.loc[gene]
                ids =  df_g['id']

                            
                if df_g.shape[0] > 1:
                    
                    ids = [str(x) for x in ids if str(x) in matrix.index.tolist()]
                    if len(ids) > 0:
                        
                        ids_cs = []
                        for id in ids:
                            #print(id)
                            ids_cs.append(np.sum(matrix.loc[id]))
                        #print(ids_cs)
                        best_id = ids[np.argmax(ids_cs)]
                        matrix_new.loc[gene] = matrix.loc[str(best_id)]
                    else: 
                        matrix_new.loc[gene] = np.zeros(matrix.shape[1])
                        sum = sum +1

                elif str(ids) in matrix.index.tolist():
                    
                    
                    matrix_new.loc[gene] =  matrix.loc[str(ids)]

                else: 
                    matrix_new.loc[gene] = np.zeros(matrix.shape[1])
                    sum = sum +1
            

        if return_df:
            return matrix_new
        else: 
            self.matrix_new = matrix_new


    def get_weighted_matrix(self):

        all_ids =  self.get_gene_ids_from_sys(self.genes_list, verbose=False)
        df_all_ids = self.df_scores.loc[all_ids]
        return df_all_ids
    
    def get_GIS_scores(self):

        scrs = self.df_scores.sum(axis=1)/self.df.sum(axis=1)
        #scrs_genes = scrs.loc[gene_set]
        return scrs
    
    def get_GIS_scores_symbol(self):

        print('Calculating GIS scores .....')
        d1 = self.get_symbol_matrix(self.df_scores, return_df=True)
        d2 = self.get_symbol_matrix(self.df, return_df=True)
        scrs = d1.sum(axis=1)/d2.sum(axis=1)
        #scrs_genes = scrs.loc[gene_set]
        return scrs
    
    def get_gene_ids_from_sys(self, genes_list, verbose=True):

        sys_df =  load_dataset(self.ids_file, idx=0)
        gis_df = load_dataset('/home/mongardi/Metagene/GO_ann/gis_scores_norm.csv')
        merged = sys_df.merge(gis_df, right_on=gis_df.index, left_index=True)
        merged = merged.sort_values('0', ascending=False).drop_duplicates('symbol')
        merged = merged.set_index('symbol')
        #print(merged.shape[0])
        #print(merged.head())

        #genes_list = [x.split('.1')[0] for x in genes_list]
        #print(len(genes_list))
        #genes_sys= list(set(genes_list) & set(merged.index))
        genes_sys = [gene for gene in genes_list if gene in merged.index]
        #print([gene for gene in genes_list if gene not in merged.index])
        if verbose:
            print('Number of selected genes with GO annotations:' ,len(genes_sys))
        else:
            print('The total number of genes with GO annotations in the dataset:' ,len(genes_sys))
        
        genes_ids =  merged.loc[genes_sys][merged.columns[0]].values
        genes_ids = [str(x) for x in genes_ids]
        return genes_ids
    
    def get_dict_gene_sys_from_ids(self, genes_list, verbose=True):

        sys_df =  load_dataset(self.ids_file, idx=0)
        gis_df = load_dataset('/home/mongardi/Metagene/GO_ann/gis_scores_norm.csv')
        merged = sys_df.merge(gis_df, right_on=gis_df.index, left_index=True)
        merged = merged.sort_values('0', ascending=False).drop_duplicates('symbol')
        merged = merged.set_index('symbol')
        #print(merged.shape[0])
        #print(merged.head())

        #genes_list = [x.split('.1')[0] for x in genes_list]
        #print(len(genes_list))
        #genes_sys= list(set(genes_list) & set(merged.index))
        genes_sys = [gene for gene in genes_list if gene in merged.index]
        #print([gene for gene in genes_list if gene not in merged.index])
        if verbose:
            print('Number of selected genes with GO annotations:' ,len(genes_sys))
        else:
            print('The total number of genes with GO annotations in the dataset:' ,len(genes_sys))
        
        genes_ids =  merged.loc[genes_sys][merged.columns[0]].values
        genes_ids = [str(x) for x in genes_ids]
        return genes_ids
    
    def get_sum_annotations(self, genes_subset):
        
        ids = self.get_gene_ids_from_sys(genes_subset)
        tot_anns = sum(self.df.loc[ids].sum(axis=1))
        return tot_anns

   

    def get_specific_anns(self, genes_subset):

        ids = self.get_gene_ids_from_sys(genes_subset)
        d_go_anns = load_dict(self.go_file)
        all_anns = []
        for id in ids:
            all_anns.append(self.df.loc[id][self.df.loc[id] > 0].index)
        all_anns =  list(chain.from_iterable(all_anns))
        anns_u = np.unique(all_anns)
        anns_u = [x for x in anns_u if x.startswith('GO')]
        specific_anns = []
        for ann in anns_u:
            specific = 1
            for g in anns_u:
                if ann in d_go_anns[g]:
                    specific = 0
                    break;
            if specific ==1:
                specific_anns.append(ann)
        #print(specific_anns)
        print('Number of specific annotations: ', len(specific_anns))
        return specific_anns
    
    def get_terms_info(self, terms, IC_scores_file):

        scores = load_dataset(IC_scores_file)
        d_score_terms = {k:scores.loc[k].values[0] for k in terms}

        return d_score_terms
    
    
    def get_significant_terms(self, genes_subset, universe='all', alternative = 'two-sided',
                              correction='B', alpha = 0.05):


        if universe == 'all': 
            df_symbols = self.matrix_new
            print(df_symbols.shape)
            ids = genes_subset
        
        elif universe == 'annotated':
            df_symbols = self.matrix_new[self.matrix_new.sum(axis=1) != 0]
            print(df_symbols.shape)
            ids = [x for x in genes_subset if x in df_symbols.index]

        all_anns = []
        for id in ids:
            all_anns.append(df_symbols.loc[id][df_symbols.loc[id] > 0].index)
        all_anns =  list(chain.from_iterable(all_anns))
        print('Number of total annotations:', len(all_anns))
        anns_u = np.unique(all_anns)
        anns_u = [x for x in anns_u if x.startswith('GO')]
        test = {}
        
        s = len(ids)
        #print(len(all_ids))
        #print(s)
        g = df_symbols.shape[0] - s
        print('Total: ', len(self.genes_list), '=> Check sum: ', s + g)
        #print(g)
        print('Number of total unique annotations:', len(anns_u), '/', self.df.shape[1])
        self.df_ids = df_symbols.loc[ids]
        self.df_all_ids = df_symbols

        for ann in tqdm.tqdm(anns_u):
            #start_time = time.time()
            n_s_a = self.df_ids[ann].sum()
            n_s_na = s - n_s_a
            n_o_a = self.df_all_ids[ann].sum() - n_s_a
            n_o_na = g - n_o_a
            data = [[n_s_a, n_s_na ],[n_o_a, n_o_na]]
            #print(data)
            #print(data)
            #print("--- %s seconds ---" % (time.time() - start_time))
            #start_time = time.time()
            if data[1][0] >= 0:
                # _, p_value,_,_ = chi2_contingency(data, correction= False)
                res = fisher_exact(data, alternative=alternative)
                # res = fisher_exact(data, alternative='greater')
                p_value = res.pvalue
                #print("--- %s seconds ---" % (time.time() - start_time))
            else:
                #print(ann)
                #print(data)
                data[1][0] = 0
                res = fisher_exact(data, alternative=alternative)
                p_value = res.pvalue
                #p_value = None

            test[ann] = p_value
        n_test = len([k for k in test if test[k] is not None])
        print('Number of test performed: ', n_test)
        if correction == 'B':

            test_corrected = {k:test[k]*n_test for k in test if test[k] is not None}
            significant = {k:test_corrected[k] for k in test_corrected if test_corrected[k] < alpha}

        if correction == 'BH':
             
            test_true = {k:test[k]for k in test if test[k] is not None}
            #print(test_true)
            test_true_sorted = dict(sorted(test_true.items(), key=lambda x:x[1]))
            #print(test_true_sorted)
            test_corrected = {k:test_true_sorted[k]* n_test /(i+1) for i,k in enumerate(test_true_sorted)}        
            #print(test_corrected)                
            significant = {k:test_corrected[k] for k in test_corrected if test_corrected[k] < alpha}
            #print(significant)
        print('Number of significant annotations with ' + correction + ': ',len(significant))
        print('Percentage: ', 100*len(significant)/n_test)
        return significant, n_test