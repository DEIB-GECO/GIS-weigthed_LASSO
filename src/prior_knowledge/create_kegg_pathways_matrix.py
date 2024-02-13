import numpy as np
import pandas as pd
import tqdm
import json
from scipy.stats import fisher_exact
from itertools import chain



class KEGG_pathways_matrix:
    def __init__(self, pathways_file, filename, ids_file, genes_list):
        self.pathways_file = pathways_file
        self.filename = filename
        self.ids_file = ids_file
        self.genes_list = genes_list

    def load_txt(self, file):
        content = []
        with open(file)as f:
            for line in f:
                content.append(line.strip())
        return content


    def load_dataset(self, idx=0):

        df = pd.read_csv(self.filename, low_memory=False)
        df_ind = df.set_index(df.columns[0].split())
        return df_ind
    
    def load_sys_dataset(self, idx=0):

        df = pd.read_csv(self.ids_file)
        #print(df.head())
        df_ind = df.set_index(df.columns[idx])
        #print(df_ind.head())
        return df_ind
    
    def build_empty_matrix(self):
        
        pathways= self.load_txt(self.pathways_file)
        self.df = self.load_dataset()
        genes = np.unique(self.df.index) 
        pm = pd.DataFrame(np.zeros([len(genes), len(pathways)]),index = genes, columns= pathways)
        return pm 
    
    def build_matrix(self):
        self.pm = self.build_empty_matrix()
        for gene in  tqdm.tqdm(self.pm.index):
            paths = self.df.loc[gene][self.df.loc[gene].notnull()]
            if len(paths) > 0:
                self.pm.loc[gene][paths] = 1.0
        #return self.pm 

    def get_sum(self, gene_list):

        sum_genes = self.pm.sum(axis=1)
        total_sum = sum_genes.sum()
        #print('Total number of pathway counts: ', total_sum)
        #print(sum_genes)
        return sum_genes
    
    def get_gene_ids_from_sys(self, genes_list, verbose=True):

        sys_df =  self.load_sys_dataset(idx=0)
        pathways_df = self.get_sum(genes_list).to_frame()
        #print(pathways_df)
        merged = sys_df.merge(pathways_df, right_on=pathways_df.index, left_index=True)
        merged = merged.sort_values(0, ascending=False).drop_duplicates('symbol')
        merged = merged.set_index('symbol')

        #genes_list = [x.split('.1')[0] for x in genes_list]
        #genes_sys= list(set(genes_list) & set(merged.index))
        genes_sys = [gene for gene in genes_list if gene in merged.index]
        if verbose:
            print('Number of selected genes with KEGG annotations:' ,len(genes_sys))
        else:
            print('The total number of genes with KEGG annotations in the dataset:' ,len(genes_sys))
        genes_ids =  merged.loc[genes_sys][merged.columns[0]].values
        genes_ids = [x for x in genes_ids]
        return genes_ids
    
    '''def get_gene_ids_from_sys(self, genes_list):

        sys_df =  self.load_sys_dataset(idx=-1)
        intersection = list(set(genes_list) & set(sys_df.index))
        # check 
        if len(intersection) != len(np.unique(intersection)):
            # to do: select idx with more annotations
            print('To implement')
            exit

        genes_ids =  sys_df.loc[intersection]['id'].values
        #genes_ids = [str(x) for x in genes_ids]
        #print(genes_ids)
        genes_ids = list(set(genes_ids) & set(self.pm.index))
        print('Number of genes with reactome pathways:' ,len(genes_ids))
        return genes_ids '''
    
    def get_matching_ids(self, genes_subset, matching_set):
        ids = self.get_gene_ids_from_sys(genes_subset, verbose=False)
        matching_set = [int(x) for x in matching_set]
        matching_set_ids = [x for x in ids if x in matching_set]
        print('Number of matches: ', len(matching_set_ids))
        #print(len(genes_subset))
        print('Percentages of matching genes: ', len(matching_set_ids)/len(genes_subset))
        return matching_set_ids
    
    
    def get_matching_sys(self, genes_subset, matching_set):
       
        matching_set_ids = [x for x in genes_subset if x in matching_set]
        print('Number of matches: ', len(matching_set_ids))
        #print(len(genes_subset))
        print('Percentages of matching genes: ', len(matching_set_ids)/len(genes_subset))
        return matching_set_ids
    

    def get_symbol_matrix(self, m):

        #sys_df =  self.load_dataset('/home/mongardi/Metagene/GO_ann/genes_and_ids_all_red.csv', idx=1)
        sys_df =  self.load_sys_dataset(idx=1)
        matrix_new = pd.DataFrame(np.zeros((len(self.genes_list), m.shape[1])),index = self.genes_list, columns= m.columns)
        #matrix = m.drop_duplicates()
        matrix = m
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
                
                    ids = [x for x in ids if x in matrix.index.tolist()]
                    if len(ids) > 0:
                    
                        ids_cs = []
                        for id in ids:

                            ids_cs.append(np.sum(matrix.loc[id]))
                        #print(ids_cs)
                        best_id = ids[np.argmax(ids_cs)]
                        matrix_new.loc[gene] = matrix.loc[best_id]
                    else: 
                        matrix_new.loc[gene] = np.zeros(matrix.shape[1])
                        sum = sum +1

                elif ids in matrix.index.tolist():
                    
                    matrix_new.loc[gene] =  matrix.loc[ids]

                else: 
                    matrix_new.loc[gene] = np.zeros(matrix.shape[1])
                    sum = sum +1
            
        self.matrix_new = matrix_new
    
    def get_significant_pathways(self, genes_subset, universe='all', alternative = 'two-sided', 
                                 correction='B', alpha = 0.05):

        
        if universe == 'all': 
            df_symbols = self.matrix_new
            print(df_symbols.shape)
            ids = genes_subset
        
        elif universe == 'annotated':
            df_symbols = self.matrix_new[self.matrix_new.sum(axis=1) != 0]
            print(df_symbols.shape)
            ids = [x for x in genes_subset if x in df_symbols.index]

       
        all_ps = []
        for id in ids:
            all_ps.append(df_symbols.loc[id][df_symbols.loc[id] > 0].index)
        all_ps =  list(chain.from_iterable(all_ps))
        #print('Number of total pathway counts: ', self.pm.loc[ids].sum(axis=1).sum())
        print('Number of total pathways:', len(all_ps))
        ps_u = np.unique(all_ps)
        test = {}
        s = len(ids)
        #print(len(all_ids))
        #print(s)
        g = df_symbols.shape[0] - s
        #print(g)
        print('Total: ', len(self.genes_list), '=> Check  sum: ', s + g)
        print('Number of total unique pathways:', len(ps_u), '/', self.pm.shape[1])
        pm_ids = df_symbols.loc[ids]
        pm_all_ids = df_symbols
        for p in tqdm.tqdm(ps_u):
            #n_s_a = self.pm.loc[ids][p].sum()
            n_s_a = pm_ids[p].sum()
            n_s_na = s - n_s_a
            #n_o_a = self.pm.loc[all_ids][p].sum() - n_s_a
            n_o_a = pm_all_ids[p].sum() - n_s_a
            n_o_na = g - n_o_a
            data = [[n_s_a, n_s_na ],[n_o_a, n_o_na]]
            if data[1][0] != -1:
                # _, p_value,_,_ = chi2_contingency(data, correction= False)
                res = fisher_exact(data, alternative=alternative)
                p_value = res.pvalue
            else:
                print(p)
                print(data)
                data[1][0] = 0
                res = fisher_exact(data, alternative=alternative)
                p_value = res.pvalue
                #p_value = None
        
            test[p] = p_value
        #print(test)
        n_test = len([k for k in test if test[k] is not None])
        print('Number of test performed: ', n_test)
        if correction == 'B':

            test_corrected = {k:test[k]*n_test for k in test if test[k] is not None}
            significant = {k:test_corrected[k] for k in test_corrected if test_corrected[k] < alpha}

        if correction == 'BH':
             
            test_true = {k:test[k]for k in test if test[k] is not None}
            test_true_sorted = dict(sorted(test_true.items(), key=lambda x:x[1]))
            test_corrected = {k:test_true_sorted[k]* n_test /(i+1) for i,k in enumerate(test_true_sorted)}                        
            significant = {k:test_corrected[k] for k in test_corrected if test_corrected[k] < alpha}

        print('Number of significant pathways with ' + correction + ': ',len(significant))
        print('Percentage: ', 100*len(significant)/n_test)
        return significant, n_test