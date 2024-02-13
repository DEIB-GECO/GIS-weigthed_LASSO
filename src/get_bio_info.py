import numpy as np
import pandas as pd 
from Parser import * 
from utils.utils import *

from prior_knowledge.create_go_matrix import go_matrix
from prior_knowledge.create_hpo_matrix import hpo_matrix
from prior_knowledge.create_kegg_pathways_matrix import KEGG_pathways_matrix
from prior_knowledge.create_reactome_pathways_matrix import Reactome_pathways_matrix

#args = get_args()

def build_df_results(n_classes):
    df = pd.DataFrame(np.zeros((n_classes, 10)), columns= ['Subtype', 'N', 
                                                            'GO', 'GO sig', 'KEGG', 'KEGG sig', 'REACTOME', 'REACTOME sig', 'HPO', 'HPO sig'])
    return df    

def get_bio_stats(args, genes_list, training_results, results_dire, labels):
    
  
    # GO 
    get_go_anns_matrix = go_matrix(args.go_file_dir, args.go_alt_file_dir, 
                                   args.genes_ids_file_dir, genes_list)
    get_go_anns_matrix.build_matrix(args.go_annotations_file_dir)
    get_go_anns_matrix.get_symbol_matrix(get_go_anns_matrix.df)

    # KEGG

    get_kegg_pathways_matrix = KEGG_pathways_matrix(args.kegg_hsa_pathways_file_dir, 
                                                    args.kegg_gene_pathways_file_dir,
                                                    args.genes_ids_file_dir,
                                                    genes_list)
    get_kegg_pathways_matrix.build_matrix()
    get_kegg_pathways_matrix.get_symbol_matrix(get_kegg_pathways_matrix.pm)

    # Reactome

    get_reactome_paths_matrix = Reactome_pathways_matrix(args.reactome_hsa_pathways_file_dir, 
                                                         args.genes_ids_file_dir, 
                                                         genes_list)
    get_reactome_paths_matrix.build_matrix()
    get_reactome_paths_matrix.get_symbol_matrix(get_reactome_paths_matrix.pm)
    
    # hpo 
    get_hpo_anns_matrix = hpo_matrix(args.hpo_file_dir, args.hpo_alt_file_dir, 
                                     args.genes_ids_file_dir, genes_list)
    get_hpo_anns_matrix.build_matrix(args.hpo_annotations_file_dir)
    get_hpo_anns_matrix.get_symbol_matrix(get_hpo_anns_matrix.df)
    df_results = build_df_results(training_results.shape[0]+1)

    for i in range(training_results.shape[0]):
        genes_subset = training_results.iloc[i][training_results.iloc[i]!= 0].index
        store_significant = {}
        results = []
        results.append(labels[i])
        #print(genes_subset)
        #print(genes_list)
        print('--------------------------')
        print(labels[i])
        print('--------------------------')
        print('Number of selected genes:', len(genes_subset))
        results.append(len(genes_subset))

        # GO
        print('------------GO------------')
        print('----significant annotations----')
        significant_terms, n_test = get_go_anns_matrix.get_significant_terms(genes_subset, correction='BH')
        results.append(n_test)
        results.append(len(significant_terms))
        store_significant['GO'] = significant_terms
        if len(significant_terms) > 0:
            info_terms = get_go_anns_matrix.get_terms_info(significant_terms, args.IC_go_dir)
            print(get_stats(list(info_terms.values())))
        else: 
            info_terms = []
        save_dictionary(info_terms, results_dire + '/significant_terms_info_'+ labels[i] +'.txt')

        print()

        # KEGG
        print('------------KEGG------------')
        significant_paths, n_test = get_kegg_pathways_matrix.get_significant_pathways(genes_subset, correction='BH')
        results.append(n_test)
        results.append(len(significant_paths))
        store_significant['KEGG'] = significant_paths
        print()

        # Reactome 
        print('------------Reactome------------')
        significant_paths, n_test =  get_reactome_paths_matrix.get_significant_pathways(genes_subset, correction='BH')
        results.append(n_test)
        results.append(len(significant_paths))
        store_significant['REACTOME'] = significant_paths

        # Reactome 
        print('------------hpo------------')
        significant_terms, n_test =  get_hpo_anns_matrix.get_significant_terms(genes_subset, correction='BH')
        results.append(n_test)
        results.append(len(significant_terms))
        store_significant['HPO'] = significant_terms
        save_dictionary(store_significant, results_dire + '/significant_terms_'+ labels[i] +'.txt')

        df_results.iloc[i] = results

    # overall 
    store_significant = {}
    results = []
    results.append('all')
    #print(genes_subset)
    #print(genes_list)
    print('--------------------------')
    print('all')
    print('--------------------------')
    print('Number of selected genes:', len(training_results.columns.tolist()))
    results.append(len(training_results.columns.tolist()))

    # GO
    print('------------GO------------')
    print('----significant annotations----')
    significant_terms, n_test = get_go_anns_matrix.get_significant_terms(training_results.columns.tolist(), correction='BH')
    results.append(n_test)
    results.append(len(significant_terms))
    store_significant['GO'] = significant_terms
    if len(significant_terms) > 0:
        info_terms = get_go_anns_matrix.get_terms_info(significant_terms, '/home/mongardi/Metagene/GO_ann/IC_scores_normalized.csv')
        print(get_stats(list(info_terms.values())))
    else: 
        info_terms = []
    save_dictionary(info_terms, results_dire + '/significant_terms_info_all.txt')

    print()

    # KEGG
    print('------------KEGG------------')
    significant_paths, n_test = get_kegg_pathways_matrix.get_significant_pathways(training_results.columns.tolist(), correction='BH')
    results.append(n_test)
    results.append(len(significant_paths))
    store_significant['KEGG'] = significant_paths
    print()

    # Reactome 
    print('------------Reactome------------')
    significant_paths, n_test =  get_reactome_paths_matrix.get_significant_pathways(training_results.columns.tolist(), correction='BH')
    results.append(n_test)
    results.append(len(significant_paths))
    store_significant['REACTOME'] = significant_paths

    # Reactome 
    print('------------hpo------------')
    significant_terms, n_test =  get_hpo_anns_matrix.get_significant_terms(training_results.columns.tolist(), correction='BH')
    results.append(n_test)
    results.append(len(significant_terms))
    store_significant['HPO'] = significant_terms

    save_dictionary(store_significant, results_dire + '/significant_terms_all.txt')

    print(len(results))
    df_results.iloc[len(training_results.index)] = results
    df_results.to_csv(results_dire + '/enrichment_results.csv', index_label=False)    


