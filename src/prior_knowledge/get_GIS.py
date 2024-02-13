import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

from create_go_matrix import go_matrix
from create_reactome_pathways_matrix import Reactome_pathways_matrix
from create_hpo_matrix import hpo_matrix
from utils.utils import * 

sys.path.append("/home/mongardi/Metagene_repo/src")
from Parser import get_args


args = get_args()
# set directory
dire = args.save_files_dir

# GO 
get_go_anns_matrix = go_matrix(args.go_file_dir, args.go_alt_file_dir, 
                                args.genes_ids_file_dir, [])
get_go_anns_matrix.build_matrix(args.go_annotations_file_dir)
get_go_anns_matrix.build_weighted_matrix(args.go_annotations_file_dir, os.path.join(dire, 'IC_scores_normalized_go.csv'))

scrs = get_go_anns_matrix.df_scores.sum(axis=1)/get_go_anns_matrix.df.sum(axis=1)
scrs.to_csv(os.path.join(dire, 'gene_scores_norm_go.csv'))

# reactome

get_reactome_paths_matrix = Reactome_pathways_matrix(args.reactome_hsa_pathways_file_dir, 
                                                         args.genes_ids_file_dir, 
                                                         [])
get_reactome_paths_matrix.build_matrix()
get_reactome_paths_matrix.get_symbol_matrix(get_reactome_paths_matrix.pm)
human_pathways = get_reactome_paths_matrix.pm.columns.tolist()
IC_scores_reactome = load_dataset(os.path.join(dire,'IC_scores_normalized_reactome.csv'))
IC_score_reactome_sorted = {k:IC_scores_reactome.loc[k].values[0] for k in get_reactome_paths_matrix.pm.columns}
df_scores_reactome = get_reactome_paths_matrix.pm.mul(IC_score_reactome_sorted)
# change indexes reactome
new_indexes = np.array([str(x) for x in  get_reactome_paths_matrix.pm.index])
df_scores_reactome_new = df_scores_reactome.set_index(new_indexes)
df_reactome_new = get_reactome_paths_matrix.pm.set_index(new_indexes)

srcs = df_scores_reactome_new.sum(axis=1)/ df_reactome_new.sum(axis=1)
scrs.to_csv(os.path.join(dire, 'gene_scores_norm_reactome.csv'))

# hpo 
get_hpo_anns_matrix = hpo_matrix(args.hpo_file_dir, args.hpo_alt_file_dir, 
                                    args.genes_ids_file_dir, [])
get_hpo_anns_matrix.build_matrix(args.hpo_annotations_file_dir)
get_hpo_anns_matrix.build_weighted_matrix(args.hpo_annotations_file_dir, os.path.join(dire, 'IC_scores_normalized_hpo.csv'))

scrs = get_hpo_anns_matrix.df_scores.sum(axis=1)/get_hpo_anns_matrix.df.sum(axis=1)
scrs.to_csv(os.path.join(dire, 'gene_scores_norm_hpo.csv'))

# Combined
# go + reactome 
merged_matrix = pd.concat([get_go_anns_matrix.df_scores, df_scores_reactome_new], axis=1)
print(merged_matrix.shape)
merged_matrix_binary =  pd.concat([get_go_anns_matrix.df, df_reactome_new], axis=1)
print(merged_matrix_binary.shape)
merged_matrix_new = merged_matrix.fillna(0)
merged_matrix_binary_new = merged_matrix_binary.fillna(0)
scrs = merged_matrix_new.sum(axis=1)/merged_matrix_binary_new.sum(axis=1)
scrs.to_csv(os.path.join(dire, 'gene_scores_norm_go_reactome.csv'))

# go + hpo
merged_matrix = pd.concat([get_go_anns_matrix.df_scores, get_hpo_anns_matrix.df_scores], axis=1)
print(merged_matrix.shape)
merged_matrix_binary =  pd.concat([get_go_anns_matrix.df,  get_hpo_anns_matrix.df], axis=1)
print(merged_matrix_binary.shape)
merged_matrix_new = merged_matrix.fillna(0)
merged_matrix_binary_new = merged_matrix_binary.fillna(0)
scrs = merged_matrix_new.sum(axis=1)/merged_matrix_binary_new.sum(axis=1)
scrs.to_csv(os.path.join(dire, 'gene_scores_norm_go_hpo.csv'))

# go + reactome + hpo
merged_matrix = pd.concat([get_go_anns_matrix.df_scores, df_scores_reactome_new, get_hpo_anns_matrix.df_scores], axis=1)
print(merged_matrix.shape)
merged_matrix_binary =  pd.concat([get_go_anns_matrix.df, df_reactome_new,  get_hpo_anns_matrix.df], axis=1)
print(merged_matrix_binary.shape)
merged_matrix_new = merged_matrix.fillna(0)
merged_matrix_binary_new = merged_matrix_binary.fillna(0)
scrs = merged_matrix_new.sum(axis=1)/merged_matrix_binary_new.sum(axis=1)
scrs.to_csv(os.path.join(dire, 'gene_scores_norm_go_reactome_hpo.csv'))