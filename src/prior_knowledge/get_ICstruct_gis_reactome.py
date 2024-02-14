import tqdm
import numpy as np
import pandas as pd
import sys
import os
import copy
import matplotlib.pyplot as plt


from create_reactome_pathways_matrix import Reactome_pathways_matrix


from utils.utils import *
sys.path.append("src")
from Parser import get_args
# set directory

args = get_args()
dire =  args.save_files_dir
get_reactome_paths_matrix = Reactome_pathways_matrix(args.reactome_hsa_pathways_file_dir, 
                                                         args.genes_ids_file_dir, 
                                                         [])
get_reactome_paths_matrix.build_matrix()
get_reactome_paths_matrix.get_symbol_matrix(get_reactome_paths_matrix.pm)
human_pathways = get_reactome_paths_matrix.pm.columns.tolist()


reactome_dag = build_reactome_dag(human_pathways, args.reactome_relations_file_dir)
d_unf = copy.deepcopy(reactome_dag)
d_counts = compute_dag_counts(reactome_dag)
coverage = compute_coverage(reactome_dag, d_counts)
specifity = compute_specificity(d_unf)

IC_score_df = compute_IC_struct(reactome_dag, specifity, coverage)
IC_score_df.to_csv(os.path.join(dire, 'IC_scores_normalized_reactome.csv'))

# Reactome
