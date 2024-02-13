import numpy as np
import pandas as pandas
import sys
import os
import copy
import matplotlib.pyplot as plt

from create_hpo_matrix import hpo_matrix
from utils.create_dag_hpo import unfold_dag
from utils.utils import *


sys.path.append("/home/mongardi/Metagene_repo/src")
from Parser import get_args

args = get_args()
# set directory
dire =  args.save_files_dir

# build dag
file = os.path.join(dire, 'hp.obo')
hpo_dag = unfold_dag(file=file)
d_unf = copy.deepcopy(hpo_dag)

d_counts = compute_dag_counts(hpo_dag)
coverage = compute_coverage(hpo_dag, d_counts)
specifity = compute_specificity(d_unf)

IC_score_df = compute_IC_struct(hpo_dag, specifity, coverage)
IC_score_df.to_csv(os.path.join(dire, 'IC_scores_normalized_hpo.csv'))

# compute gene annotation scores

get_hpo_anns_matrix = hpo_matrix(args.hpo_file_dir, args.hpo_alt_file_dir, 
                                    args.genes_ids_file_dir, [])
get_hpo_anns_matrix.build_matrix(args.hpo_annotations_file_dir)
get_hpo_anns_matrix.build_weighted_matrix(args.hpo_annotations_file_dir, os.path.join(dire, 'IC_scores_normalized_hpo.csv'))

scrs = get_hpo_anns_matrix.df_scores.sum(axis=1)/get_hpo_anns_matrix.df.sum(axis=1)
scrs.to_csv(os.path.join(dire, 'gene_scores_norm_hpo.csv'))