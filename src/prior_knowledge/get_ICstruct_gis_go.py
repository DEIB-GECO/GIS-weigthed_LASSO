import tqdm
import numpy as np
import pandas as pd
import sys
import os
import copy
import matplotlib.pyplot as plt

from create_go_matrix import go_matrix
from utils.create_dag_go import unfold_dag
from utils.utils import *
sys.path.append("src")
from Parser import get_args
# set directory

args = get_args()
dire =  args.save_files_dir

# build dag
file = os.path.join(dire, 'go.obo')
go_dag = unfold_dag(file=file)
d_unf = copy.deepcopy(go_dag)

d_class = get_sub_ontology_go(file)
coverage = compute_coverage_go(go_dag, d_class)
specifity = compute_specificity(d_unf, d_class=d_class, go=True)

IC_score_df = compute_IC_struct(go_dag, specifity, coverage)
IC_score_df.to_csv(os.path.join(dire, 'IC_scores_normalized_go.csv'))

