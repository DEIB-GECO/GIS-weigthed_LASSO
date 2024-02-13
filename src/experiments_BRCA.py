import numpy as np
import pandas as pd
import os
from sklearn.utils import shuffle
from run_experiments import run_experiments
from Parser import get_args

args =  get_args()

# load dataset

df = pd.read_csv(args.BRCA_dataset_dir, index_col=args.BRCA_dataset_index_col)
df.head()
df = shuffle(df, random_state=42)
df.head()

X = df.drop(['expert_PAM50_subtype'], axis=1)
y = y = df['expert_PAM50_subtype']
X.head()
y.head()

results_dire = os.path.join(args.BRCA_results_dir, 'multi_lasso')
selected_features = run_experiments(args, X, y, results_dire, tol =1e-4, lambdas=0.5, n_splits = 10, is_GIS=False, bio_info=True)

results_dire = os.path.join(args.BRCA_results_dir, 'multi_lasso_gis_go_reactome_hpo')
selected_features = run_experiments(args, X, y, results_dire, tol =1e-4, lambdas=0.5, n_splits = 10, is_GIS=True, bio_info=True, go=True, reactome=True, hpo=True)


