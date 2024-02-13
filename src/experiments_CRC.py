import numpy as np
import pandas as pd
import os
from sklearn.utils import shuffle
from run_experiments import run_experiments
from Parser import get_args

args =  get_args()

# load dataset
df = pd.read_csv(args.CRC_dataset_dir, index_col=args.CRC_dataset_index_col)

df = shuffle(df, random_state=42)

X = df.drop(['class'], axis=1)
y = df['class']


rpm = X.div(X.sum(axis=1).values, axis=0) *1e6
rpm_log = np.log2(rpm + 1)

#results_dire = os.path.join(args.CRC_results_dir, 'multi_lasso_test')
#selected_features = run_experiments(args, rpm_log, y, results_dire, tol =1e-4, lambdas=0.5, n_splits = 1, is_GIS=False, bio_info=True)

results_dire = os.path.join(args.CRC_results_dir, 'multi_lasso_gis_go_reactome_hpo_test_paper')
selected_features = run_experiments(args, rpm_log, y, results_dire, tol =1e-4, lambdas=0.5, n_splits = 10, is_GIS=True, bio_info=True, go=True, reactome=True, hpo=True)
