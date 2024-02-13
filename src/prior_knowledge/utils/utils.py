import numpy as np
import tqdm
import math
import pandas as pd


def load_dataset(filename, idx=0):

    df = pd.read_csv(filename)
    #print(df.head())
    df_ind = df.set_index(df.columns[idx])
    #print(df_ind.head())
    return df_ind


def get_sub_ontology_go(file):
    d_class={}
    name = None
    ancestors = {}
    with open(file) as f:
        for line in f:
            if line.startswith("id:"):
                if name is not None:
                    d_class[name] = ancestors
                
                name = line.strip()[4:]
            if line.startswith("namespace: "):
                
                if line[11:].startswith("biological_process"):
                    ancestors = 'BP'
                if line[11:].startswith("molecular_function"):
                    ancestors = 'MF'
                if line[11:].startswith("cellular_component"):
                    ancestors = 'CC'              
    return d_class

def get_anc(d_unf, node, d_class= None, go=False):
    
    # base condition
    if len(d_unf[node]) == 0:
        return node, 1
    node_arr = []
    depth_arr = []
    # do something

    if go:
        for anc in d_unf[node]:
            
            if d_class[anc] == d_class[node]:
                n, d = get_anc(d_unf, anc, d_class, go) 
            
                node_arr.append(n)
                depth_arr.append(d)
            
        max_depth_index = np.argmax(depth_arr)
   
    else:
        for anc in d_unf[node]:
            
        
            n, d = get_anc(d_unf, anc) 
            
            node_arr.append(n)
            depth_arr.append(d)
            
        max_depth_index = np.argmax(depth_arr)
    # return
    return node_arr[max_depth_index], depth_arr[max_depth_index] + 1


def compute_dag_counts(dag):
    d_counts={}
    for ann in tqdm.tqdm(dag):
        
        count=0
        for g in dag:
            if ann in dag[g]:
                count +=1
        d_counts[ann] = count

    return d_counts

def compute_coverage(dag, d_counts):
    coverage = {}

    for ann in tqdm.tqdm(dag):
        coverage[ann]= (1 - (math.log(d_counts[ann]+1, 10)/math.log(len(dag), 10)))

    return coverage

def compute_coverage_go(dag, d_class):
    
    ## check summing to total
    values = [d_class[x] for x in dag]
    n_bp = values.count('BP')
    n_mf = values.count('MF')
    n_cc = values.count('CC')  
    if n_bp + n_mf + n_cc != len(dag):
        raise TypeError("Check function!")
    
    d_counts={}

    for ann in tqdm.tqdm(dag):
        
        count=0
        for g in dag:
            if ann in dag[g] and d_class[ann]== d_class[g]:
                count +=1
        d_counts[ann] = count
    
    coverage = {}

    for ann in tqdm.tqdm(dag):

        ont_class = d_class[ann]
        if ont_class == 'BP':
            counts = n_bp
        if ont_class == 'MF':
            counts = n_mf
        if ont_class == 'CC':
            counts = n_cc
    
        coverage[ann]= (1 - (math.log(d_counts[ann]+1, 10)/math.log(counts, 10)))

    return coverage



def compute_specificity(dag, d_class= None, go=False, normalized=True):
    specificity={}

    for node in tqdm.tqdm(dag):
    
        start, max_depth = get_anc(dag, node, d_class=d_class, go=go)
        specificity[node] = max_depth


    if normalized:

        if go:
            specificity_normalized = {}
            max_mf = np.max([specificity[x] for x in specificity if d_class[x] == 'MF'])
            max_bp = np.max([specificity[x] for x in specificity if d_class[x] == 'BP'])
            max_cc = np.max([specificity[x] for x in specificity if d_class[x] == 'CC'])

            for term in specificity:

                if d_class[term] == 'MF':
                    specificity_normalized[term] = specificity[term]/max_mf
                
                
                if d_class[term] == 'BP':
                    specificity_normalized[term] = specificity[term]/max_bp
                    
                
                if d_class[term] == 'CC':
                    specificity_normalized[term] = specificity[term]/max_cc

        else:

            max_depth_dag = max(list(specificity.values()))
            print('Max depth ontology:', max_depth_dag)

            specificity_normalized = {k: specificity[k]/max_depth_dag for k in specificity}

        return specificity_normalized

    else: 

        return specificity

def compute_IC_struct(dag, specificity, coverage):

    IC_score = {k:specificity[k]*coverage[k] for k in dag}
    IC_score_df = pd.DataFrame.from_dict(IC_score, orient='index')

    return IC_score_df


def build_reactome_dag(human_pathways, file):
    paths_dag = {path:[] for path in human_pathways}
    with open(file) as f:
        for line in f:
            f = line.strip().split()[0]
            c = line.strip().split()[1]
            if f in human_pathways:

                if f in paths_dag:
                    paths_dag[f].append(c)
                else:
                    paths_dag[f] = []
                    paths_dag[f].append(c)
                
                if c not in paths_dag:
                    paths_dag[c] = []
    return paths_dag