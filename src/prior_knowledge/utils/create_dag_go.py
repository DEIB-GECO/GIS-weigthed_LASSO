
import tqdm
import re
import numpy as np
import pandas as pd
import os
import json



def get_ancestors(pattern = r"GO:[0-9]+", file="go.obo"):
    d_dag= {}
    pattern = pattern
    with open(file) as f:
        name = None
        ancestors = {}
        for line in f:
            if line.startswith("id:"):
                if name is not None and  name.startswith('GO'):
                    d_dag[name] = ancestors
                ancestors = set()
                name = line.strip()[4:]

            if line.startswith("is_a: "):
                ancestors.add(line[6:16])
                
            if line.startswith("relationship: ") and ('part_of' in line):
                print('-------------')
                print(name)
                ancestors.add(re.search(pattern, line).group())
                print(re.search(pattern, line).group())
                print('-------------')
    
    # check len
    go_ids = [x for x in d_dag if x.startswith('GO')]
    print(len(d_dag))
    print(len(go_ids))
    if len(go_ids) != len(d_dag):
        raise TypeError("Check function!")
    
    return d_dag

def unfold_dag(pattern = r"GO:[0-9]+", file="go.obo"):

    d_dag = get_ancestors(pattern, file)
    for i in range(20):
        total_init = sum([len(d_dag[g]) for g in d_dag])

        for id in tqdm.tqdm(d_dag):
            ancs = list(d_dag[id])
            for a in ancs:
                d_dag[id] = d_dag[id].union(d_dag.get(a, set()))
        total_end = sum([len(d_dag[g]) for g in d_dag])

        if total_init==total_end:
            print(('Unfolding completed!'))
            break
    
    d_dag = {k:list(d_dag[k]) for k in d_dag}
    return d_dag

def get_alternative_ids(pattern = r"GO:[0-9]+", file="go.obo"):
    alt = {}
    pattern = pattern
    with open(file) as f:
        name = None
        alts = []
        for line in f:
            if line.startswith("id:"):
                if name is not None:
                    for a in alts:
                        alt[a] = name
                    alts = []
                name = line.strip()[4:]
            if line.startswith("alt_id: "): 
                alts.append(line[8:18])
    return alt 

#def get_specific_disease_phenotypes(disease_name):


def save_dict(dictionary,filename):
    with open(filename + ".txt", "w") as fp:
        json.dump(dictionary, fp)
        print("Done writing dict into .txt file")

if __name__ == "__main__":

    import sys
    sys.path.append("/home/mongardi/Metagene_repo/src")
    from Parser import get_args
    args = get_args()
    results_dire =  args.save_files_dir

    get_dag = True
    
    if get_dag:
        go_dag = unfold_dag(file = args.go_obo_file_dir)
        # save unfolded ontology
        save_dict(go_dag, os.path.join(results_dire, 'go_dag'))

        alternative_ids = get_alternative_ids(file = args.go_obo_file_dir)
        save_dict(alternative_ids, os.path.join(results_dire, 'go_alternative_ids'))
    
    