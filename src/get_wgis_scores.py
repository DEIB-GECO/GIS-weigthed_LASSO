
import numpy as np
import pandas as pd
import tqdm
from utils.utils import * 



def get_GIS_scores(args, genes_list, notebook=False, dire_1=None, dire_2=None):

    if notebook:
        sys_df = load_dataset(dire_1, idx=0)
        gis_df = load_dataset(dire_2)
    else:
        sys_df = load_dataset(args.genes_ids_file_dir, idx=0)
        gis_df = load_dataset(args.go_gis_dir)
    print(type(gis_df.index[0]))
    print(type(sys_df.index[0]))
    merged = sys_df.merge(gis_df, right_on=gis_df.index, left_index=True)
    merged = merged.sort_values('0', ascending=False).drop_duplicates('symbol')
    merged = merged.set_index('symbol')
    # build gis dict
    d_gis = {}
    sum = 0
    for gene in tqdm.tqdm(genes_list):
        if gene not in merged.index:
            #if gene.split('.1')[0] in merged.index:
                
                #d_gis[gene] = merged.loc[gene.split('.1')[0]]['0']
            
            #else:
            #print(gene)
            sum = sum +1
            d_gis[gene] = 0.0

        else: 
            d_gis[gene] = merged.loc[gene]['0']

    print('Number of genes with no GO annotations: ', sum)
    return d_gis


def get_GIS_scores_reactome(args, genes_list, notebook=False, dire_1=None, dire_2=None):

    if notebook:
        sys_df = load_dataset(dire_1, idx=0)
        gis_df = load_dataset(dire_2)
    else:
        sys_df = load_dataset(args.genes_ids_file_dir, idx=0)
        gis_df = load_dataset(args.reactome_gis_dir)
    print(type(gis_df.index[0]))
    print(type(sys_df.index[0]))
    merged = sys_df.merge(gis_df, right_on=gis_df.index, left_index=True)
    merged = merged.sort_values('0', ascending=False).drop_duplicates('symbol')
    merged = merged.set_index('symbol')
    # build gis dict
    d_gis = {}
    sum = 0
    for gene in tqdm.tqdm(genes_list):
        if gene not in merged.index:
            #if gene.split('.1')[0] in merged.index:
                
                #d_gis[gene] = merged.loc[gene.split('.1')[0]]['0']
            
            #else:
            #print(gene)
            sum = sum +1
            d_gis[gene] = 0.0

        else: 
            d_gis[gene] = merged.loc[gene]['0']

    print('Number of genes with no Reactome annotations: ', sum)
    return d_gis

def get_GIS_scores_go_reactome(args, genes_list,  notebook=False, dire_1=None, dire_2=None):

    if notebook:
        sys_df = load_dataset(dire_1, idx=0)
        gis_df = load_dataset(dire_2)
    else:
        sys_df = load_dataset(args.genes_ids_file_dir, idx=0)
        gis_df = load_dataset(args.go_reactome_gis_dir)
    print(type(gis_df.index[0]))
    print(type(sys_df.index[0]))
    merged = sys_df.merge(gis_df, right_on=gis_df.index, left_index=True)
    merged = merged.sort_values('0', ascending=False).drop_duplicates('symbol')
    merged = merged.set_index('symbol')
    # build gis dict
    d_gis = {}
    sum = 0
    for gene in tqdm.tqdm(genes_list):
        if gene not in merged.index:
            #if gene.split('.1')[0] in merged.index:
                
                #d_gis[gene] = merged.loc[gene.split('.1')[0]]['0']
            
            #else:
            #print(gene)
            sum = sum +1
            d_gis[gene] = 0.0

        else: 
            d_gis[gene] = merged.loc[gene]['0']

    print('Number of genes with no go o Reactome annotations: ', sum)
    return d_gis


def get_GIS_scores_go_reactome_hpo(args, genes_list,  notebook=False, dire_1=None, dire_2=None):

    if notebook:
        sys_df = load_dataset(dire_1, idx=0)
        gis_df = load_dataset(dire_2)
    else:
        sys_df = load_dataset(args.genes_ids_file_dir, idx=0)
        gis_df = load_dataset(args.go_reactome_hpo_gis_dir)
    print(type(gis_df.index[0]))
    print(type(sys_df.index[0]))
    merged = sys_df.merge(gis_df, right_on=gis_df.index, left_index=True)
    merged = merged.sort_values('0', ascending=False).drop_duplicates('symbol')
    merged = merged.set_index('symbol')
    # build gis dict
    d_gis = {}
    sum = 0
    for gene in tqdm.tqdm(genes_list):
        if gene not in merged.index:
            #if gene.split('.1')[0] in merged.index:
                
                #d_gis[gene] = merged.loc[gene.split('.1')[0]]['0']
            
            #else:
            #print(gene)
            sum = sum +1
            d_gis[gene] = 0.0

        else: 
            d_gis[gene] = merged.loc[gene]['0']

    print('Number of genes with no go or Reactome or hpo annotations: ', sum)
    return d_gis

def get_GIS_scores_hpo(args, genes_list,  notebook=False, dire_1=None, dire_2=None):

    if notebook:
        sys_df = load_dataset(dire_1, idx=0)
        gis_df = load_dataset(dire_2)
    else:
        sys_df = load_dataset(args.genes_ids_file_dir, idx=0)
        gis_df = load_dataset(args.hpo_gis_dir)
    print(type(gis_df.index[0]))
    print(type(sys_df.index[0]))
    merged = sys_df.merge(gis_df, right_on=gis_df.index, left_index=True)
    merged = merged.sort_values('0', ascending=False).drop_duplicates('symbol')
    merged = merged.set_index('symbol')
    # build gis dict
    d_gis = {}
    sum = 0
    for gene in tqdm.tqdm(genes_list):
        if gene not in merged.index:
            #if gene.split('.1')[0] in merged.index:
                
                #d_gis[gene] = merged.loc[gene.split('.1')[0]]['0']
            
            #else:
            #print(gene)
            sum = sum +1
            d_gis[gene] = 0.0

        else: 
            d_gis[gene] = merged.loc[gene]['0']

    print('Number of genes with no hpo annotations: ', sum)
    return d_gis

    
    
def compute_gis_lasso(score, k, v):      
    return (1/(1+ score**k))**v

def get_GIS_scores_lasso_func(args, genes_list, k=1, v=1, go=True, reactome=False, hpo=False, notebook=False, dire_1=None, dire_2=None):

    if go and reactome and hpo:
        gis_scores = get_GIS_scores_go_reactome_hpo(args, genes_list, notebook=notebook, dire_1=dire_1, dire_2=dire_2)
    elif go and reactome:
        gis_scores = get_GIS_scores_go_reactome(args, genes_list,  notebook=notebook, dire_1=dire_1, dire_2=dire_2)
    elif go: 
        gis_scores = get_GIS_scores(args, genes_list,  notebook=notebook, dire_1=dire_1, dire_2=dire_2)
    elif reactome:
        gis_scores = get_GIS_scores_reactome(args, genes_list,  notebook=notebook, dire_1=dire_1, dire_2=dire_2)
    elif go:
        gis_scores = get_GIS_scores_hpo(args,genes_list, notebook=notebook, dire_1=dire_1, dire_2=dire_2)
        
    gis_scores_lasso = [compute_gis_lasso(gis_scores[x], k, v) for x in gis_scores]
    return gis_scores_lasso

def get_pairs(a, b):
    pairs = []
    for i in a:
        for j in b:
            pairs.append([i, j])
    
    return pairs