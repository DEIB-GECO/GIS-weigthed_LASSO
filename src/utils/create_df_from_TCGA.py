import numpy as np
import pandas as pd
import json
import os 
import tqdm


dire_prior ='/home/mongardi/Metagene_repo/data/prior_knowledge'

def load_txt(filename):
    content = []
    with open(filename)as f:
        for line in f:
            content.append(line.strip())
    return content

def get_genes_ids(ids_list, files_dir, no_mirna=True):

    sample = ids_list[0]
    genes = []
    ids = []
    sample_file = os.path.join(files_dir,sample) + '.bed'
    with open(sample_file)as f:
        for line in f:
            genes.append(line.strip().split()[-5])
            ids.append(line.strip().split()[-6])
    ids_new = [None if x == 'gene' else x for x in ids]
    mirna_genes = [x for x in genes if 'MIR' in x.upper()]
    #mirna_genes_indexes = [i for i, x in enumerate(df.columns) if 'MIR' in x.upper()]
    #mirna_genes_ids = [x for x,y in zip(ids_new, genes) if 'MIR' in y.upper()]
    
    if no_mirna: 
        print('N of initial genes: ', len(genes))
        print('N of mirna genes removed: ', len(mirna_genes))
        genes = [x for x in genes if x not in mirna_genes]
        ids_new = [x for x,y in zip(ids_new, genes) if y not in mirna_genes]
    
    print('N of genes: ', len(genes))
    return genes, ids_new

def save_dict(dictionary,filename):
    with open(filename + ".txt", "w") as fp:
        json.dump(dictionary, fp)
        print("Done writing dict into .txt file")

def create_df_from_ids(ids_list, files_dir, no_mirna=False):
    
    store_samples = {}
    sample_ids = ids_list
    genes_all, gene_ids= get_genes_ids(ids_list, files_dir, no_mirna = False)
    genes, gene_ids= get_genes_ids(ids_list, files_dir, no_mirna = no_mirna)
    mask_genes = [g in genes for g in genes_all]
    print(sum(mask_genes))
    #mirna_genes = [x for x in df.columns if 'MIR' in x.upper()]
    #mirna_genes_indexes = [i for i, x in enumerate(df.columns) if 'MIR' in x.upper()]
    greedy_samples = [] 

    for sample in tqdm.tqdm(sample_ids):
        sample_file = os.path.join(files_dir,sample) + '.bed'
        expression_vector = []
        with open(sample_file)as f:
            for line in f:
                
                expression_vector.append(int(line.strip().split()[-3]))
                
        if len(expression_vector) > len(genes_all): 
            # remove greedy samples
            expression_vector = expression_vector[:len(genes_all)]
        
        if len(expression_vector) == len(genes_all):

            #d_expr = {g:e for g,e in zip(genes_all, expression_vector)}
            expression_vector = np.array(expression_vector)[mask_genes]
            is_greedy = is_greedy_sample(expression_vector)
            if is_greedy:
                greedy_samples.append(sample)
            else:  
                store_samples[sample] = expression_vector.tolist()
        
    print('N of initial samples: ', len(sample_ids))
    print('N of final samples: ', len(store_samples))
    print('N of greedy samples : ', len(greedy_samples))

    #print(store_samples)
    df = pd.DataFrame.from_dict(store_samples, orient='index',columns=genes)
    print(df.head())
    print(df.shape)
    return df, genes, gene_ids

def is_greedy_sample(vt, ratio=0.2):
    
    v = vt.copy()
    v_ratio = sum(v[::-1][:5])/sum(v)
    return v_ratio > ratio

def find_greedy_samples(df, ratio=0.2):

    ordered_genes = df.copy()
    ordered_genes.values.sort()
    ordered_genes = ordered_genes.iloc[:,::-1]
    tot_reads_per_sample = ordered_genes.sum(axis=1)
    top_5_genes_per_sample = ordered_genes.iloc[:, :5].sum(axis=1)
    ratio_geni_ingordi = top_5_genes_per_sample / tot_reads_per_sample
    samples_geni_ingordi =  ratio_geni_ingordi[ ratio_geni_ingordi > ratio].index
    print('N of greedy genes: ', len(samples_geni_ingordi))
    return samples_geni_ingordi

def find_not_expressed_genes(df, ratio=0.2):

    non_zero_entries_per_gene = df.astype(bool).sum(axis=0).sort_values()
    #print(df.shape[0])
    #print(non_zero_entries_per_gene[:20]/df.shape[0])
    #print(non_zero_entries_per_gene[-20:]/df.shape[0])
    expressed_genes = non_zero_entries_per_gene[non_zero_entries_per_gene/df.shape[0] >= ratio].index
    not_expressed_genes = non_zero_entries_per_gene[non_zero_entries_per_gene/df.shape[0] < ratio].index
    print('N of not expressed genes: ', len(not_expressed_genes))
    return not_expressed_genes

def find_not_expressed_genes_threshold(df, threshold = 1, ratio=0.2):

    non_zero_entries_per_gene = (df>=threshold).sum(axis=0).sort_values()
    #print(df.shape[0])
    #print(non_zero_entries_per_gene[:20]/df.shape[0])
    #print(non_zero_entries_per_gene[-20:]/df.shape[0])
    expressed_genes = non_zero_entries_per_gene[non_zero_entries_per_gene/df.shape[0] >= ratio].index
    not_expressed_genes = non_zero_entries_per_gene[non_zero_entries_per_gene/df.shape[0] < ratio].index
    print('N of not expressed genes: ', len(not_expressed_genes))
    return not_expressed_genes

def load_dataset(filename):

    df = pd.read_csv(filename)
    df_ind = df.set_index(df.columns[0])

    return df_ind

def shuffle(df, reset = False, seed=42):

    # reset only if idx are numbers
    if reset: 
        df = df.sample(frac=1).reset_index(drop=True)
    else:
         df = df.sample(frac=1)
    return df

def get_samples(df_labels, df_expr, balanced = True, seed=42, is_shuffle=True):
    
    #print(df_expr.head())
    #df_labels= df_labels.loc[df_expr.index]
    df_h = df_labels[df_labels[df_labels.columns[1]]== True]
    df_c = df_labels[df_labels[df_labels.columns[1]]== False]
    labels_counts = df_h.value_counts(df_h.columns[0])
    print(df_h.shape)
    if balanced:

        #print(df_h.value_counts(df_h.columns[0]))
        min_count = np.min(df_h.value_counts(df_h.columns[0]))
        df_h = df_h.groupby(df_h.columns[0], group_keys=False).apply(lambda x: x.sample(min_count, random_state=seed))
        df_c = df_c.groupby(df_c.columns[0], group_keys=False).apply(lambda x: x.sample(min_count, random_state=seed))
    
    df_final = pd.concat([df_h, df_c])
    
    #print(df_final.head())
    #print(df_final.shape)
    df_expr_final = df_expr.loc[df_final.index]
   
    df_final = pd.concat([df_final, df_expr_final], axis=1)
    if is_shuffle:
        df_final = shuffle(df_final)
    return df_final


def preprocessing(df_labels, files_dir, no_mirna= True,
                min_value = None, class_label = None, 
                balanced = True, seed=42, is_shuffle=False, coding= False, go_ann = False, is_RPM= True):



    print('shape:', df_labels.shape)
   
    if min_value is not None:
        
        class_counts = df_labels[df_labels[df_labels.columns[1]]==True].value_counts(df_labels.columns[0])
        print(class_counts)
        outputclasses = class_counts[class_counts >= min_value].index
        #samples = df_labels[df_labels[df_labels.columns[0]].isin(outputclasses) & df_labels[df_labels.columns[1]]==True].index
        samples = df_labels[df_labels[df_labels.columns[0]].isin(outputclasses)].index
        #samples = [x for x in samples if x in df_expr.index]
        #df_expr = df_expr.loc[samples]
        df_labels= df_labels.loc[samples]
        #print(df_expr.shape)

    if class_label is not None:
        #samples = df_labels[df_labels[df_labels.columns[0]] == class_label].index
        samples = df_labels[df_labels[df_labels.columns[0]].isin(class_label)].index
        #df_expr = df_expr.loc[samples]
        df_labels= df_labels.loc[samples]
    
    print(df_labels.shape)
    ids_list = df_labels.index.tolist()
    df_expr, _, _ = create_df_from_ids(ids_list, files_dir, no_mirna=no_mirna)

    print('df labels', df_labels.shape)
    print('----------Removing not expressed genes-----------')
    # remove genes not expressed in 80% of samples
    not_expressed_genes = find_not_expressed_genes_threshold(df_expr, threshold=4)
    df_reduced = df_expr.drop(not_expressed_genes, axis=1)
    print(df_reduced.shape)
    # get df with samples

    df_final = get_samples(df_labels, df_reduced , balanced = balanced, seed=seed, is_shuffle=is_shuffle)
    df_final = df_final.loc[:, ~df_final.columns.duplicated()]
    
    if coding:
        pc_filename = os.path.join(dire_prior, 'genes_coding.txt')
        pc_genes = load_txt(pc_filename)

        #intersection = list(set(df_final.columns) & set(pc_genes))
        intersection = [x for x in df_final.columns if x in pc_genes]

        df_final = pd.concat([df_final[df_final.columns[:2]], df_final[intersection]], axis=1)
        print('final df shape: ', df_final.shape)
        
    if go_ann:
        # to be fixed
        pc_filename = os.path.join(dire_prior, 'genes_coding_anns.txt')
        pc_genes = load_txt(pc_filename)

        #intersection = list(set(df_final.columns) & set(pc_genes))
        intersection = [x for x in pc_genes if x in df_final.columns]
        #intersection = [x for x in df_final.columns if x in pc_genes]
        print('intersection', len(intersection))
  
    
        
        df_final = pd.concat([df_final[df_final.columns[:2]], df_final[intersection]], axis=1)
        print('final df shape: ', df_final.shape)
        
    if is_RPM:
        print('RPM')
        #rpm = X / X.sum(axis=0) *1e6
        #rmp = np.log2(rpm + 1)
        
    return df_final


if __name__ == "__main__":

    import os 
    dire_data = '/home/mongardi/Metagene/Cancer/data'
    dire_results = '/home/mongardi/Metagene_repo/data/data_kidney'
    txt_file =  os.path.join(dire_data, 'samples_source_ids.txt')
    files_dir =  os.path.join(dire_data, 'files')
    labels_df_filename =  os.path.join(dire_data, 'tcga_classes.csv')
    df_labels = load_dataset(labels_df_filename)

    df_preprocessed = preprocessing(df_labels, files_dir , balanced = True,class_label= ['Kidney Renal Clear Cell Carcinoma'], coding=  True)
    print(df_preprocessed.shape)
    df_preprocessed.to_csv(os.path.join(dire_results,'Kidney_df_tr_coding_new.csv'))
 
