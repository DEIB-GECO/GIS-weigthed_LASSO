def get_args():
  
    import argparse
    import os

    default_dire = 'data/prior_knowledge'
    default_dire_ids = 'data/prior_knowledge'
    parser = argparse.ArgumentParser()

    # prior knowledge
    parser.add_argument('--go_file_dir', type=str, default=os.path.join(default_dire, 'go_dag.txt'))
    parser.add_argument('--go_alt_file_dir', type=str, default= os.path.join(default_dire, 'go_alternative_ids.txt'))
    parser.add_argument('--go_obo_file_dir', type=str, default= os.path.join(default_dire, 'go.obo'))
    parser.add_argument('--genes_ids_file_dir', type=str, default= os.path.join(default_dire_ids, 'genes_and_ids_all_red.csv'))
    parser.add_argument('--go_annotations_file_dir', type=str, default= os.path.join(default_dire, 'go_annotations_new.csv'))
    parser.add_argument('--kegg_gene_pathways_file_dir', type=str, default= os.path.join(default_dire, 'kegg_gene_pathways_names.csv'))
    parser.add_argument('--kegg_hsa_pathways_file_dir', type=str, default= os.path.join(default_dire, 'kegg_hsa_pathways.txt'))
    parser.add_argument('--reactome_hsa_pathways_file_dir', type=str, default= os.path.join(default_dire, 'genes_to_pathways_hsa.csv'))
    parser.add_argument('--reactome_relations_file_dir', type=str, default= os.path.join(default_dire, 'ReactomePathwaysRelation.txt'))
    parser.add_argument('--hpo_file_dir', type=str, default= os.path.join(default_dire, 'hpo_dag.txt'))
    parser.add_argument('--hpo_alt_file_dir', type=str, default= os.path.join(default_dire, 'hpo_alternative_ids.txt'))
    parser.add_argument('--hpo_annotations_file_dir', type=str, default= os.path.join(default_dire, 'genes_to_phenotype_red.txt'))
    parser.add_argument('--hpo_obo_file_dir', type=str, default=os.path.join(default_dire, 'hp.obo'))
    
    # folder where to save the update prior knowledge
    parser.add_argument('--save_files_dir', type=str, default= 'data/prior_knowledge_updated')
    
    parser.add_argument('--IC_go_dir', type=str, default= os.path.join(default_dire, 'IC_scores_normalized_go.csv'))
    
    #GIS
    parser.add_argument('--go_gis_dir', type=str, default= os.path.join(default_dire, 'gene_scores_norm_go.csv'))
    parser.add_argument('--reactome_gis_dir', type=str, default= os.path.join(default_dire, 'gene_scores_norm_reactome.csv'))
    parser.add_argument('--hpo_gis_dir', type=str, default= os.path.join(default_dire, 'gene_scores_norm_hpo.csv'))
    parser.add_argument('--go_reactome_gis_dir', type=str, default= os.path.join(default_dire, 'gene_scores_norm_go_reactome.csv'))
    parser.add_argument('--go_hpo_gis_dir', type=str, default= os.path.join(default_dire, 'gene_scores_norm_go_hpo.csv'))
    parser.add_argument('--go_reactome_hpo_gis_dir', type=str, default= os.path.join(default_dire, 'gene_scores_norm_go_reactome_hpo.csv'))

    # datasets 
    parser.add_argument('--BRCA_dataset_dir', type=str, default= 'data/data_BRCA/BRCA_dataset.csv')
    parser.add_argument('--BRCA_dataset_index_col', type=str, default='sample_id')
    parser.add_argument('--BRCA_results_dir', type=str, default= 'results/BRCA')
    parser.add_argument('--CRC_dataset_dir', type=str, default= 'data/data_CRC/CRC_dataset.csv')
    parser.add_argument('--CRC_dataset_index_col', type=str, default='Unnamed: 0')
    parser.add_argument('--CRC_results_dir', type=str, default= 'results/CRC')

    args = parser.parse_args()
    return args

