import cytopus
import pandas as pd
import numpy as np

def overlap_coefficient(set_a,set_b):
    '''
    calculate the overlap coefficient between two sets
    '''
    min_len = min([len(set_a),len(set_b)])
    intersect_len = len(set_a.intersection(set_b))
    overlap = intersect_len/min_len
    return overlap

def label_marker_genes(marker_genes, gs_label_dict, threshold = 0.4):
    '''
    label an array of marker genes using a KnowledgeBase or a dictionary derived from the KnowledgeBase
    returns a dataframe of overlap coefficients for each gene set annotation and marker gene
    marker_genes: array factors x marker genes or a KnowledgeBase object
    label an array containing marker genes by its overlap with a dictionary of gene sets from the knowledge base:
    KnowledgeBase.celltype_process_dict
    '''
    import numpy as np
    import warnings

    if isinstance(gs_label_dict,cytopus.kb.kb_queries.KnowledgeBase):
        #collapse annotation dict
        gs_dict = {}
        key_list = []
        for key, value in gs_label_dict.celltype_process_dict.items():
            for k,v in value.items():
                if k not in key_list:
                    gs_dict[k]=v
                    key_list.append(k)
    else:
        gs_dict = gs_label_dict
    if threshold <=0.5:
        warnings.warn('threshold <= 0.5 can resolve in ties, consider using higher threshold or manually check for ties in  overlap_df output')
    overlap_df = pd.DataFrame()
    for i, v in pd.DataFrame(marker_genes).T.iteritems():
        for gs_name, gs in gs_dict.items():
            overlap_df.loc[i,gs_name] =  overlap_coefficient(set(gs),set(v))
    marker_gene_labels = [] #gene sets
    for marker_set in overlap_df.index:
        max_overlap = overlap_df.loc[marker_set].sort_values().index[-1]
        if overlap_df.loc[marker_set].sort_values().values[-1] >threshold:
            marker_gene_labels.append(max_overlap)
        else:
            marker_gene_labels.append(marker_set)
    overlap_df.index = marker_gene_labels     
        
    return overlap_df


#include function to label leukocytes