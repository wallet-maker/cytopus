from cytopus.knowledge_base import KnowledgeBase
import pandas as pd

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
    
    marker_genes: numpy.array or list of lists, factors x marker genes
    gs_label_dict: cytopus.KnowledgeBase or dict, with gene set names (str) as keys and gene sets (list) as values
    threshold: float, if overlap coefficient > than threshold the factor will be labeled with the gene set name with 
    maximum overlap coefficient
    
    returns: pandas.DataFrame, with overlap coefficients of factors (rows) and gene sets (columns), indices are relabeled 
    to the gene set with the maximum overlap coefficient
    '''
    import numpy as np

    if isinstance(gs_label_dict,KnowledgeBase):
        #collapse annotation dict
        gs_dict = {}
        key_list = []
        for key, value in gs_label_dict.celltype_process_dict.items():
            for k,v in value.items():
                if k not in key_list:
                    gs_dict[k]=v
                    key_list.append(k)
    elif isinstance(gs_label_dict, dict):
            for v in gs_label_dict.values():
                if isinstance(v,dict):
                    raise ValueError('gs_label_dict is a nested dictionary. gs_label_dict must be a flat/non-nested dictionary with gene set names as keys (str) amd gene sets (lists of strings) as values')
            gs_dict = gs_label_dict
    else:
        raise ValueError('gs_label_dict must be a dictionary or a cytopus.kb.queries.KnowledgeBase object')

    overlap_df = pd.DataFrame()
    for i, v in pd.DataFrame(marker_genes).T.items():
        overlap_temp = []
        gs_names_temp = []
        for gs_name, gs in gs_dict.items():
            gene_set = set(gs)
            marker_set = set(v)
            #check and remove for nans
            if 'nan' in gene_set:
                gene_set.remove('nan')
            if 'nan' in marker_set:
                marker_set.remove('nan')
            if len(gene_set) > 0 and len(marker_set)>0:
                overlap_temp.append(overlap_coefficient(set(gene_set),set(marker_set)))
            else:
                overlap_temp.append(np.nan)
            gs_names_temp.append(gs_name)
        overlap_df_temp = pd.DataFrame(overlap_temp, columns=[i],index=gs_names_temp).T
        overlap_df = pd.concat([overlap_df,overlap_df_temp])
    marker_gene_labels = [] #gene sets
    for marker_set in overlap_df.index:
        max_overlap = overlap_df.loc[marker_set].sort_values().index[-1]
        if overlap_df.loc[marker_set].sort_values().values[-1] >threshold:
            marker_gene_labels.append(max_overlap)
        else:
            marker_gene_labels.append(marker_set)
    overlap_df.index = marker_gene_labels     
        
    return overlap_df


def get_celltype(adata, celltype_key,factor_list=None,Spectra_cell_scores= 'SPECTRA_cell_scores'):
    '''
    For a list of factors check in which cell types they are expressed
    adata: anndata.AnnData, containing cell type labels in adata.obs[celltype_key]
    celltype_key: str, key for adata.obs containing the cell type labels
    factor_list: list, list of keys for factor loadings in .obs, if none use factor loadings in adata.obsm['SPECTRA_factors']
    return: dictionary mapping factor names and celltypes
    Spectra_cell_scores: str, key for Spectra cell scores in adata.obsm
    '''
    
    if factor_list!= None:
        factors= adata.obs[factor_list]
        factors['celltype'] = list(adata.obs[celltype_key])
    else:
        factors = pd.DataFrame(adata.obsm[Spectra_cell_scores])
        factors['celltype'] = list(adata.obs[celltype_key])
    
    #create factor:celltype dict
    grouped_df = factors.groupby('celltype').mean()
    #get factor names for global (expressed in all cells) and cell type spec factors
    global_factor_names = grouped_df.T[(grouped_df!=0).all()].index
    specific_factor_names= [x for x in grouped_df.columns if x not in global_factor_names]
    #add global factors to dict
    factor_names_global = {x:'global' for x in global_factor_names}

    #get celltype for celltype spec factors
    grouped_df_spec = grouped_df[specific_factor_names]

    for i in grouped_df_spec.columns:
        factor_names_global[i] = grouped_df_spec[i].sort_values(ascending=False).index[0]
    return factor_names_global
    

def get_gmt(gs_dict,save=False,path=None):
    '''
    transform a dictionary into a .gmt file
    gs_dict: dict, gene set dictionary with format {'gene set name':['Gene_a','Gene_b','Gene_c',...]}
    save: bool, if True saves .gmt file to path
    path: str, path to save .gmt file
    '''
    import numpy as np
    import pandas as pd
    #retrieve all genes from dict
    genes = []
    for k,v in gs_dict.items():
        genes = genes+v
    genes = list(set(genes))
    
    #pad the lists in gs_dict to equal lengths
    max_length = max(map(len, gs_dict.values()))

    for k,v in gs_dict.items():
        if len(v)<max_length:
            gs_dict[k]+= [np.nan]*(max_length-len(v))

    #transform into df
    gs_df = pd.DataFrame(gs_dict).T
    
    if save:
        gs_df.to_csv(path,sep='\t',header=False)
        print('print saving to:',path) 
    else:
        return gs_df
    
    

