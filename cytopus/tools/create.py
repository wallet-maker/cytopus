from cytopus.knowledge_base import KnowledgeBase
import networkx as nx

def construct_kb(celltype_edges, geneset_gene_edges,geneset_celltype_edges,annotation_dict,metadata_dict=None,save=False, save_path=None):
    '''
    construct a cytopus.kb.KnowledgeBase object
    celltype_edges: list, list of tuples storing the edges of the cell type hierarchy as ('child', 'parent')
    geneset_gene_edges: list, list of tuples storing the edges connecting every gene_set with every gene as ('gene_set','gene')
    geneset_celltype_edges: list, list of tuples storing the edges connecting every gene sets with its cell type as ('gene_set','celltype')
    annotation_dict: dict, containing the gene set names as keys and their annotation names (cellular_process or cellular_identity) as values
    metadata_dict: dict, nested dict containing the gene set names as keys and a dict storing their attributes_categories as keys and corresponding attributes as values
    save: bool, if True saves the data to the path provided in save_path
    save_path: str, path to save the data to (.txt file)
    '''

    #get genes, genesets, celltypes
    genes = list(set([x[1] for x in geneset_gene_edges]))
    genes = [(x,{'class':'gene'}) for x in genes]
    gene_sets = list(set([x[0] for x in geneset_gene_edges]))
    celltypes = list(set([x[0] for x in celltype_edges]).union(set([x[1] for x in celltype_edges])))
    celltypes = [(x,{'class':'cell_type'})for x in celltypes]

    #some sanity checks
    celltypes_in_hierarchy = set([x[0] for x in celltypes])
    celltypes_of_genesets = set([x[1] for x in geneset_celltype_edges])
    set_dif = celltypes_of_genesets - celltypes_in_hierarchy
    if  set_dif != set():
        print('WARNING: missing cell types:',set_dif,'in the cell type hierarchy. Please append cell type hierarchy.')
    else:
        print('all cell types in gene set are contained in the cell type hierarchy')
    
    genesets_in_celltype_edges = set([x[0] for x in geneset_celltype_edges])
    genesets_in_gene_edges = set([x[0] for x in geneset_gene_edges])

    if  genesets_in_celltype_edges != genesets_in_gene_edges:
        print('WARNING: Gene sets in geneset_celltype_edges and geneset_gene_edges are not identical')

    #set edge attributes (important for queries)
    geneset_gene_edges = [x + ({'class':'gene_OF'},) for x in geneset_gene_edges]
    celltype_edges = [x + ({'class':'SUBSET_OF'},) for x in celltype_edges]


    #sort processes and identities
    processes = []
    identities = []

    for i in gene_sets:
        if annotation_dict[i] == 'cellular_process':
            processes.append(i)
        elif annotation_dict[i] == 'cellular_identity':
            identities.append(i)
        else:
            raise(ValueError('all gene sets annotation names should be either cellular_process or cellular_identity'))

    geneset_gene_edges_processes = [x for x in geneset_gene_edges if x[0] in processes]
    geneset_gene_edges_identities = [x for x in geneset_gene_edges if x[0] in identities]
    geneset_celltype_edge_processes = [x + ({'class':'process_OF'},) for x in geneset_celltype_edges if x[0] in processes]
    geneset_celltype_edge_identities = [x + ({'class':'identity_OF'},) for x in geneset_celltype_edges if x[0] in identities]
        
    #construct graph
    G = nx.DiGraph()
    G.add_nodes_from(genes)
    G.add_nodes_from(gene_sets)
    G.add_nodes_from(identities)
    G.add_nodes_from(celltypes)
    G.add_edges_from(geneset_gene_edges_processes)
    G.add_edges_from(geneset_gene_edges_identities)
    G.add_edges_from(celltype_edges)
    G.add_edges_from(geneset_celltype_edge_processes)
    G.add_edges_from(geneset_celltype_edge_identities)

    #set node metadata
    if isinstance(metadata_dict,dict):
        nx.set_node_attributes(G, metadata_dict)
    else:
        print('No metadata dictionary provided (optional), skipping metadata assignment.')
    if save:
        if not isinstance(save_path,str):
            print('WARNING: Please provide save_path if you want to save the data. Skipping saving step.')
        else:
            import pickle
            with open(save_path, 'wb') as f:
                pickle.dump(G, f)
            print('Pickled and saved to:',save_path)
    return KnowledgeBase(graph=G)
    
