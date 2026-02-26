Notebook code extracts
=======================

The repository contains a set of Jupyter notebooks under `notebooks/ <https://github.com/wallet-maker/cytopus/tree/main/notebooks>`_ that
illustrate typical workflows.  Rather than merely listing their filenames, the
following sections reproduce the Python code cells directly so you can review
or copy the examples without opening the notebooks themselves.

`Cytopus_utils_tutorial.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/Cytopus_utils_tutorial.ipynb>`_
-------------------------------------------------------

.. code-block:: python

    # install and load the package
    %pip install cytopus

    # pygraphviz is required for plotting
    !apt install libgraphviz-dev
    !pip install pygraphviz

    # pyvis is required for plotting interactive graphs
    !pip install pyvis

    # install scanpy (required for loading example data used in this tutorial)
    !pip install scanpy

    # Import packages & load KnowledgeBase
    import cytopus as cp
    import pandas as pd
    import scanpy as sc

    G = cp.KnowledgeBase()

    # write files required to re-build cytopus KnowledgeBase to .csv
    cp.tl.hierarchy_to_csv(G.get_celltype_hierarchy(),filename='hierarchy.csv',
                            header_name=['Parent','Child'])
    cp.tl.geneset_to_csv(G.processes, filename='processes.csv',
                         header_name=['gene_set_name','gene_name'])
    cp.tl.geneset_to_csv(G.identities, filename='identities.csv',
                         header_name=['gene_set_name','gene_name'])
    cp.tl.metadata_to_csv(G.graph, 'metadata.csv', specific_class = False,
                           class_value=None)

    # Export gene sets from KnowledgeBase as .gmt files
    gp_dict = G.processes
    cell_dict = G.identities
    cp.label.get_gmt(cell_dict,save=True,path='cell_identities.gmt')
    cp.label.get_gmt(gp_dict,save=True,path='cellular_processes.gmt')

    # load example data
    adata = sc.read('data/adata_spectra.h5ad')
    adata.uns['SPECTRA_markers']

    # reload KnowledgeBase and compute cell type specificity
    G = cp.KnowledgeBase()
    cell_type_specificity = cp.label.get_celltype(adata,
                 celltype_key='cell_type_annotations',
                 factor_list=None,
                 Spectra_cell_scores='SPECTRA_cell_scores')

    # label marker genes and adjust indices
    overlap_df = cp.label.label_marker_genes(adata.uns['SPECTRA_markers'],
                                             G.processes, threshold=0.2)
    new_index = []
    for i,v in enumerate(overlap_df.index):
        new_index.append(cell_type_specificity[i]+'-X-'+str(v))
    overlap_df.index = new_index
    new_index = []
    for h,i in enumerate(overlap_df.index):
        new_index.append(str(h)+'-X-'+str(i))
    overlap_df.index = new_index

`Utils_tutorial.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/Utils_tutorial.ipynb>`_
--------------------------------------------

.. code-block:: python

    import cytopus as cp

    G = cp.KnowledgeBase()

    # write files required to re-build cytopus KnowledgeBase to .csv
    cp.tl.hierarchy_to_csv(G.get_celltype_hierarchy(),filename='hierarchy.csv',
                            header_name=['Parent','Child'])
    cp.tl.geneset_to_csv(G.processes, filename='processes.csv',
                         header_name=['gene_set_name','gene_name'])
    cp.tl.get_gmt(G.identities,save=True,path='cell_identities.gmt')
    cp.tl.get_gmt(G.processes,save=True,path='gene_programs.gmt')

`KnowledgeBase_queries_colaboratory.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/KnowledgeBase_queries_colaboratory.ipynb>`_
--------------------------------------------------------------------------------------------------

.. code-block:: python

    pip install cytopus
    # pygraphviz is required for plotting
    !apt install libgraphviz-dev
    !pip install pygraphviz
    # pyvis is required for plotting interactive graphs
    !pip install pyvis

    import cytopus as cp
    import IPython

    G = cp.KnowledgeBase()
    len(G.identities.keys())

    query_celltype = 'M'
    gs_of_interest = [x[0] for x in G.filter_edges(attribute_name='class',
                        attributes=['process_OF'], target=query_celltype)
                      if x[1]==query_celltype]
    print('Gene sets directly related to celltype',query_celltype,'are:',
          gs_of_interest)

    G.celltypes
    G.plot_celltypes()
    G.plot_graph_interactive(attributes=['cell_type','cellular_process'],
                              colors=['red','blue'], save_path='graph.html')
    !ls

    list(G.graph.successors('cDC'))
    list(G.graph.predecessors('cDC'))
    G.get_identities(['TNK'],include_subsets=True).keys()
    G.processes.keys()
    G.identities.keys()

    query_celltype = 'Treg'
    gs_of_interest = [x[0] for x in G.filter_edges(attribute_name='class',
                        attributes=['process_OF'], target=query_celltype)
                      if x[1]==query_celltype]
    print('Gene sets directly related to celltype',query_celltype,'are:',
          gs_of_interest)

    print('here external gene set from GO:')
    print(G.graph.nodes['all_macroautophagy_regulation_positive'])
    print('here manually curated gene set by the package author:')
    print(G.graph.nodes['Mac_LPS_response'])

    # dictionary conversion examples for Spectra (omitted for brevity)

`KnowledgeBase_construct.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/KnowledgeBase_construct.ipynb>`_
------------------------------------------------------------------------

.. code-block:: python

    import pandas as pd
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    from networkx.drawing.nx_agraph import graphviz_layout

    # path constants (set accordingly)
    DATA_DIR =  #path to directory containing .csv files for knowledge base
    gene_sets_path = DATA_DIR + '/Cytopus_1.31nc_gene-sets_x_genes.csv'
    metadata_path = DATA_DIR + '/Cytopus_1.31nc_versions_metadata.csv'
    cellular_hierarchies_path = DATA_DIR + '/Cytopus_1.31nc_hierarchies.csv'

    %pip install cytopus

    # prepare data
    cellular_hierarchies = pd.read_csv(cellular_hierarchies_path)
    celltype_edges = list(zip(list(cellular_hierarchies['child']),
                              list(cellular_hierarchies['parent'])))

    gene_sets = pd.read_csv(gene_sets_path)
    geneset_gene_edges = list(zip(list(gene_sets['gene_set_name']),
                                  list(gene_sets['gene_name'])))

    metadata = pd.read_csv(metadata_path,index_col='gene_set_name')
    geneset_celltype_edges = list(zip(list(metadata.index),
                                      list(metadata['cell_type_name'])))

    annotation_dict = metadata['annotation_name'].to_dict()
    metadata_columns = ['version_id', 'author', 'license',
           'license_link', 'license_type', 'gene_set_type', 'gene_set_topic']
    metadata_dict = metadata[metadata_columns].to_dict('index')

    # construct the KnowledgeBase
    import cytopus as cp
    G = cp.create.construct_kb(celltype_edges, geneset_gene_edges,
                   geneset_celltype_edges, annotation_dict,
                   metadata_dict=metadata_dict,
                   save=True, save_path=DATA_DIR+'Cytopus_1.31nc.txt')
    G
    G.plot_celltypes()

`Hierarchical_annotation_tutorial.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/Hierarchical_annotation_tutorial.ipynb>`_
------------------------------------------------------------------------------------------------------

.. code-block:: python

    # Install dependencies (example for Debian/Ubuntu)
    !apt install libgraphviz-dev
    !pip install pygraphviz
    !pip install cytopus
    !pip install scanpy
    !pip install networkx

    import networkx as nx
    import cytopus as cp
    import scanpy as sc
    from matplotlib import rcParams

    G = cp.KnowledgeBase()

    hierarchy_dict = cp.tl.hierarchy.get_hierarchy_dict(G)
    H = cp.tl.hierarchy.Hierarchy(hierarchy_dict)
    H.plot_celltypes(figsize=[30,30])

    import importlib.resources as pkg_resources
    import cytopus.data
    with pkg_resources.path(cytopus.data, 'adata_spectra.h5ad') as file_path:
        adata = sc.read_h5ad(file_path)

    adata.obs[['annotation_level_1', 'annotation_level_2', 'annotation_level_3']]
    H.add_cells(adata, obs_columns=['annotation_level_1',
                                     'annotation_level_2',
                                     'annotation_level_3'])
    H.query_ancestors(query_node='leukocyte', adata=adata,
                      obs_key='hierarchical_query')
    H.annotations.keys()
    adata.obs['hierarchical_query'].head()
    FIGSIZE = (5, 5)
    rcParams["figure.figsize"] = FIGSIZE
    sc.pl.umap(adata,color='hierarchical_query')

    H.query_ancestors(query_node='M', adata=adata, obs_key='M')
    H.query_ancestors(query_node='DC', adata=adata, obs_key='DC')
    sc.pl.umap(adata,color='M')
    sc.pl.umap(adata,color='DC')
    adata_myeloid = adata[~adata.obs['M'].isna()]
    adata_myeloid
