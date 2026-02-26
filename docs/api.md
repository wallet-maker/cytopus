# API & Core Logic Reference

This page documents the public interfaces and core functions available in the repository (derived directly from the code).

## Package entry
- `import cytopus as cp`
- Exposed symbols (from `cytopus/__init__.py` and subpackages):
  - `cp.KnowledgeBase` (class)
  - `cp.get_data` (function)
  - `cp.tl` subpackage with `label`, `create`, and `hierarchy` modules

## `get_data(filename)`
- Location: `cytopus/knowledge_base/kb_queries.py`
- Signature: `get_data(filename)`
- Purpose: Return the filesystem path to a packaged data file using `pkg_resources.resource_filename('cytopus', 'data/' + filename)`.

## `KnowledgeBase` class (core)
- Location: `cytopus/knowledge_base/kb_queries.py`
- Constructor signature: `KnowledgeBase(graph=None)`
  - `graph` can be:
    - `None` (loads default packaged file `Cytopus_1.31nc.txt` via `get_data`),
    - a `networkx.DiGraph` instance,
    - or a string path to a pickled graph (unpickled with `pickle.load`).

- Important attributes set during initialization:
  - `self.graph` — NetworkX DiGraph containing nodes and edges.
  - `self.celltypes` — list of node names where node attribute `class` == `'cell_type'`.
  - `self.processes` — dict mapping process gene-set names to lists of genes (cleaned of NaN values).
  - `self.identities` — dict mapping identity gene-set names to lists of genes.

- Key methods:
  - `filter_nodes(attributes, attribute_name=None, origin=None, target=None)` — filter nodes by node attribute values.
  - `filter_edges(attributes, attribute_name=None, origin=None, target=None)` — filter edges by edge attribute values.
  - `get_processes(gene_sets)` — return a dict mapping gene_set -> genes for specified gene_sets by reading `'gene_OF'` edges.
  - `get_celltype_hierarchy(node='all-cells', invert=False)` — wrapper around `extract_hierarchy` to return nested dict hierarchy rooted at `node`.
  - `get_celltype_processes(celltypes, global_celltypes=[None], get_parents=True, get_children=True, parent_depth=1, child_depth=None, fill_missing=True, parent_depth_dict=None, child_depth_dict=None, inplace=True)` — core query returning nested dict mapping queried cell types to available process gene sets and their gene members. If `inplace=True`, stores result to `self.celltype_process_dict`.
  - `get_identities(celltypes_identities, include_subsets=False)` — return identity gene sets per cell type.
  - `plot_celltypes(figure_size=[30,30], node_size=1000, edge_width=1, arrow_size=20, edge_color='k', node_color='#8decf5', label_size=20)` — plot cell-type subgraph using Graphviz layout and matplotlib.
  - `plot_graph_interactive(attributes=['cell_type','cellular_process'], colors=['red','blue'], save_path='graph.html')` — export interactive HTML using `pyvis`.

## `cytopus.tl.create.construct_kb`
- Location: `cytopus/tl/create.py`
- Signature: `construct_kb(celltype_edges, geneset_gene_edges, geneset_celltype_edges, annotation_dict, metadata_dict=None, save=False, save_path=None)`
- Purpose: Build a `networkx.DiGraph` from lists of edges and annotations and return `KnowledgeBase(graph=G)`.
- Important behavior:
  - Expects `annotation_dict` mapping gene set name -> `'cellular_process'` or `'cellular_identity'`.
  - Adds edge attribute `'class'` with values such as `'gene_OF'` and `'SUBSET_OF'` for use by filters.
  - Optionally pickles the constructed graph if `save=True` and `save_path` is provided.

## `cytopus.tl.hierarchy` module
- Location: `cytopus/tl/hierarchy.py`
- Key functions:
  - `get_hierarchy_dict(G)` — build nested dict for the cell-type hierarchy from a `KnowledgeBase` object `G`.
  - `create_hierarchical_graph(data, type_label)` — build a NetworkX DiGraph from a nested dict `data` and set node attribute `type`.
  - `get_nodes_of_type(graph, node_type)` — return nodes with `graph.nodes[node]['type'] == node_type`.
- `Hierarchy` class:
  - Constructor: `Hierarchy(hierarchy_dict)` — build internal graph and sets node types.
  - `plot_celltypes(...)` — plot cell types using Graphviz layout.
  - `add_cells(adata, obs_columns=None)` — attach cell barcodes from an AnnData object to the hierarchical graph.
  - `query_ancestors(query_node, adata=None, obs_key='hierarchical_query')` — retrieve barcodes belonging to `query_node` and its subsets; optionally annotate `adata.obs[obs_key]`.

## `cytopus.tl.label` module
- Location: `cytopus/tl/label.py`
- Key functions:
  - `overlap_coefficient(set_a, set_b)` — returns |A ∩ B| / min(|A|, |B|).
  - `label_marker_genes(marker_genes, gs_label_dict, threshold=0.4)` — label marker gene lists against gene sets contained in a `KnowledgeBase` or a dict mapping gene sets to gene lists. Returns a `pandas.DataFrame` of overlap coefficients per marker set vs gene set, with marker set indices relabeled to the gene-set with maximum overlap above `threshold`.
  - `get_celltype(adata, celltype_key, factor_list=None, Spectra_cell_scores='SPECTRA_cell_scores')` — map factor loadings to cell types using mean factor scores per cell type.
  - `get_gmt(gs_dict, save=False, path=None)` — transform a gene-set dict to a GMT-like `pandas.DataFrame` or save to file.
  - CSV export helpers: `flatten_hierarchical_dict`, `hierarchy_to_csv`, `geneset_to_csv`, `metadata_to_csv`.

## Notes about API completeness
- The repository does not provide formal docstrings for all parameters; the above signatures and behaviors are derived from source code bodies and README usage examples.
- Some functions expect `pandas`, `numpy`, and `anndata` objects to be available at runtime; where they are referenced, those libraries are required to use the related functions.
