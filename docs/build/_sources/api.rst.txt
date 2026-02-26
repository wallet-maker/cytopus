API & Core Logic Reference
==========================

This page documents the public interfaces and core functions available in the repository (derived directly from the code).

Package entry
-------------

- `import cytopus as cp`
- Exposed symbols (from `cytopus/__init__.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/__init__.py>`_ and subpackages):
  * `cp.KnowledgeBase` (class)
  * `cp.get_data` (function)
  * `cp.tl` subpackage with `label`, `create`, and `hierarchy` modules

`get_data(filename)`
---------------------

- Location: `cytopus/knowledge_base/kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_
- Signature: `get_data(filename)`
- Purpose: Return the filesystem path to a packaged data file using `pkg_resources.resource_filename('cytopus', 'data/' + filename)`.

`KnowledgeBase` class (core)
----------------------------

- Location: `cytopus/knowledge_base/kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_
- Constructor signature: `KnowledgeBase(graph=None)`

  * `graph` can be one of the following:

    - ``None`` (loads default packaged file ``Cytopus_1.31nc.txt`` via ``get_data``),
    - a ``networkx.DiGraph`` instance,
    - or a string path to a pickled graph (unpickled with ``pickle.load``).

- Important attributes set during initialization:
  * `self.graph` — NetworkX DiGraph containing nodes and edges.
  * `self.celltypes` — list of node names where node attribute `class` == `'cell_type'`.
  * `self.processes` — dict mapping process gene-set names to lists of genes (cleaned of NaN values).
  * `self.identities` — dict mapping identity gene-set names to lists of genes.

- Key methods:
  * `filter_nodes(attributes, attribute_name=None, origin=None, target=None)` — filter nodes by node attribute values.
  * `filter_edges(attributes, attribute_name=None, origin=None, target=None)` — filter edges by edge attribute values.
  * `get_processes(gene_sets)` — return a dict mapping gene_set -> genes for specified gene_sets by reading `'gene_OF'` edges.
  * `get_celltype_hierarchy(node='all-cells', invert=False)` — wrapper around `extract_hierarchy` to return nested dict hierarchy rooted at `node`.
  * `get_celltype_processes(celltypes, global_celltypes=[None], get_parents=True, get_children=True, parent_depth=1, child_depth=None, fill_missing=True, parent_depth_dict=None, child_depth_dict=None, inplace=True)` — core query returning nested dict mapping queried cell types to available process gene sets and their gene members. If `inplace=True`, stores result to `self.celltype_process_dict`.
  * `get_identities(celltypes_identities, include_subsets=False)` — return identity gene sets per cell type.
  * `plot_celltypes(figure_size=[30,30], node_size=1000, edge_width=1, arrow_size=20, edge_color='k', node_color='#8decf5', label_size=20)` — plot cell-type subgraph using Graphviz layout and matplotlib.
  * `plot_graph_interactive(attributes=['cell_type','cellular_process'], colors=['red','blue'], save_path='graph.html')` — export interactive HTML using `pyvis`.

`cytopus.tl.create.construct_kb`
================================

- Location: `cytopus/tl/create.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/create.py>`_
- Signature: `construct_kb(celltype_edges, geneset_gene_edges, geneset_celltype_edges, annotation_dict, metadata_dict=None, save=False, save_path=None)`
- Purpose: Build a `networkx.DiGraph` from lists of edges and annotations and return `KnowledgeBase(graph=G)`.
- Important behavior:
  * Expects `annotation_dict` mapping gene set name -> `'cellular_process'` or `'cellular_identity'`.
  * Adds edge attribute `'class'` with values such as `'gene_OF'` and `'SUBSET_OF'` for use by filters.
  * Optionally pickles the constructed graph if `save=True` and `save_path` is provided.

`cytopus.tl.hierarchy` module
================================

- Location: `cytopus/tl/hierarchy.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/hierarchy.py>`_
- Key functions:
  * `get_hierarchy_dict(G)` — build nested dict for the cell-type hierarchy from a `KnowledgeBase` object `G`.
  * `create_hierarchical_graph(data, type_label)` — build a NetworkX DiGraph from a nested dict `data` and set node attribute `type`.
  * `get_nodes_of_type(graph, node_type)` — return nodes with `graph.nodes[node]['type'] == node_type`.
- `Hierarchy` class:
  * Constructor: `Hierarchy(hierarchy_dict)` — build internal graph and sets node types.
  * `plot_celltypes(...)` — plot cell types using Graphviz layout.
  * `add_cells(adata, obs_columns=None)` — attach cell barcodes from an AnnData object to the hierarchical graph.
  * `query_ancestors(query_node, adata=None, obs_key='hierarchical_query')` — retrieve barcodes belonging to `query_node` and its subsets; optionally annotate `adata.obs[obs_key]`.

`cytopus.tl.label` module
==========================

- Location: `cytopus/tl/label.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/label.py>`_
- Key functions:
  * `overlap_coefficient(set_a, set_b)`

    - Purpose: compute the overlap coefficient between two sets.
    - Signature: ``overlap_coefficient(set_a, set_b)``
    - Returns: float :math:`\frac{|A \cap B|}{\min(|A|,|B|)}`

  * `label_marker_genes(marker_genes, gs_label_dict, threshold=0.4)`

    - Purpose: Label arrays of marker genes using the KnowledgeBase or a flat dictionary of gene sets.
    - Signature: ``label_marker_genes(marker_genes, gs_label_dict, threshold=0.4)``
    - Parameters:
      - ``marker_genes``: list or array-like (factors x marker genes). Each column (or list) corresponds to a factor/marker set.
      - ``gs_label_dict``: either a `KnowledgeBase` object or a flat dict mapping gene set names to lists of genes.
      - ``threshold``: float threshold for assigning a label based on maximum overlap (default 0.4).
    - Returns: a `pandas.DataFrame` where rows correspond to marker sets and columns to gene sets; indices are relabeled to the assigned gene set when overlap > threshold.

  * `get_celltype(adata, celltype_key, factor_list=None, Spectra_cell_scores='SPECTRA_cell_scores')`

    - Purpose: Determine which cell types express each factor (Spectra-style factor scores).
    - Signature: ``get_celltype(adata, celltype_key, factor_list=None, Spectra_cell_scores='SPECTRA_cell_scores')``
    - Notes: If ``factor_list`` is provided, factor values are read from ``adata.obs``; otherwise values are taken from ``adata.obsm[Spectra_cell_scores]``. Returns a dict mapping factor names to either `'global'` or a specific cell type.

  * `get_gmt(gs_dict, save=False, path=None)`

    - Purpose: Convert a gene-set dictionary to a GMT-like table and optionally save it.
    - Signature: ``get_gmt(gs_dict, save=False, path=None)``

  * CSV export helpers: ``flatten_hierarchical_dict``, ``hierarchy_to_csv``, ``geneset_to_csv``, ``metadata_to_csv`` — utilities that write hierarchy and geneset content to CSV files.

Notes about API completeness
-----------------------------

- The repository does not provide formal docstrings for all parameters; the above signatures and behaviors are derived from source code bodies and README usage examples.
- Some functions expect `pandas`, `numpy`, and `anndata` objects to be available at runtime; where they are referenced, those libraries are required to use the related functions.

Detailed function & class reference
===================================

Below are expanded, developer-oriented references for the primary modules. These are written from the repository source and include signatures, parameter descriptions, return values, and small usage examples.

`KnowledgeBase` (expanded)
---------------------------

Location: `cytopus/knowledge_base/kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_

Class: ``KnowledgeBase``

Constructor
~~~~~~~~~~~

``KnowledgeBase(graph=None)``

- ``graph``: ``None``, a path string to a pickled NetworkX DiGraph, or a ``networkx.DiGraph`` instance. If ``None``, the constructor loads the packaged default KB file ``Cytopus_1.31nc.txt`` via ``get_data``.

On initialization the object sets the following attributes:

- ``self.graph`` — the underlying ``networkx.DiGraph``.
- ``self.celltypes`` — list of nodes where node attribute ``class == 'cell_type'`` (retrieved via ``filter_nodes``).
- ``self.processes`` — dict mapping process gene set names to lists of genes (cleaned of NaNs).
- ``self.identities`` — dict mapping identity gene set names to lists of genes.

Important methods
~~~~~~~~~~~~~~~~~

- ``filter_nodes(attributes, attribute_name=None, origin=None, target=None)``
  - Filters nodes by node attribute values. Returns list of node keys.

- ``filter_edges(attributes, attribute_name=None, origin=None, target=None)``
  - Filters edges by edge attribute values. Returns list of edge tuples.

- ``get_processes(gene_sets)``
  - Input: list of gene set names
  - Returns: dict mapping gene set name -> list of gene names (reads edges marked with edge attribute ``class == 'gene_OF'``)

- ``get_identities(celltypes_identities, include_subsets=False)``
  - Input: list of cell type names
  - Returns: dict mapping cell type -> identity gene set (reads ``identity_OF`` edges)

- ``get_celltype_hierarchy(node='all-cells', invert=False)``
  - Returns a nested dictionary representation of the hierarchy rooted at ``node`` (uses ``extract_hierarchy``).

- ``get_celltype_processes(celltypes, global_celltypes=[None], get_parents=True, get_children=True, parent_depth=1, child_depth=None, fill_missing=True, parent_depth_dict=None, child_depth_dict=None, inplace=True)``
  - Core query that returns (or stores in ``self.celltype_process_dict``) a nested dict mapping requested cell types to their process gene sets and member genes. This method:
    - Builds a subgraph limited to known cell types
    - Traverses parents/children using BFS (controlled by ``parent_depth`` and ``child_depth``)
    - Collects edges labeled ``process_OF`` to find which gene sets connect to the discovered cell types
    - Maps gene sets -> genes via ``get_processes`` and constructs a nested celltype -> {geneset: [genes]} dict

- ``plot_celltypes(...)`` and ``plot_graph_interactive(...)``
  - Visualization helpers. ``plot_celltypes`` uses Graphviz layout via ``networkx.drawing.nx_agraph.graphviz_layout`` (requires ``pygraphviz``/Graphviz). ``plot_graph_interactive`` uses ``pyvis`` to export an interactive HTML graph.

Example usage
~~~~~~~~~~~~~

````python
import cytopus as cp
G = cp.KnowledgeBase()           # loads packaged KB by default
print(len(G.celltypes))          # number of cell type nodes
G.get_celltype_processes(['B'])  # populate G.celltype_process_dict
print(G.celltype_process_dict['B'])
````

`tl.create.construct_kb` (expanded)
------------------------------------

Location: `cytopus/tl/create.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/create.py>`_

``construct_kb(celltype_edges, geneset_gene_edges, geneset_celltype_edges, annotation_dict, metadata_dict=None, save=False, save_path=None)``

- ``celltype_edges``: list of tuples ``('child', 'parent')`` representing the cell type hierarchy.
- ``geneset_gene_edges``: list of tuples ``('gene_set','gene')`` linking gene sets to genes.
- ``geneset_celltype_edges``: list of tuples ``('gene_set','celltype')`` linking gene sets to their cell types.
- ``annotation_dict``: dict mapping gene set name -> annotation type string; allowed values: ``'cellular_process'`` or ``'cellular_identity'``.
- ``metadata_dict``: optional dict mapping node names to attribute dicts.
- ``save``/``save_path``: optionally pickle the constructed NetworkX graph to disk.

Behavior
~~~~~~~~~
- Adds node attributes (``class``) for genes and cell types.
- Adds edge attributes (``class``) labeling edges as ``'gene_OF'``, ``'process_OF'``, ``'identity_OF'``, and ``'SUBSET_OF'`` to encode semantics.
- Returns a ``KnowledgeBase(graph=G)`` wrapping the constructed graph.

Example usage
~~~~~~~~~~~~~

````python
from cytopus.tl.create import construct_kb
# Provide lists of edges and annotation_dict per the function signature
# then:
kb = construct_kb(celltype_edges, geneset_gene_edges, geneset_celltype_edges, annotation_dict)
````

`tl.hierarchy` (expanded)
-------------------------

Location: `cytopus/tl/hierarchy.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/hierarchy.py>`_

Key utilities:

- ``get_hierarchy_dict(G)``: build a nested dictionary representing the cell-type hierarchy from a ``KnowledgeBase`` instance (filters nodes where ``class == 'cell_type'``).
- ``create_hierarchical_graph(data, type_label)``: create a NetworkX DiGraph from a nested dict representation.
- ``get_nodes_of_type(graph, node_type)`` and ``get_node_labels(graph, node_type)``: helpers for traversal and ordering.

`Hierarchy` class
~~~~~~~~~~~~~~~~~

``Hierarchy(hierarchy_dict)``

- Constructor builds an internal ``networkx.DiGraph`` with node attribute ``type='cell_type'``.
- ``add_cells(adata, obs_columns=None)``: attach cell barcodes from an ``anndata.AnnData`` object into the hierarchical graph by reading annotations from ``adata.obs`` and adding nodes of type ``'cell'`` and edges from the most granular cell-type node to the barcode node.
- ``query_ancestors(query_node, adata=None, obs_key='hierarchical_query')``: for a given cell-type node, collect barcodes for that cell type and its subsets; optionally annotate ``adata.obs[obs_key]`` with assigned labels.

Example usage
~~~~~~~~~~~~~

````python
from cytopus.tl.hierarchy import Hierarchy, get_hierarchy_dict
G = cp.KnowledgeBase()
hier = Hierarchy(get_hierarchy_dict(G))
# Attach cells from an AnnData object (requires `anndata` package)
hier.add_cells(adata, obs_columns=['annotation'])
hier.query_ancestors('B', adata=adata, obs_key='hierarchical_query')
````

`tl.label` (expanded)
---------------------

Location: `cytopus/tl/label.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/label.py>`_

Functions and responsibilities summarized above. Important details:

- ``label_marker_genes`` supports passing either a ``KnowledgeBase`` or a flat dict of gene sets. When a ``KnowledgeBase`` is passed, it collapses ``G.celltype_process_dict`` into a flat gene-set dictionary for labeling.
- ``get_gmt`` constructs a table-like GMT-friendly representation (pandas DataFrame) and can save to a file.

Example usage
~~~~~~~~~~~~~

````python
import cytopus as cp
from cytopus.tl.label import label_marker_genes
G = cp.KnowledgeBase()
G.get_celltype_processes(['B','T'], global_celltypes=['all-cells'])
marker_sets = [['CD19','MS4A1'], ['CD3D','CD3E']]
df = label_marker_genes(marker_sets, G, threshold=0.3)
print(df)
````

Example usage
~~~~~~~~~~~~~

````python
from cytopus.tl.hierarchy import Hierarchy, get_hierarchy_dict
G = cp.KnowledgeBase()
hier = Hierarchy(get_hierarchy_dict(G))
# Attach cells from an AnnData object (requires `anndata` package)
hier.add_cells(adata, obs_columns=['annotation'])
hier.query_ancestors('B', adata=adata, obs_key='hierarchical_query')
````

`tl.label` (expanded)
---------------------

Location: `cytopus/tl/label.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/label.py>`_

Functions and responsibilities summarized above. Important details:

- ``label_marker_genes`` supports passing either a ``KnowledgeBase`` or a flat dict of gene sets. When a ``KnowledgeBase`` is passed, it collapses ``G.celltype_process_dict`` into a flat gene-set dictionary for labeling.
- ``get_gmt`` constructs a table-like GMT-friendly representation (pandas DataFrame) and can save to a file.

Example usage
~~~~~~~~~~~~~

````python
import cytopus as cp
from cytopus.tl.label import label_marker_genes
G = cp.KnowledgeBase()
G.get_celltype_processes(['B','T'], global_celltypes=['all-cells'])
marker_sets = [['CD19','MS4A1'], ['CD3D','CD3E']]
df = label_marker_genes(marker_sets, G, threshold=0.3)
print(df)
````
