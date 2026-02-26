Architecture & Internal Workflow
=================================

This section explains how components interact and how data flows through the system (derived strictly from repository code).

High-level components
---------------------

- `KnowledgeBase <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_ (in `cytopus/knowledge_base/kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_): central class wrapping a NetworkX DiGraph and exposing query and plotting methods.
- `tl <https://github.com/wallet-maker/cytopus/tree/main/cytopus/tl>`_ tooling (in `cytopus/tl/ <https://github.com/wallet-maker/cytopus/tree/main/cytopus/tl>`_): helper modules to construct KBs, manipulate hierarchies, and label/export gene sets.
- `data/ <https://github.com/wallet-maker/cytopus/tree/main/cytopus/data>`_ files: packaged pickled NetworkX graphs and example `AnnData` dataset used as defaults.

Graph model and semantics
-------------------------

- The KnowledgeBase stores elements in a NetworkX `DiGraph` with typed nodes and typed edges.
- Node attribute `class` is used to distinguish types such as `'gene'` and `'cell_type'`.
- Edge attribute `class` encodes relationships:
  - `'gene_OF'` — gene_set -> gene edges
  - `'process_OF'` — gene_set -> celltype edges for processes
  - `'identity_OF'` — gene_set -> celltype edges for identities
  - `'SUBSET_OF'` — child -> parent cell type edges (hierarchy)

These attributes are used by filtering utilities (`filter_nodes`, `filter_edges`) to select nodes/edges for queries.

Typical execution flow (example)
---------------------------------

1. User imports the package and constructs a `KnowledgeBase`::

    import cytopus as cp
    G = cp.KnowledgeBase()

- If `graph` argument is omitted, `KnowledgeBase.__init__` uses `get_data("Cytopus_1.31nc.txt")` to obtain the path to the packaged KB and unpickles the NetworkX DiGraph.
- If a NetworkX DiGraph is passed, it is used directly; if a path string is passed, it is unpickled and used.

2. On initialization, `KnowledgeBase` computes helper properties:
    * `self.celltypes` via `filter_nodes(attribute_name='class', attributes=['cell_type'])`.
    * `self.processes` via `get_processes(...)` (maps process gene sets to gene lists) and cleaned of NaN values.
    * `self.identities` via `get_identities(...)`.

3. Querying for cell-type-specific gene sets:

- `get_celltype_processes(celltypes, global_celltypes=[None], get_parents=True, get_children=True, ...)`:
  * Builds a `subgraph_view` limited to nodes identified as `cell_type`.
  * Traverses parents (predecessors) and children (successors) using BFS to requested depths.
  * Collects `process_OF` edges linking gene sets to the discovered cell types.
  * Uses `get_processes(...)` to map each gene set to its genes (via `gene_OF` edges).
  * Aggregates results into a nested dict mapping `celltype -> {gene_set_name: [genes]}` and optionally merges global cell types into a `global` key.
  * Stores result under `self.celltype_process_dict` if `inplace=True`.

4. Constructing a KB from raw edges:

- `tl.create.construct_kb(celltype_edges, geneset_gene_edges, geneset_celltype_edges, annotation_dict, metadata_dict=None, save=False, save_path=None)`
  * Validates and classifies gene sets into `processes` vs `identities` using `annotation_dict` values (`'cellular_process'` or `'cellular_identity'`).
  * Adds node and edge attributes for class semantics.
  * Builds a `networkx.DiGraph`, sets node metadata if provided, optionally pickles the graph, and returns `KnowledgeBase(graph=G)`.

5. Hierarchy and annotation of single cells:

- `tl.hierarchy.Hierarchy` can be built from a nested dictionary representation of the cell-type hierarchy.
- `Hierarchy.add_cells(adata, obs_columns=None)` attaches cell barcodes as nodes with type `'cell'` and edges to their most granular annotation by reading `adata.obs` fields.
- `Hierarchy.query_ancestors(query_node, adata=None, obs_key='hierarchical_query')` retrieves barcodes for a cell type and its subsets and optionally writes assigned labels into `adata.obs`.

Module responsibilities
-----------------------

- `cytopus/knowledge_base/kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_ — loading KB data, filtering utilities, constructing process/identity dictionaries, plotting (matplotlib/pygraphviz/pyvis wrappers), and main `KnowledgeBase` API.
- `cytopus/tl/create.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/create.py>`_ — graph construction from lists/metadata; returns `KnowledgeBase` instances.
- `cytopus/tl/hierarchy.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/hierarchy.py>`_ — nested hierarchy <-> graph conversion utilities and the `Hierarchy` class for cell attachments and queries.
- `cytopus/tl/label.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/label.py>`_ — overlap computations and labeling utilities focused on mapping marker gene lists or factor loadings to KB gene sets; also CSV/GMT export helpers.

Data flow summary
-----------------

- Source of truth: NetworkX DiGraph stored in `self.graph` of `KnowledgeBase`.
- Derived views/dictionaries: `self.celltypes`, `self.processes`, `self.identities`, `self.celltype_process_dict` (after calling `get_celltype_processes`).
- External inputs: pickled graph files in `cytopus/data/`, or user-provided graph or raw edges passed into `construct_kb`.

Design notes
------------

- The code uses node/edge attributes to encode domain semantics rather than separate classes; filtering functions use these attributes to select relevant nodes/edges.
- Plotting functions are optional and attempt to import optional dependencies, printing user-facing messages if missing.
- There is no CLI or automated test harness included in the repository.
