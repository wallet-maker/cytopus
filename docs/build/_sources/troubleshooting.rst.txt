Troubleshooting & FAQ
======================

This page lists potential issues and debugging tips derived strictly from repository code and README information.

Common issues (based on repository content)
-------------------------------------------

- ImportError for missing packages:
  * Cause: The code references packages not included in `setup.py` `install_requires` (e.g., `pandas`, `numpy`, `matplotlib`, `anndata`, `pyvis`, `pygraphviz`).
  * Fix: Install required packages manually, e.g., `pip install numpy pandas matplotlib anndata` and optionally `pip install pyvis`. For `pygraphviz`, the README suggests installing via conda: `conda install --channel conda-forge pygraphviz`.

- Graphviz/layout errors when plotting:
  * Cause: `networkx.drawing.nx_agraph.graphviz_layout` requires Graphviz and possibly `pygraphviz` or `pydot` installed.
  * Fix: Install `pygraphviz` (conda-forge recommended) or ensure Graphviz is installed system-wide.

- `KnowledgeBase` fails to load when a custom `graph` path is given:
  * Cause: `KnowledgeBase.__init__` expects a pickled NetworkX DiGraph when `graph` is a string and attempts `pickle.load`.
  * Fix: Ensure the provided `graph` file is a pickled networkx.DiGraph object. Alternatively, pass a live `networkx.DiGraph` instance directly.

- Empty results from `get_celltype_processes`:
  * Cause: Requested cell types may not exist in `self.celltypes` or traversal depth settings exclude relevant nodes.
  * Fix: Verify that the requested cell types exist in `G.celltypes` and adjust `parent_depth`, `child_depth`, or `include_subsets` flags as appropriate.

Debugging tips
--------------

- Inspect graph node/edge attributes to understand structure::

    G = cp.KnowledgeBase()
    list(G.graph.nodes(data=True))[:10]
    list(G.graph.edges(data=True))[:10]

- Use `G.filter_nodes(attribute_name='class', attributes=['cell_type'])` to list cell types recognized by the KnowledgeBase.
- Use `G.filter_edges(attribute_name='class', attributes=['process_OF'])` to list process gene-set to cell-type edges.

Not defined in repository
-------------------------

- There are no explicit FAQ items, error codes, or detailed logging instructions provided in the repository. If problems persist, inspect the source files under `cytopus/` and the notebooks in `notebooks/` for examples and further clues.
