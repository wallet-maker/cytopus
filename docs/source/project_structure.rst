Project Structure
=================

This section lists the files and directories present in the repository and explains their purpose. All information is taken from the repository contents.

Top-level files and folders
---------------------------

- `README.md <https://github.com/wallet-maker/cytopus/blob/main/README.md>`_ — Project README and quickstart usage examples (used as reference in this documentation).
- `setup.py <https://github.com/wallet-maker/cytopus/blob/main/setup.py>`_ — Packaging and install configuration. Declares `networkx>2.7` in `install_requires` and includes `data/*.txt` and `data/*.h5ad` as package data.
- `LICENSE <https://github.com/wallet-maker/cytopus/blob/main/LICENSE>`_ — project license file (present in repository root).
- `notebooks/ <https://github.com/wallet-maker/cytopus/tree/main/notebooks>`_ — collection of Jupyter notebooks used as tutorials.
- `cytopus/ <https://github.com/wallet-maker/cytopus/tree/main/cytopus>`_ — main Python package directory.

`cytopus/` package contents (code and data)
-------------------------------------------

- `cytopus/__init__.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/__init__.py>`_ — top-level package exports. Imports `tl` package and exposes `KnowledgeBase` and `get_data` from `knowledge_base`.
- `cytopus/knowledge_base/ <https://github.com/wallet-maker/cytopus/tree/main/cytopus/knowledge_base>`_ — KnowledgeBase implementation and helper functions:
  * `__init__.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/__init__.py>`_ — re-exports `KnowledgeBase` and `get_data`.
  * `kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_ — main implementation of the `KnowledgeBase` class, `get_data`, and helper functions such as `extract_hierarchy`.
- `cytopus/tl/ <https://github.com/wallet-maker/cytopus/tree/main/cytopus/tl>`_ — utility tools and helpers for constructing, labeling, and handling hierarchies:
  * `__init__.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/__init__.py>`_ — imports `label`, `create`, `hierarchy`.
  * `create.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/create.py>`_ — `construct_kb(...)` helper to build a Graph and return `KnowledgeBase`.
  * `label.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/label.py>`_ — labeling utilities (overlap coefficient, `label_marker_genes`, exports to GMT, CSV helpers).
  * `hierarchy.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/hierarchy.py>`_ — hierarchy utilities and a `Hierarchy` class for hierarchical graphs and cell attachments.
- `cytopus/data/ <https://github.com/wallet-maker/cytopus/tree/main/cytopus/data>`_ — packaged data files included by `setup.py`:
  * `Cytopus_1.31nc.txt <https://github.com/wallet-maker/cytopus/blob/main/cytopus/data/Cytopus_1.31nc.txt>`_ (default KnowledgeBase file used when `KnowledgeBase()` is called without arguments)
  * `Cytopus_1.23.txt <https://github.com/wallet-maker/cytopus/blob/main/cytopus/data/Cytopus_1.23.txt>`_
  * `Cytopus_1.22.txt <https://github.com/wallet-maker/cytopus/blob/main/cytopus/data/Cytopus_1.22.txt>`_
  * `Cytopus_1.2.txt <https://github.com/wallet-maker/cytopus/blob/main/cytopus/data/Cytopus_1.2.txt>`_
  * `adata_spectra.h5ad <https://github.com/wallet-maker/cytopus/blob/main/cytopus/data/adata_spectra.h5ad>`_ (example AnnData file used in tutorials)

Jupyter notebooks (under `notebooks/`)
---------------------------------------

- `Cytopus_utils_tutorial.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/Cytopus_utils_tutorial.ipynb>`_
- `Hierarchical_annotation_tutorial.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/Hierarchical_annotation_tutorial.ipynb>`_
- `KnowledgeBase_construct.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/KnowledgeBase_construct.ipynb>`_
- `KnowledgeBase_queries_colaboratory.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/KnowledgeBase_queries_colaboratory.ipynb>`_
- `Utils_tutorial.ipynb <https://github.com/wallet-maker/cytopus/blob/main/notebooks/Utils_tutorial.ipynb>`_

Code file counts (from repository scan)
---------------------------------------

- Python files: 8 (including `setup.py` and package Python modules)
  * `setup.py <https://github.com/wallet-maker/cytopus/blob/main/setup.py>`_
  * `cytopus/__init__.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/__init__.py>`_
  * `cytopus/tl/__init__.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/__init__.py>`_
  * `cytopus/tl/create.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/create.py>`_
  * `cytopus/tl/hierarchy.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/hierarchy.py>`_
  * `cytopus/tl/label.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/label.py>`_
  * `cytopus/knowledge_base/__init__.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/__init__.py>`_
  * `cytopus/knowledge_base/kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_
- Jupyter notebooks: 5
- Data files in `cytopus/data/`: 5

Purpose of key files
--------------------

- `cytopus/knowledge_base/kb_queries.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/knowledge_base/kb_queries.py>`_ — Core class `KnowledgeBase` and helper functions to load, filter, and query the graph.
- `cytopus/tl/create.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/create.py>`_ — Build a graph from raw edges and metadata and return a `KnowledgeBase` instance.
- `cytopus/tl/hierarchy.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/hierarchy.py>`_ — Utilities to transform hierarchies to nested dicts, construct hierarchical graphs, and the `Hierarchy` class that can attach cell barcodes from AnnData.
- `cytopus/tl/label.py <https://github.com/wallet-maker/cytopus/blob/main/cytopus/tl/label.py>`_ — Functions to compute overlap coefficients, label marker gene sets against KnowledgeBase gene sets, and export gene sets/hierarchies to GMT/CSV.

Entry point
-----------

- There is no executable CLI or console script defined in the repository. The library entry point for usage is the Python importable package `cytopus` (i.e., `import cytopus as cp`). Typical usage then constructs a `KnowledgeBase()` object.
