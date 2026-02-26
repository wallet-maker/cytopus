# Project Structure

This section lists the files and directories present in the repository and explains their purpose. All information is taken from the repository contents.

## Top-level files and folders
- `README.md` — Project README and quickstart usage examples (used as reference in this documentation).
- `setup.py` — Packaging and install configuration. Declares `networkx>2.7` in `install_requires` and includes `data/*.txt` and `data/*.h5ad` as package data.
- `LICENSE` — project license file (present in repository root).
- `notebooks/` — collection of Jupyter notebooks used as tutorials.
- `cytopus/` — main Python package directory.

## `cytopus/` package contents (code and data)
- `cytopus/__init__.py` — top-level package exports. Imports `tl` package and exposes `KnowledgeBase` and `get_data` from `knowledge_base`.
- `cytopus/knowledge_base/` — KnowledgeBase implementation and helper functions:
  - `__init__.py` — re-exports `KnowledgeBase` and `get_data`.
  - `kb_queries.py` — main implementation of the `KnowledgeBase` class, `get_data`, and helper functions such as `extract_hierarchy`.
- `cytopus/tl/` — utility tools and helpers for constructing, labeling, and handling hierarchies:
  - `__init__.py` — imports `label`, `create`, `hierarchy`.
  - `create.py` — `construct_kb(...)` helper to build a Graph and return `KnowledgeBase`
  - `label.py` — labeling utilities (overlap coefficient, `label_marker_genes`, exports to GMT, CSV helpers)
  - `hierarchy.py` — hierarchy utilities and a `Hierarchy` class for hierarchical graphs and cell attachments
- `cytopus/data/` — packaged data files included by `setup.py`:
  - `Cytopus_1.31nc.txt` (default KnowledgeBase file used when `KnowledgeBase()` is called without arguments)
  - `Cytopus_1.23.txt`
  - `Cytopus_1.22.txt`
  - `Cytopus_1.2.txt`
  - `adata_spectra.h5ad` (example AnnData file used in tutorials)

## Jupyter notebooks (under `notebooks/`)
- `Cytopus_utils_tutorial.ipynb`
- `Hierarchical_annotation_tutorial.ipynb`
- `KnowledgeBase_construct.ipynb`
- `KnowledgeBase_queries_colaboratory.ipynb`
- `Utils_tutorial.ipynb`

## Code file counts (from repository scan)
- Python files: 8 (including `setup.py` and package Python modules)
  - `setup.py`
  - `cytopus/__init__.py`
  - `cytopus/tl/__init__.py`
  - `cytopus/tl/create.py`
  - `cytopus/tl/hierarchy.py`
  - `cytopus/tl/label.py`
  - `cytopus/knowledge_base/__init__.py`
  - `cytopus/knowledge_base/kb_queries.py`
- Jupyter notebooks: 5
- Data files in `cytopus/data/`: 5

## Purpose of key files
- `cytopus/knowledge_base/kb_queries.py` — Core class `KnowledgeBase` and helper functions to load, filter, and query the graph.
- `cytopus/tl/create.py` — Build a graph from raw edges and metadata and return a `KnowledgeBase` instance.
- `cytopus/tl/hierarchy.py` — Utilities to transform hierarchies to nested dicts, construct hierarchical graphs, and the `Hierarchy` class that can attach cell barcodes from AnnData.
- `cytopus/tl/label.py` — Functions to compute overlap coefficients, label marker gene sets against KnowledgeBase gene sets, and export gene sets/hierarchies to GMT/CSV.

## Entry point
- There is no executable CLI or console script defined in the repository. The library entry point for usage is the Python importable package `cytopus` (i.e., `import cytopus as cp`). Typical usage then constructs a `KnowledgeBase()` object.
