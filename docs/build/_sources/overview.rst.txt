Project Overview
================

**Project name:** Cytopus

**Purpose:**

Cytopus is a Python package providing a KnowledgeBase for single-cell genomics annotations. The repository contains a packaged graph-form KnowledgeBase, utilities to query and visualize that KnowledgeBase, and helper tools to construct a KnowledgeBase from raw annotations.

**Problem it solves:**

- Stores gene sets and a cell-type hierarchy as a NetworkX graph and exposes domain-specific queries to retrieve process gene sets and identity gene sets for cell types.
- Enables exporting gene-set dictionaries for downstream analysis (e.g., Spectra) and offers utilities to label marker genes and annotate factor analysis outputs.

**Core objectives:**

- Provide a programmatic KnowledgeBase object wrapping a NetworkX DiGraph.
- Offer tools to construct, query, export, and visualize the KnowledgeBase.
- Ship an example KnowledgeBase and example AnnData for tutorials.

**Technology stack (from repository):**

- Language: Python (package with `setup.py`)
- Primary libraries referenced: `networkx` (required), `matplotlib`, `numpy`, `pandas`, `anndata` (used by some utilities), `pyvis`/`pygraphviz` (optional, for visualization)
- Data: networkx pickled graph files under `cytopus/data/`
- Tutorials: Jupyter notebooks under `notebooks/`

**High-level architecture:**

- `cytopus.knowledge_base.KnowledgeBase` — core wrapper class for a NetworkX DiGraph containing nodes (genes, gene sets, cell types) and typed edges.
- `cytopus.tl` subpackage — tools for constructing KnowledgeBases (`create.py`), managing hierarchies (`hierarchy.py`), and labeling/export helpers (`label.py`).
- Packaged data in `cytopus/data/` serves as default input for `KnowledgeBase()` when no graph is supplied.

See the other pages for installation, file-level documentation, and API details.
