# Configuration & Packaged Data

This page summarizes configuration-relevant files and packaged data (based only on repository contents).

## Packaged data
The package includes prebuilt KnowledgeBase files under `cytopus/data/`. The default used by `KnowledgeBase()` when no `graph` is provided is `Cytopus_1.31nc.txt`.

Files included in `cytopus/data/`:
- `Cytopus_1.31nc.txt` — default KnowledgeBase file referenced by `KnowledgeBase.__init__`.
- `Cytopus_1.23.txt`
- `Cytopus_1.22.txt`
- `Cytopus_1.2.txt`
- `adata_spectra.h5ad` — example AnnData object used by tutorials.

`setup.py` configures package data inclusion with:

```python
package_data={'cytopus': ['data/*.txt','data/*.h5ad']}
```

## Configuration files
- `setup.py` — packaging configuration and `install_requires` (declares `networkx>2.7`).
- `mkdocs.yml` (added in repository root by this documentation) — site configuration for MkDocs-based builds (not part of the original repository but provided to enable immediate documentation hosting).

## Environment variables or external configuration
- The repository contains no explicit environment variable requirements or configuration files (e.g., `.env`, settings files) defining runtime behavior.
- Any configuration of optional plotting libraries or Graphviz installation is expected to be handled at system level (e.g., installing `pygraphviz` via conda) as referenced in `README.md`.

## How configuration impacts execution
- Absence of optional packages (`pyvis`, `pygraphviz`, `matplotlib`, `pandas`, `numpy`, `anndata`) will limit functionality of plotting, CSV/GMT export, and AnnData integration; the core `KnowledgeBase` relies on `networkx` and pickled graph files.
- The `KnowledgeBase` constructor defaults to `Cytopus_1.31nc.txt` when no graph is supplied.

## Notes about missing/unspecified information
- The repository does not define configuration files for logging, environment, or runtime flags.
- No explicit Python version is specified in repository files.
