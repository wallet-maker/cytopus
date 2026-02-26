Installation Guide
===================

This page documents installation and setup steps based only on repository contents.

System requirements
-------------------

- Python (no explicit version specified in repository). The package uses `setuptools` and `networkx>2.7` in `setup.py`.

Dependencies (declared in repository)
--------------------------------------

- `networkx>2.7` (declared in `setup.py` `install_requires`).

Repository code references additional packages that are required for certain functionality (not listed in `setup.py` but used in code and tutorial notebooks):

- `numpy`
- `pandas`
- `matplotlib`
- `anndata`
- `pyvis` (optional; used in `kb_queries.plot_graph_interactive`)
- `pygraphviz` (optional; recommended for some plotting; README shows conda install suggestion)

Note: The repository explicitly mentions optional plotting libraries in `README.md` and the code attempts to import `pyvis`/`pygraphviz` where applicable.

Install from PyPI (as per README)
----------------------------------

The repository README states an installation command for users (this requires the package to be published on PyPI)::

    pip install cytopus

Install from source (as per README)
------------------------------------

The README documents installing from the GitHub source::

    pip install git+https://github.com/wallet-maker/cytopus.git

Install locally from this repository
------------------------------------

From the repository root run one of the following commands::

    # editable install (recommended for development)
    pip install -e .

    # normal install
    pip install .

Optional plotting / visualization
---------------------------------

- Install `pygraphviz` via conda (README suggestion)::

    conda install --channel conda-forge pygraphviz

- Install `pyvis` via pip for interactive HTML visualizations:

```
pip install pyvis
```

Verify installation
-------------------

A minimal verification in a Python shell (derived from README and code):

```
import cytopus as cp
G = cp.KnowledgeBase()
print(G.celltypes)   # should list cell type nodes
print(list(G.processes.keys())[:5])  # show some process gene set names
```

If those calls succeed and do not raise import errors, the core package and the default KnowledgeBase file are available.

Notes and missing information
-----------------------------

- The repository does not specify an exact Python version requirement.
- Not all runtime imports (e.g., `pandas`, `numpy`, `anndata`) are present in `setup.py` `install_requires`; install them as needed for full functionality.
