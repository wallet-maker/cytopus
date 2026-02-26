# Scripts & Automation

This repository contains code modules and Jupyter notebooks but does not include automated test scripts, CI/CD configuration, or formal build scripts. This page summarizes relevant items.

## Packaging / Install
- `setup.py` is the packaging entry and defines `install_requires` and `package_data`. Use `pip install .` or `pip install -e .` from the repository root to install the package locally.

## Tutorials / Notebooks
Jupyter notebooks are included under `notebooks/` and serve as manual tutorials and interactive examples:
- `Cytopus_utils_tutorial.ipynb`
- `Hierarchical_annotation_tutorial.ipynb`
- `KnowledgeBase_construct.ipynb`
- `KnowledgeBase_queries_colaboratory.ipynb`
- `Utils_tutorial.ipynb`

## Automation and testing
- No `tests/` directory or testing harness exists in the repository.
- No CI/CD configuration files (e.g., GitHub Actions, Travis CI, GitLab CI) are present.

## Deployment readiness
- To publish documentation to Read the Docs or MkDocs, the provided `mkdocs.yml` (created for this documentation) and `docs/` directory are ready for build. The repository itself contains no automated docs build configuration.

## Example manual verification steps (derived from README)
1. Install the package locally:

```bash
pip install -e .
```

2. Run a quick Python verification (in a Python REPL):

```python
import cytopus as cp
G = cp.KnowledgeBase()
print(len(G.celltypes))
print(list(G.processes.keys())[:5])
```

3. Optional: open notebooks for guided examples and inspecting `adata_spectra.h5ad`.

## Notes about missing automation
- The repository does not include scripts to auto-generate documentation, run tests, or perform linting — these would need to be added separately if desired.
