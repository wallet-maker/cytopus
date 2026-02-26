# Cytopus Documentation Build

This directory contains the source files for the Cytopus documentation, built with Sphinx and the Read the Docs theme.

## Building locally

### Using Sphinx

From the repository root:

```bash
# Install dependencies
pip install sphinx sphinx_rtd_theme

# Build HTML documentation
cd docs
sphinx-build -b html source build

# Output will be in docs/build/
```

### Viewing the documentation

Open `docs/build/index.html` in a web browser to view the generated documentation.

## File structure

- `source/` — Sphinx reStructuredText source files
  - `conf.py` — Sphinx configuration
  - `index.rst` — Landing page
  - `*.rst` — Individual documentation pages
- `build/` — Generated HTML output (created after build)
- `requirements-docs.txt` — Python dependencies for building docs

## Deployment to Read the Docs

The `.readthedocs.yaml` file at the repository root contains the configuration for Read the Docs. The repository is ready to be built on Read the Docs by:

1. Pushing the repository to GitHub
2. Connecting the repository to Read the Docs
3. Triggering a build

The documentation will automatically build on every push to the main branch.

## Documentation structure

The documentation is organized into the following sections:

- **Getting Started**: Overview, Installation, Project Structure
- **User Guide**: Architecture, API Reference, Configuration
- **Development**: Scripts & Automation, Contributing, Troubleshooting

All documentation is generated directly from the repository contents, with no external assumptions.
