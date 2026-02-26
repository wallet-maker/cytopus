# Read the Docs Deployment Guide

This document provides instructions for deploying the Cytopus documentation to Read the Docs.

## Prerequisites

- Repository hosted on GitHub (or GitLab, Bitbucket)
- Read the Docs account at https://readthedocs.org/
- Admin access to the GitHub repository

## Deployment Steps

### 1. Connect Repository to Read the Docs

1. Log in to https://readthedocs.org/ with your GitHub account
2. Click "Import a Project"
3. Authenticate with GitHub if prompted
4. Search for and select the `wallet-maker/cytopus` repository
5. Click the "Create" button

### 2. Verify Configuration

Read the Docs will automatically detect the `.readthedocs.yaml` file. Verify:

- **Python version**: 3.10 (as specified in `.readthedocs.yaml`)
- **Documentation source**: `docs/source/conf.py`
- **Build format**: HTML (and optionally PDF/EPUB)

### 3. Configure Build Settings (Optional)

In the Read the Docs admin panel:

- Set default version: `latest`
- Enable PDF builds if desired
- Configure email notifications for build failures

### 4. Trigger Initial Build

- Push a commit to the repository or click "Build" in the Read the Docs admin
- Monitor the build progress in the "Builds" tab
- Once complete, the documentation will be available at `https://cytopus.readthedocs.io/`

## Automatic Builds

After initial setup, documentation will automatically rebuild whenever you:

- Push commits to the main branch
- Create pull requests (with preview builds)
- Create new releases or tags

## File Structure Overview

The documentation is fully configured with:

- **`.readthedocs.yaml`** — Read the Docs configuration
- **`docs/source/conf.py`** — Sphinx configuration
- **`docs/source/*.rst`** — Documentation pages
- **`docs/source/_static/custom.css`** — Custom styling
- **`docs/requirements-docs.txt`** — Build dependencies

## Troubleshooting Build Failures

If a build fails:

1. Check the build log in the Read the Docs admin panel
2. Verify all Python dependencies are in `docs/requirements-docs.txt`
3. Ensure all `.rst` files have proper reStructuredText syntax
4. Check that the Sphinx configuration in `conf.py` is correct

Common issues:

- **Missing dependencies**: Add to `docs/requirements-docs.txt`
- **Import errors**: Use `autodoc` directive or import stubs in `conf.py`
- **RST syntax errors**: Validate with `sphinx-build -b html source build` locally

## Building Locally

For testing before deployment:

```bash
cd docs
pip install -r requirements-docs.txt
sphinx-build -b html source build
# Open build/index.html in a browser
```

## Documentation Features

The generated documentation includes:

- **Responsive design** — Works on desktop, tablet, and mobile
- **Full-text search** — Searchable across all pages
- **Versioning** — Multiple documentation versions for different releases
- **PDF export** — Downloadable PDF version
- **Dark mode support** — User-selectable theme

## Documentation Structure

The documentation is organized into three main sections:

### Getting Started
- Project Overview
- Installation Guide
- Project Structure

### User Guide
- Architecture & Workflow
- API Reference
- Configuration & Data

### Development
- Scripts & Automation
- Contributing Guidelines
- Troubleshooting & FAQ

## Customization

To further customize the documentation:

1. **Theme options**: Edit `html_theme_options` in `docs/source/conf.py`
2. **Styling**: Edit `docs/source/_static/custom.css`
3. **Logo/Branding**: Add a logo image to `docs/source/_static/` and configure in `conf.py`
4. **Additional extensions**: Add Sphinx extensions to `extensions` list in `conf.py`

## Support

For Read the Docs help and documentation:
- Visit: https://docs.readthedocs.io/
- Community: https://readthedocs.org/support/

For Sphinx documentation:
- Visit: https://www.sphinx-doc.org/
