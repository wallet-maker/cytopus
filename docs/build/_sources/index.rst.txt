Cytopus - Single Cell Omics KnowledgeBase
=========================================

A Python package providing a KnowledgeBase for single-cell genomics annotations.

**Cytopus** stores gene sets and a cell-type hierarchy as a NetworkX graph and exposes domain-specific queries to retrieve process gene sets and identity gene sets for cell types. The package enables exporting gene-set dictionaries for downstream analysis (e.g., Spectra) and offers utilities to label marker genes and annotate factor analysis outputs.

--------

**Project Status**

Cytopus is a Python package for querying and managing a KnowledgeBase of single-cell genomics annotations.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   overview
   installation
   project_structure

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   architecture
   api
   configuration
   notebook_codes

.. toctree::
   :maxdepth: 1
   :caption: Development

   scripts
   contributing
   troubleshooting

--------

Key Features
============

- **Graph-based KnowledgeBase**: Store cell types, gene sets, and cellular processes as a NetworkX DiGraph with semantic relationships.
- **Flexible Querying**: Retrieve gene sets for specific cell types with parent/child traversal using BFS.
- **Hierarchical Annotation**: Attach and query single cells within a hierarchical cell-type structure using AnnData objects.
- **Export Utilities**: Export to GMT, CSV, and Spectra-compatible formats.
- **Multiple KnowledgeBase Versions**: Packaged with several versions of the KnowledgeBase (1.2, 1.22, 1.23, 1.31).
- **Interactive Visualization**: Plot cell-type hierarchies and KnowledgeBase graph structures (with optional pygraphviz/pyvis).

Quick Start
===========

Install the package::

    pip install cytopus

Load the default KnowledgeBase::

    import cytopus as cp
    G = cp.KnowledgeBase()
    
    # List all cell types
    print(G.celltypes)
    
    # List all cellular processes
    print(list(G.processes.keys())[:10])

Query for cell-type-specific gene sets::

    # Get processes for B, T, and M (myeloid) cells
    G.get_celltype_processes(['B', 'T', 'M'], global_celltypes=['all-cells'])
    
    # Access the result
    print(G.celltype_process_dict)

Getting Started (README excerpt)
================================

The following content is taken from the repository `README.md` and presented here to help users get started quickly.

Overview
--------

Package to query our single cell genomics KnowledgeBase.

The KnowledgeBase is provided in graph format based on the networkx package. Central to the KnowledgeBase is a cell type hierarchy and **cellular_processess** which correspond to the cell types in this hierarchy. Cell types are supported by gene sets indicative of their **cellular identities**. Moreover, the KnowledgeBase contains metadata about the gene sets such as author ship, the gene set topic etc..

Basic usage example::

    import cytopus as cp
    # retrieve default KnowledgeBase (human only):
    G = cp.KnowledgeBase()
    # list of all cell types in KnowledgeBase
    print(G.celltypes)
    # dictionary of all cellular processes in KnowledgeBase
    print(G.processes)

Notes
-----

Some plotting functions require `pygraphviz` or `pyvis`. Install either or both as needed:

::

    # pygraphviz using conda
    conda install --channel conda-forge pygraphviz

    # pyvis using pip
    pip install pyvis

Learn More
==========

- Read the :doc:`Project Overview <overview>` for a high-level introduction.
- Check the :doc:`Installation Guide <installation>` for setup instructions.
- Explore the :doc:`API Reference <api>` for a complete function reference.
- Review the :doc:`Architecture & Workflow <architecture>` to understand how the system works.
- See :doc:`Troubleshooting <troubleshooting>` for solutions to common issues.
