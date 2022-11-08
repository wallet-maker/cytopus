# Cytopus :octopus:

## A language for single cell omics biology

![Image of Cytopus](https://github.com/wallet-maker/cytopus/blob/main/cytopus_v1.1_stable_graph.png)

currently version 1.2 (stable)

## Overview:

Package to query our single cell genomics KnowledgeBase.

KnowledgeBase is provided in graph format based on the networkx package. Central to the KnowledgeBase is a cell type hierarchy and **cellular processess** which correspond to these cell types. Cell types are supported by gene sets indicative of their **cellular identities**. 

The KnowledgeBase can be queried to retrieve gene sets for specific cell types and organize them in a dictionary format for downstream use with the **Soectra** package https://github.com/dpeerlab/spectra: 

{cell_type:{gene_set:[gene_A,gene_B,...]}}

A tutorial can be found under https://github.com/wallet-maker/cytopus/notebooks/KnowledgeBase_queries.ipynb

## you can submit gene sets to be added to the KnowledgeBase here:

https://docs.google.com/forms/d/e/1FAIpQLSfWU7oTZH8jI7T8vFK0Nqq2rfz6_83aJIVamH5cogZQMlciFQ/viewform?usp=sf_link

All submissions will be reviewed by 2 (computational) biologists and if needed revised before they will be added to the database. This will ensure consistency of the annotations and avoid gene set duplication. Authorship will be acknoledged in the KnowledgeBase for all submitted gene sets added. to the KnowledgeBase
