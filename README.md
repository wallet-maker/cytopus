# Cytopus :octopus: [![DOI](https://zenodo.org/badge/389175717.svg)](https://zenodo.org/badge/latestdoi/389175717)


## Single cell omics biology annotations

![Image of Cytopus](https://github.com/wallet-maker/cytopus/blob/main/img/cytopus_v1.1_stable_graph.png)


## Overview:

Package to query our single cell genomics KnowledgeBase.

If you use Cytopus :octopus: or its gene sets please cite the original Cytopus publication [![DOI](https://zenodo.org/badge/389175717.svg)](https://zenodo.org/badge/latestdoi/389175717). 

For details see the [license](https://github.com/wallet-maker/cytopus/blob/Cytopus_1.3/LICENSE)

The KnowledgeBase is provided in graph format based on the networkx package. Central to the KnowledgeBase is a cell type hierarchy and **cellular_processess** which correspond to the cell types in this hierarchy. Cell types are supported by gene sets indicative of their **cellular identities**. Moreover, the KnowledgeBase contains metadata about the gene sets such as author ship, the gene set topic etc.. 

The KnowledgeBase can be queried to retrieve gene sets for specific cell types and organize them in a dictionary format for downstream use with the [Spectra](https://github.com/dpeerlab/spectra) package: 

Please cite the original Cytopus publication [![DOI](https://zenodo.org/badge/389175717.svg)](https://zenodo.org/badge/latestdoi/389175717). 


## Installation

install from pypi:

```
pip install cytopus
```

install from source:

```
pip install git+https://github.com/wallet-maker/cytopus.git
```

Some plotting functions require pygraphviz or pyvis. Install either or both:

pygraphviz using conda:
```
conda install --channel conda-forge pygraphviz
```

pyvis using pip
```
pip install pyvis
```

## Tutorial

### Quickstart - Querying the Knowledge Base:

Retrieve default KnowledgeBase (human only):

```
import cytopus as cp
G = cp.kb.KnowledgeBase()
```
Retrieve custom KnowledgeBase (documentation to build KnowledgeBase object [here](https://github.com/wallet-maker/cytopus/blob/Cytopus_1.3/notebooks/KnowledgeBase_construct.ipynb)):
```
file_path = '~/dir1/dir2/knowledgebase_file.txt'
G = cp.kb.KnowledgeBase(file_path)
```
Access data in KnowledgeBase:
```
#list of all cell types in KnowledgeBase
G.celltypes
#dictionary of all cellular processes in KnowledgeBase as a dictionary {'process_1':['gene_a','gene_e','gene_y',...],'process_2':['gene_b','gene_u',...],...}
G.processes
#dictionary of all cellular identities in KnowledgeBase as a dictionary {'identity_1':['gene_j','gene_k','gene_z',...],'identity_2':['gene_y','gene_p',...],...}
G.identities
#dictionary with gene set properties (for cellular processes or identities)
G.graph.nodes['gene_set_name']
```

Plot the cell type hierarchy stored in the KnowledgeBase as a directed graph with edges pointing into the direction of the parents:
```
G.plot_celltypes()
```


![Image of Cell type hierarchy](https://github.com/wallet-maker/cytopus/blob/main/img/celltype_hierarchy_1.2.png)



Prepare a nested dictionary assigning cell types to their cellular processes and cellular processes to their corresponding genes. This dictionary can be used as an input for Spectra.

First, select the cell types which you want to retrieve gene sets for. 
These cell types can be selected from the cell type hierarchy (see .plot_celltypes() method above)
```
celltype_of_interest = ['M','T','B','epi']
```

Second, select the cell types which you want merge gene sets and set them as global gene sets for the Spectra package. These gene sets should be valid for all cell types in the data. 
```
##e.g. if you are working with different human cells
global_celltypes = ['all-cells']
##e.g. if you are working with human leukocytes
global_celltypes = ['all-cells','leukocyte']
##e.g. if you are working with B cells
global_celltypes = ['all-cells','leukocyte','B']
```

Third retrieve dictionary of format {celltype_a:{process_a:[gene_a,gene_b,...],...},...}.
Decide whether you want to merge gene sets for all children or all parents (unusual) of the selected cell types.
```
G.get_celltype_processes(celltype_of_interest,global_celltypes = global_celltypes,get_children=True,get_parents =False)
```

Fourth, dictionary will be stored in the KnowledgeBase
```
G.celltype_process_dict
```

### Detailed tutorial for Querying the Knowledge Base:
Learn how to explore the Knowledge Base and retrieve a dicitionary which can be used for [Spectra](https://github.com/dpeerlab/spectra):
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/wallet-maker/cytopus/blob/main/notebooks/KnowledgeBase_queries_colaboratory.ipynb)

### Detailed tutorial for Generating a cytopus Knowledge Base object:
Learn how to create a Knowledge Base object from gene sets annotations and cell type hierarchies stored in .csv files:
[here](https://github.com/wallet-maker/cytopus/blob/Cytopus_1.3/notebooks/KnowledgeBase_construct.ipynb)

### Utils tutorial - Labeling Factor Analysis Outputs (Spectra):
Learn how to label marker genes from factor analysis, determine factor cell type specificity and export the Knowledge Base content as .gmt files for other applications:
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/wallet-maker/cytopus/blob/main/notebooks/Cytopus_utils_tutorial.ipynb)


## you can submit gene sets to be added to the KnowledgeBase here:

https://docs.google.com/forms/d/e/1FAIpQLSfWU7oTZH8jI7T8vFK0Nqq2rfz6_83aJIVamH5cogZQMlciFQ/viewform?usp=sf_link

All submissions will be reviewed and if needed revised before they will be added to the database. This will ensure consistency of the annotations and avoid gene set duplication. Authorship will be acknowledged in the KnowledgeBase for all submitted gene sets which pass review and are added to the KnowledgeBase. You can also create entirely new KnowledgeBase objects with this package.

## Citation and Usage 

For gene sets from external sources you must also abide to the licenses of the original gene sets. To make this easier we have stored these in the Knowledge Base object:

```
import cytopus as cp
G = cp.kb.KnowledgeBase()
gene_set_of_interest = 'all_macroautophagy_regulation_positive'
print(G.graph.nodes[gene_set_of_interest])
```

