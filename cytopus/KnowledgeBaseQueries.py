import pandas as pd
import numpy as np
import networkx as nx
from pyvis.network import Network
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
import pickle

class KnowledgeBase:
    def __init__(self, graph_path):
        '''
        load KnowledgeBase from file
        retrieve all cell types in KnowledgeBase
        create dictionary for cellular processes in KnowledgeBase
        '''
        # load KnowledgeBase from pickled file
        with open(graph_path, 'rb') as f:  # notice the r instead of w
            self.graph = pickle.load(f) 
        #retrieve all cell types from data
        self.celltypes = self.filter_nodes(attribute_name = 'class',attributes= ['cell_type'],origin=None,target=None)
        
        #create gene set : gene dict for all cellular processes
        self.processes = self.get_processes(gene_sets = list(set([x[0] for x in self.filter_edges(attribute_name = 'class', attributes = ['process_OF'],target=self.celltypes)])))
        #self.processes = #self.filter_nodes(attribute_name = 'class',attributes= ['processes'],origin=None,target=None)
        print(self)
        #create gene set : gene dict for all cellular identities
        self.identities = self.get_identities(self.filter_nodes(attribute_name = 'class',attributes=['cell_type'],origin=None,target=None))
        
    
    def __str__(self):
        print(f"KnowledgeBase object containing {len(self.celltypes)} cell types and {len(self.processes)} cellular processes")
        return ""
    
    def filter_nodes(self, attributes,attribute_name=None, 
                     origin=None,target=None):
        '''
        filter nodes in networkx graph for their attributes
        G: networkx graph
        attribute_name: attribute key in node dictionary
        attributes: attribute list to select
        origin: list of node origin of node
        target: list of node end/target
        '''
        node_list = [] 
        if attribute_name == None:
            node_list = self.graph.nodes
        else:
            for x in self.graph.nodes:
                if attribute_name in self.graph.nodes[x].keys():
                    if self.graph.nodes[x][attribute_name] in attributes:
                        node_list.append(x)
        if origin!=None:
            node_list = [x for x in node_list if x[0] in origin]
        if target!=None:
            node_list = [x for x in node_list if x[1] in target]
        return node_list


    def filter_edges(self,  attributes,attribute_name=None,origin=None,target=None, ):
        '''
        filter edges in networkx graph for their attributes
        G: networkx graph
        attribute_name: attribute key in node dictionary
        attributes: attribute list to select
        origin: list of node origin of edge
        target: list of node end/target
        '''
        edge_list = [] 
        if attribute_name == None:
            edge_list = self.graph.edges
        else:
            for x in self.graph.edges:
                if attribute_name in self.graph.edges[x].keys():
                    if self.graph.edges[x][attribute_name] in attributes:
                        edge_list.append(x)
        if origin!=None:
            edge_list = [x for x in edge_list if x[0] in origin]
        if target!=None:
            edge_list = [x for x in edge_list if x[1] in target]
        return edge_list

    def get_processes(self,gene_sets): 
        '''
        create dictionary gene sets for cellular processes : genes
        self: KnowledgeBase object (networkx)
        gene_sets: list of gene sets for cellular processes
        '''
        #get genes per gene set    
        genes = self.filter_nodes(attribute_name = 'class',attributes =['gene'])
        #dictionary geneset : genes
        gene_edges = self.filter_edges( attribute_name ='class', attributes = ['gene_OF'],origin=gene_sets,target=genes)
        gene_set_dict = {}
        for i in gene_edges:
            if i[0] in gene_set_dict.keys():
                gene_set_dict[i[0]].append(i[1])#
            else:
                gene_set_dict[i[0]]= [i[1]]#
        return gene_set_dict
    
    def get_celltype_processes(self,celltypes,global_celltypes=None,get_parents =True,get_children =True):
        '''
        get gene sets for specific cell types
        self: KnowledgeBase object (networkx)
        celltypes: list of celltypes to retrieve
        global_celltypes: list of celltypes to set as 'global' for Spectra
        get_parent: also retrieve gene sets for the parents of the cell types in celltypes
        get_children: also retrieve gene sets for the parents of the cell types in celltypes
        '''
        import itertools

        ## limit to celltype subgraph to retrieve relevant celltypes

        node_list_plot = self.celltypes

        def filter_node(n1):
            return n1 in node_list_plot

        view = nx.subgraph_view(self.graph, filter_node=filter_node)

        for x in list(set(celltypes+global_celltypes)):
            if x not in list(view.nodes):
                raise ValueError('Not all cell types are contained in the Immune Knowledge base')

        all_celltypes_parents = {}

        for i in celltypes:
            if i in view.nodes:
                all_celltypes_parents[i]=  [n for n in nx.traversal.bfs_tree(view, i)]
            else:
                all_celltypes_parents[i]=  [i]
                print('cell type of interest',i,'is not in the input graph')

        all_celltypes_children = {}

        for i in celltypes:
            if i in view.nodes:
                all_celltypes_children[i]=  [n for n in nx.traversal.bfs_tree(view, i,reverse=True)]
            else:
                all_celltypes_children[i]=  [i]
                print('cell type of interest',i,'is not in the input graph')

        if get_parents ==True and  get_children==True:
            all_celltypes = list(itertools.chain.from_iterable(list(all_celltypes_children.values())+list(all_celltypes_parents.values())))
        elif get_parents==True and  get_children==False:
            all_celltypes = list(itertools.chain.from_iterable(list(all_celltypes_parents.values())))
        elif get_parents == False and  get_children==True:
            all_celltypes = list(itertools.chain.from_iterable(list(all_celltypes_children.values())))
        else:
            all_celltypes = []
        all_celltypes = list(set(all_celltypes +  global_celltypes + celltypes))
        #get process genesets connected to these celltypes
        gene_set_edges =self.filter_edges(attribute_name = 'class', attributes = ['process_OF'],target=all_celltypes)  
        
        #dictionary gene set  : cell type
        gene_set_celltype_dict = dict(gene_set_edges)

        #dictionary cell type: [gene_set1, gene_set2,...]
        celltype_gene_set_dict = {}

        for key,value in gene_set_celltype_dict.items():
            if value in celltype_gene_set_dict.keys():
                celltype_gene_set_dict[value].append(key)
            else:
                celltype_gene_set_dict[value] = [key]
        
        #get dict gene sets for cellular processes : genes
        gene_set_dict = self.get_processes(gene_sets = list(set([x[0] for x in gene_set_edges])))
        
        #construct dictionary
        process_dict = {}

        for key,value in celltype_gene_set_dict.items():
            process_dict[key] = {}
            for gene_set in value:
                process_dict[key][gene_set] = gene_set_dict[gene_set]
        if global_celltypes != None:
            global_gs = {}
            for i in global_celltypes:
                if i in process_dict.keys():
                    global_gs = global_gs | process_dict[i]
                    del process_dict[i]
                else:
                    print('did not find',i,'in cell type keys to set as global')
            process_dict['global'] = global_gs

        else:
            print('you have to add a "global" key to run Spectra. E.g. set one key as "global"')

        ## merge relevant children and parents into cell type specific keys

        process_dict_merged = {}

        if get_children:
            for key,value in all_celltypes_children.items():
                merged_dict = {}
                for cell_type in value:
                    if cell_type in process_dict.keys():
                        merged_dict = merged_dict | process_dict[cell_type]
                process_dict_merged[key]=merged_dict 

        if get_parents:
            for key,value in all_celltypes_parents.items():
                merged_dict = {}
                for cell_type in value:
                    if cell_type in process_dict.keys():
                        merged_dict = merged_dict | process_dict[cell_type]
                process_dict_merged[key]=merged_dict 
        process_dict_merged['global'] = process_dict['global']
            
        ## eck if cell types contain shared children or parents
        import itertools
        from collections import Counter
        if get_children:
            shared_children = []
            for key,value in Counter(list(itertools.chain.from_iterable(list(all_celltypes_children.values())))).items():
                if value >1:
                    shared_children.append(key)
            if shared_children != []:

                print('cell types of interest share the following children:',shared_children,'Generally, this is not desirable.')
        if get_parents:
            shared_parents = []
            for key,value in Counter(list(itertools.chain.from_iterable(list(all_celltypes_parents.values())))).items():
                if value >1:
                    shared_parents.append(key)
            if shared_parents != []:
                print('cell types of interest share the following parents:',shared_parents,'This may be desired.')

        self.celltype_process_dict = process_dict_merged
        self.processes = gene_set_dict
        
    def get_identities(self, celltypes):
        '''
        self: KnowledgeBase object (networkx)
        celltypes: list of cell types to retrieve identity gene sets for
        '''
        identity_edges = self.filter_edges( attribute_name ='class', attributes = ['identity_OF'],target=celltypes)
        
        #construct dictionary geneset:gene
        gene_edges = self.filter_edges( attribute_name ='class', attributes = ['gene_OF'])
        gene_set_dict = {}
        for i in gene_edges:
            if i[0] in gene_set_dict.keys():
                gene_set_dict[i[0]].append(i[1])#
            else:
                gene_set_dict[i[0]]= [i[1]]#
                      
        #construct dictionary celltype: identity_geneset
        identity_dict = {}
        
        for edge in identity_edges:
            if edge[1] in list(self.celltypes):
                identity_gs = gene_set_dict[edge[0]]
                identity_dict[edge[1]] = identity_gs
            else:
                print(edge[1],'not contained in KnowledgeBase')
        return identity_dict
        
    def plot_celltypes(self, figure_size = [30,30], node_size = 1000, edge_width= 1, arrow_size=20, 
                       edge_color= 'k', node_color='#8decf5', label_size = 20):
        node_list_plot = self.filter_nodes(attribute_name='class', attributes = ['cell_type'])

        def filter_node(n1):
            return n1 in node_list_plot

        plt.rcParams["figure.figsize"] = figure_size
        plt.rcParams["figure.autolayout"] = True

        view = nx.subgraph_view(self.graph, filter_node=filter_node)

        pos=graphviz_layout(view)

        nodes = nx.draw_networkx_nodes(view, pos=pos,node_color=node_color,nodelist=None,node_size=node_size,label=True)
        edges = nx.draw_networkx_edges(view, pos=pos, edgelist=None, width=edge_width, edge_color=edge_color, style='solid', alpha=None, arrowstyle=None, 
                                       arrowsize=arrow_size, 
                            edge_cmap=None, edge_vmin=None, edge_vmax=None, ax=None, arrows=None, label=None, 
                            node_size=node_size, nodelist=None, node_shape='o', connectionstyle='arc3', 
                            min_source_margin=0, min_target_margin=0)
        labels = nx.draw_networkx_labels(view,pos=pos,font_size=label_size)
        print('all celltypes in knowledge base:',list(labels.keys()))
         