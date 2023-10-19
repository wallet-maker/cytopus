from os.path import dirname
import matplotlib.pyplot as plt
import pickle
import pkg_resources
import networkx as nx


def get_data(filename):
    """
    Load data from cytopus/data.
    """
    return pkg_resources.resource_filename('cytopus', 'data/' + filename)

#nested dict with celltype hierarchy
def extract_hierarchy(G, node='all-cells',invert=False):
    '''
    extract cell type hierarchy from KnowledgeBase object as a nested dictionary
    G: cytopus.kb.KnowledgeBase, with a hierarchy of celltypes
    node: str, celltype to use as starting points in the hiearchy (e.g. 'all-cells')
    invert: bool, if False the dict will contain all children below the node, if True the dict will contain all parents above the node
    '''
    node_list_plot = G.celltypes

    def filter_node(n1):
        return n1 in node_list_plot

    view = nx.subgraph_view(G.graph, filter_node=filter_node)
    if invert:
        predecessors = view.successors(node)
    else:
        predecessors = view.predecessors(node)
    if not predecessors:
        return node
    else:
        return {s: extract_hierarchy(G, s) for s in predecessors}


class KnowledgeBase:
    def __init__(self, graph=None):
        '''
        load KnowledgeBase from file
        retrieve all cell types in KnowledgeBase
        create dictionary for cellular processes in KnowledgeBase
        graph: str or networkx.DiGraph, path to pickled networkx.DiGraph object formatted for cytopus or networkx.DiGraph
        '''
        
        # Initialise default graph data
        if graph is None:
            graph = get_data("Cytopus_1.31nc.txt")
        # load KnowledgeBase from pickled file
        if isinstance(graph, nx.classes.digraph.DiGraph):
            self.graph = graph
        elif isinstance(graph, str):
            with open(graph, 'rb') as f:  # notice the r instead of w
                self.graph = pickle.load(f) 
        else:
            raise(ValueError('graph must be path (str) or networkx.classes.digraph.DiGraph object'))
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
    
    def get_celltype_hierarchy(self, node='all-cells',invert=False):
        #retrieve hierarchy
        hierarchy_dict = extract_hierarchy(self, node=node,invert=invert)
        return hierarchy_dict
        
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
    
    def get_celltype_processes(self,celltypes,global_celltypes=[None],get_parents =True, get_children =True, parent_depth=1, child_depth= None, fill_missing=True,parent_depth_dict=None, child_depth_dict=None, inplace=True):
        '''
        get gene sets for specific cell types
        self: KnowledgeBase object (networkx)
        celltypes: list of celltypes to retrieve
        global_celltypes: list of celltypes to set as 'global' for Spectra
        get_parent: also retrieve gene sets for the parents of the cell types in celltypes
        get_children: also retrieve gene sets for the parents of the cell types in celltypes
        fill_missing: add an empty dictionary for cell types not found in KnowledgeBase   
        parent_depth: steps from cell type to go up the hierarchie to retrieve gene sets linked to parents (e.g. 2 would be up to grandparents)
        parent_depth_dict: you can also set the depth for specific celltype with a dictionary {celltype1:depth1,celltype2:depth2}
        child_depth: steps from cell type to go down the hierarchie to retrieve gene sets linked to children (e.g. 2 would be down to grandchildren) 
        child_depth_dict: you can also set the depth for specific celltype with a dictionary {celltype1:depth1,celltype2:depth2}
        inplace: bool, if True save output under self.celltype_process_dict
        '''
        import itertools
        import warnings
        from collections import Counter

        ## limit to celltype subgraph to retrieve relevant celltypes

        node_list_plot = self.celltypes

        def filter_node(n1):
            return n1 in node_list_plot

        view = nx.subgraph_view(self.graph, filter_node=filter_node)

        for x in list(set(celltypes+global_celltypes)):
            if x not in list(view.nodes):
                warnings.warn('Not all cell types are contained in the Immune Knowledge base')
        if get_parents:
            all_celltypes_parents = {}
            if parent_depth_dict == None:
                parent_depth_dict = {}

            for i in celltypes:
                if i in view.nodes:#is celltype in KnowledgeBase
                    if i in parent_depth_dict.keys(): #check if query depth was manually defined
                        if parent_depth_dict[i] == None:
                            all_celltypes_parents[i]=  [i] 
                        else:
                            all_celltypes_parents[i]=  [n for n in nx.traversal.bfs_tree(view, i,depth_limit=parent_depth_dict[i])] 
                    else:
                        all_celltypes_parents[i]=  [n for n in nx.traversal.bfs_tree(view, i,depth_limit=parent_depth)] 
                elif fill_missing: 
                    all_celltypes_parents[i] = {} #if not add an empty dictionary
                    print('adding empty dictionary for cell type:',i)
                else:
                    all_celltypes_parents[i]=  [i]
                    print('cell type of interest',i,'is not in the knowledge base')
        if get_children:
            all_celltypes_children = {}
            if child_depth_dict == None:
                child_depth_dict = {}

            for i in celltypes:
                if i in view.nodes:
                    if i in child_depth_dict.keys(): #check if query depth was manually defined
                        if child_depth_dict[i] ==None:
                            all_celltypes_children[i]=  [i]
                        else:
                            all_celltypes_children[i]=  [n for n in nx.traversal.bfs_tree(view, i,reverse=True,depth_limit=child_depth_dict[i])]
                    else:
                        all_celltypes_children[i]=  [n for n in nx.traversal.bfs_tree(view, i,reverse=True,depth_limit=child_depth)]
                else:
                    all_celltypes_children[i]=  [i]
                    print('cell type of interest',i,'is not in the knowledge base')

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
        if global_celltypes != [None]:
            global_gs = {}
            for i in global_celltypes:
                if i in process_dict.keys():
                    global_gs.update(process_dict[i])
                    del process_dict[i]
                else:
                    print('did not find',i,'in cell type keys to set as global')
            process_dict['global'] = global_gs

        else:
            print('you must add a "global" key to run Spectra. E.g. set <global_celltypes> to one cell type key to be set as "global"')

        ## merge relevant children and parents into cell type specific keys

        process_dict_merged = {}
       
        
        if get_children:
            for key,value in all_celltypes_children.items():
                merged_dict = {}
                for cell_type in value:
                    if cell_type in process_dict.keys():
                        merged_dict.update(process_dict[cell_type])
                if key in process_dict_merged.keys():
                    process_dict_merged[key].update(merged_dict) 
                else:
                    process_dict_merged[key]=merged_dict 

        if get_parents:
            for key,value in all_celltypes_parents.items():
                merged_dict = {}
                for cell_type in value:
                    if cell_type in process_dict.keys():
                        merged_dict.update(process_dict[cell_type])
                if key in process_dict_merged.keys():
                    process_dict_merged[key].update(merged_dict)
                else:
                    process_dict_merged[key]=merged_dict 
                
        if get_children==False and get_parents==False:
            process_dict_merged =process_dict 
                
        if global_celltypes != [None]:
            process_dict_merged['global'] = process_dict['global']
            
        ## check if cell types contain shared children or parents
        if get_children:
            shared_children = []
            for key,value in Counter(list(itertools.chain.from_iterable(list(all_celltypes_children.values())))).items():
                if value >1:
                    shared_children.append(key)
            if shared_children != []:

                print('cell types of interest share the following children:',shared_children,'This may be desired.')
        if get_parents:
            shared_parents = []
            for key,value in Counter(list(itertools.chain.from_iterable(list(all_celltypes_parents.values())))).items():
                if value >1:
                    shared_parents.append(key)
            if shared_parents != []:
                print('cell types of interest share the following parents:',shared_parents,'This may be desired.')
        if inplace:
            self.celltype_process_dict = process_dict_merged
            #self.processes = gene_set_dict
        else:
            return process_dict_merged
        
    def get_identities(self, celltypes_identities,include_subsets=False):
        '''
        self: KnowledgeBase object (networkx)
        celltypes: list of cell types to retrieve identity gene sets for
        '''
        if include_subsets:
            def filter_node(n1):
                return n1 in self.celltypes
           

            celltype_view =nx.subgraph_view(self.graph, filter_node=filter_node)
            celltypes_new = []
            for i in celltypes_identities:
                nodes_of_specific_type =  [n for n in nx.traversal.bfs_tree(celltype_view, i,reverse=True)]
                celltypes_new += nodes_of_specific_type
            celltypes_identities = list(set(celltypes_new))
            
        identity_edges = self.filter_edges( attribute_name ='class', attributes = ['identity_OF'],target=celltypes_identities)
        
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
        ''''
        plot all celltypes contained in the KnowledgeBase using matplotlib and graphviz
        self: KnowledgeBase object (networkx)
        figure_size: figure size
        node_size: node size in graph
        edge_with: edge width in graph
        arrow_size: arrow size of directed edges
        edge_color: edge color
        node_color: node color
        label_size: size of node labels
        '''
        try:
            from networkx.drawing.nx_agraph import graphviz_layout
        except ModuleNotFoundError:
            print('please install graphviz')
            pass
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
        
    def plot_graph_interactive(self, attributes=['cell_type','cellular_process'],colors= ['red','blue'], save_path = 'graph.html'):
        '''
        plot excerpt from the KnowledgeBase using the pyvis package
        self: KnowledgeBase object (networkx)
        attributes: list of node classes to plot
        colors: list of colors in the order of the node classes to plot
        save_path: save path for .html file
        '''
        
        try:
            from pyvis.network import Network
        except ModuleNotFoundError:
            print('please install pyvis')
            pass          
        
        while len(attributes)!=len(colors):
            print('attributes and colors have to be same length')
            break
        while len(set(attributes))!=len(attributes):
            print('attributes have to be unique')
            break
            
        net = Network(notebook=True)
        #cell types
        node_list_plot = self.filter_nodes(attribute_name='class', attributes = attributes)
        def filter_node(n1):
                return n1 in node_list_plot
        view = nx.subgraph_view(self.graph,filter_node=filter_node)
        net.from_nx(view)

        for i in net.nodes:
            index_range = list(range(len(attributes)))
            for v in index_range:
                             if i['class'] == attributes[v]:
                               i['color']= colors[v]
               
        #show
        net.show(save_path)
         
