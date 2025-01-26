#import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
def build_nested_dict(graph, node):
    '''
    build nested dictionary from reverse view of cytopus cell type hierarchy
    graph: networkx.DiGraph.view, reverse view of Cytopus cell type hierarchy
    root: str, name of root node in the reversed view
    '''
    nested_dict = {node: {}}
    for neighbor in graph.successors(node):
        nested_dict[node].update(build_nested_dict(graph, neighbor))
    return nested_dict

def get_hierarchy_dict(G):
    '''
    reverse Cytopus cell type hierarchy and build nested hierarchy from it
    G: Cytopus.KnowledgeBase, containing cell type hierarchy

    '''
    import networkx as nx
    #get view of cell type hierarchy
    node_list_plot = G.filter_nodes(attribute_name='class', attributes = ['cell_type']) 
    def filter_node(n1):
                return n1 in node_list_plot
        
    view = nx.subgraph_view(G.graph, filter_node=filter_node)

    #reverse graph view (going from least granular to most granular cell type)
    reversed_view = view.reverse(copy=True)
    root_nodes = [n for n in reversed_view.nodes if reversed_view.in_degree(n) == 0]

    #build the nested dictionary
    hierarchy_dict = {}
    for root in root_nodes:
        hierarchy_dict.update(build_nested_dict(reversed_view, root))

    return hierarchy_dict

def create_hierarchical_graph(data, type_label):
    import networkx as nx
    G = nx.DiGraph()
    for parent, children in data.items():
        G.add_node(parent)
        if isinstance(children, dict):
            child_node = create_hierarchical_graph(children,type_label)
            G.add_nodes_from(child_node.nodes(data=True))
            G.add_edges_from([(u,v) for u,v in child_node.edges()])
            for child in children:
                G.add_edge(child, parent)
        else:
            for child in children:
                G.add_node(child)
                G.add_edge(child, parent)
    nx.set_node_attributes(G, type_label,'type')
    return G

def get_all_keys(d):
    keys = set()
    for k, v in d.items():
        keys.add(k)
        if isinstance(v, dict):
            keys |= get_all_keys(v)
    return keys

def get_nodes_of_type(graph, node_type):
    nodes = [node for node in graph.nodes() if graph.nodes[node]['type'] == node_type]
    nodes.sort(key=lambda x: x.split('.'))
    return nodes

def get_indices(df, value):
    return df.index[df.astype(str).apply(lambda x: x == value).any(axis=1)].tolist()

def get_node_labels(graph, node_type):
    import networkx as nx
    nodes = [node for node in nx.dfs_postorder_nodes(graph) if graph.nodes[node]['type'] == node_type]
    return nodes[::-1]


class Hierarchy:
    import networkx as nx
    def __init__(self, hierarchy_dict):
        '''
        load hierarchy class
        hierarchy_dict: dict, nested dict containing the cell type hierarchy
        '''
        self.graph = create_hierarchical_graph(hierarchy_dict,type_label = 'cell_type')
        print(self.__str__())
        
    def __str__(self):
        all_celltypes = get_nodes_of_type(self.graph, 'cell_type')
        return f"Hierarchy class containing {len(all_celltypes)} cell types:{all_celltypes}"

    def identities(self):
        '''
        plot cell types contained in hierarchy
        '''
        get_nodes_of_type(self.graph, node_type='cell_type')
        cell_type
    def plot_celltypes(self, node_color='#8decf5', node_size = 1000,edge_width= 1,arrow_size=20 ,edge_color= 'k',label_size = 10, figsize=[30,30]):
        '''
        plot all cell types contained in hierarchy object
        '''
        

        #plt.rcParams["figure.figsize"] = figure_size
        #plt.rcParams["figure.autolayout"] = True
        import networkx as nx
        import matplotlib.pyplot as plt
        node_list_plot = get_nodes_of_type(self.graph, 'cell_type')
        def filter_node(n1):
                    return n1 in node_list_plot

        view = nx.subgraph_view(self.graph, filter_node=filter_node)

        pos=graphviz_layout(view)
        plt.rcParams["figure.figsize"] = figsize
        nodes = nx.draw_networkx_nodes(view, pos=pos,node_color=node_color,nodelist=None,node_size=node_size,label=True)
        edges = nx.draw_networkx_edges(view, pos=pos, edgelist=None, width=edge_width, edge_color=edge_color, style='solid', alpha=None, arrowstyle=None, 
                                        arrowsize=arrow_size, 
                            edge_cmap=None, edge_vmin=None, edge_vmax=None, ax=None, arrows=None, label=None, 
                            node_size=node_size, nodelist=None, node_shape='o', connectionstyle='arc3', 
                            min_source_margin=0, min_target_margin=0)
        labels = nx.draw_networkx_labels(view,pos=pos,font_size=label_size)

    def add_cells(self, adata, obs_columns=None):
        '''
        Add cells to their most granular annotation in the hierarchy object.
        adata: anndata.AnnData, containing the cell type annotations under adata.obs.
        obs_columns: list, list of columns in adata.obs where the cell type annotations are stored (recommended).
        '''
        import warnings
        import networkx as nx

        if obs_columns is None:
            adata_sub = adata.obs
        else:
            adata_sub = adata.obs[obs_columns]

        # Get cell type annotations from adata
        adata_celltypes = []
        if obs_columns is not None:
            for column in obs_columns:
                adata_celltypes += list(set(adata_sub[column]))
        adata_celltypes = set(adata_celltypes)
        
        # Retrieve cell type nodes from the hierarchy
        celltype_nodes = get_node_labels(self.graph, 'cell_type')
        missing_celltypes = adata_celltypes - set(celltype_nodes)

        # Warn if there are missing cell types
        if missing_celltypes:
            warnings.warn(
                f"Cell types {list(missing_celltypes)} are not contained in the hierarchy. Skipping..."
            )

        # Loop over cell types and assign cells
        used_barcodes = set()
        for cell_type in celltype_nodes:
            barcodes = get_indices(adata_sub, cell_type)
            for barcode in barcodes:
                # Check if the cell is already in the hierarchy
                if barcode in self.graph:
                    # Find the current cell type assignments
                    current_annotations = [
                        edge[0] for edge in self.graph.in_edges(barcode)
                        if self.graph.nodes[edge[0]]['type'] == 'cell_type'
                    ]

                    if current_annotations:
                        # Compare granularity with the new annotation
                        current_annotation = current_annotations[0]
                        if nx.has_path(self.graph, cell_type, current_annotation):
                            # Current annotation is less granular (upstream), skip
                            continue
                        elif nx.has_path(self.graph, current_annotation, cell_type):
                            # New annotation is less granular (upstream), update
                            self.graph.remove_edge(current_annotation, barcode)
                    
                    # If neither path exists, assume unrelated; continue adding the new annotation

                # Add the cell to the hierarchy
                self.graph.add_node(barcode, type='cell')
                self.graph.add_edge(cell_type, barcode)
            
    def query_ancestors(self, query_node, adata=None, obs_key='hierarchical_query'):
        '''
        retrieves all cell barcodes belonging to the cell type and all of its subsets
        query_node: str, cell type name fir which to retrieve barcodes 
        node_type: str, node type of cell type node (here: 'cell_type')
        adata: anndata.AnnData, adata to store the cell type annotations under adata.obs[obs_key]
        obs_key: str, column label to store cell tyoe annotations under adata.obs[obs_key]
        returns: dict, containing the barcodes belonging to each annotation in self.annotations, if adata is provided they will also be stored in adata.obs[obs_key]
        '''
        import networkx as nx
        import anndata
        node_type='cell_type'
        if node_type == self.graph.nodes[query_node]['type']:
            nodes_of_specific_type = [node for node in nx.ancestors(self.graph, query_node) if self.graph.nodes[node]['type'] == node_type]
            nodes_of_specific_type.append(query_node)
            cell_nodes = {}
            for node in set(nodes_of_specific_type):
                cell_edges = [edge for edge in self.graph.edges(node) if self.graph.nodes[edge[1]]['type'] == 'cell']
                cell_nodes[node] = [edge[1] for edge in cell_edges] 
            cell_nodes_inv = {}
            for k,v in cell_nodes.items():
                for i in v:
                    cell_nodes_inv[i] = k
            if isinstance(adata,anndata._core.anndata.AnnData):
                adata.obs[obs_key]= adata.obs_names.map(cell_nodes_inv)
            self.annotations =  cell_nodes
        else:
            print('query_node:',query_node,'should be of type',node_type,'stopping...')
    def get_cells_for_cell_type(self, cell_type):
        """
        Retrieve all cells assigned to a specific cell type in the hierarchy.
        cell_type: str, name of the cell type node to query.
        returns: ls, of cell barcodes assigned to the given cell type.
        """
        import networkx as nx

        # Check if the provided node is a valid cell type
        if cell_type not in self.graph.nodes:
            raise ValueError(f"Cell type '{cell_type}' does not exist in the hierarchy.")
        if self.graph.nodes[cell_type]['type'] != 'cell_type':
            raise ValueError(f"Node '{cell_type}' is not of type 'cell_type'.")

        # Retrieve all 'cell' nodes connected to the cell type node
        cell_edges = [
            edge for edge in self.graph.edges(cell_type)
            if self.graph.nodes[edge[1]]['type'] == 'cell'
        ]
        cells = [edge[1] for edge in cell_edges]

        return cells
        