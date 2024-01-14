import pytest
import networkx as nx
import cytopus as cp

def test_loading_data():
    data = cp.get_data("Cytopus_1.31nc.txt")

def test_initialise_cytopus():
    """Test that graph is loaded correctly"""
    data = cp.get_data("Cytopus_1.31nc.txt")
    G = cp.KnowledgeBase(graph=data)
    assert G.graph.number_of_nodes() == 6113
    assert G.graph.number_of_edges() == 10785
    
    # test variables
    assert G.celltypes, f"celltypes not set"
    assert G.identities, f"identities not set"
    assert G.processes, f"processes not set"
    assert G.graph, f"graph not set"
    
    assert isinstance(G.celltypes, list), f"celltypes is not a list"
    assert isinstance(G.identities, dict), f"identities is not a dict"
    assert isinstance(G.processes, dict), f"processes is not a dict"
    assert isinstance(G.graph, nx.classes.graph.Graph), f"graph is not a networkx graph"
    
    # test methods (PLEASE ADD MORE USEFUL TESTS HERE)
    # assert that following class methods work: filter_edges, filter_nodes, get_celltype_hierarchy, get_celltype_processes, get_processes, get_identities

    assert isinstance(G.filter_edges(attribute_name="class", attributes=["SUBSET_OF"]), list), f"filter_edges does not return a list"
    assert isinstance(G.filter_nodes(attribute_name="gene_set_type", attributes=["manual_internal"]), list), f"filter_nodes does not return a list"
    assert isinstance(G.get_celltype_hierarchy(), dict), f"get_celltype_identity does not return a dict"
    pytest.raises(TypeError, G.get_celltype_hierarchy, celltypes=["B"])
    pytest.raises(ValueError, G.get_celltype_processes, celltypes=["B"])
    processes = G.get_processes(gene_sets=["all_macroautophagy_regulation_positive"])
    assert isinstance(processes, dict), f"get_processes does not return a dict"
    assert len(processes) > 0, f"get_processes does not return any processes"


def test_initialise_cytopus_default():
    """Test that default graph is loaded if no graph is provided"""
    G = cp.KnowledgeBase()
    assert G.graph.number_of_nodes() == 6113
    assert G.graph.number_of_edges() == 10785
