from cnap_lib.utils import *
import numpy as np

class Graph:
    """
    A class to represent graphs from graph theory, with methods and functions designed to work with CNA techniques
    """
    def __init__(self, adj_mat, node_positions=[], fill_diag=True, root_node=None):
        """
        Create a graph object from an adjacency matrix.


        Args:
            adj_mat (np.array): adjacency matrix for the graph
            node_positions (np.array, optional): array of node positions
            fill_diag (bool, optional): set to False to avoid setting elements on the diagonal to zero. Defaults to True
            root_node (int, optional): the index of the node to be considered as root
        """

        self.adj_mat = adj_mat

        if (fill_diag):
            np.fill_diagonal(self.adj_mat, 0)

        self.adj_list = self.make_adj_list()
        self.positions = node_positions

        self.root_node = root_node

    def make_adj_list(self):
        """
        Create an adjacency list from the adj matrix of the graph

        Returns:
            list: the adjacency list
        """

        """
        Create an adjacency list from the adj matrix of the graph

        Returns:
            list: the adjacency list
        """

        return [np.where(row)[0] for row in self.adj_mat == 1]

    def make_subgraph(self, nodes):
        """
        Create a graph object from a subset of the graph's nodes and have it inherit the connections

        Args:
            nodes (list/np.array): a list of node indices

        Returns:
            Graph: the subgraph made from the given nodes
        """

        """
        Create a graph object from a subset of the graph's nodes and have it inherit the connections

        Args:
            nodes (list/np.array): a list of node indices

        Returns:
            Graph: the subgraph made from the given nodes
        """

        nodes = np.array(nodes)
        sub_adj_mat = self.adj_mat[np.ix_(nodes, nodes)] # make an adjacency matrix for the nodes of the subgraph by extracting elements i,j where both i and j are in the subgraph
        np.fill_diagonal(sub_adj_mat, 0)    # set the diagonal to 0. Avoids some issues when extracting bonds

        # return the subgraph, handles the case where the graph has no positions
        return Graph(sub_adj_mat, self.positions[nodes]) if (self.positions is not []) else Graph(sub_adj_mat)

    def longest_chain_from_node(self, starting_node):
        """
        Find the lenght of the longest chain of connected nodes, starting from starting_node.\n
        In most cases users should use the longest_chain property instead of calling this directly

        Args:
            starting_node (int): the index of the starting node

        Returns:
            int: the lenght of the longest chain starting from the starting node, measured in number of steps
        """

        # This is a slightly modified version of Depth First Search

        max_lenght = 0
        discovered = []
        stack = []

        stack.append((starting_node, 0))
        while (len(stack) > 0):
            v, l = stack[-1]
            stack = stack[:-1]

            if (v not in discovered):
                discovered.append(v)
                max_lenght = max(max_lenght, l)

                for w in self.adj_list[v]:
                    if ((l>=2) and (w == starting_node)):
                        max_lenght = max(max_lenght, l+1)

                    if (w not in discovered):
                        stack.append((w, l+1))

        return max_lenght

    def save(self, filepath):
        np.savez_compressed(filepath, adj_mat=self.adj_mat, pos=self.positions)

    def __repr__(self):
        return f"Number of nodes: {self.number_of_nodes}, number of bonds: {self.number_of_unique_bonds}"

    # Properties

    @property
    def longest_chain(self):
        """
        A property representing the lenght of the longest chain of connected nodes in the graph

        Returns:
            int: lenght of the longest chain starting from any node in the graph
        """

        """
        A property representing the lenght of the longest chain of connected nodes in the graph

        Returns:
            int: lenght of the longest chain starting from any node in the graph
        """

        # Note: there is a bug where lenght is calculated incorrectly for closed chains
        # should be easily solved by removing the starting node from the list of explored nodes, once the algorithm is far enough along

        max_lenght = 0
        for i,_ in enumerate(self.adj_list):
            max_lenght = max(self.longest_chain_from_node(i), max_lenght)

        return max_lenght

    @property
    def unique_bonds(self):
        """
        A property representing the set of unique bonds on the graph.
        Repetitions and pair-ordering are already handled by this property

        Returns:
            np.array: an array of unique pairs (i,j) where i,j are the indices of bonded atoms
        """

        """
        A property representing the set of unique bonds on the graph.
        Repetitions and pair-ordering are already handled by this property

        Returns:
            np.array: an array of unique pairs (i,j) where i,j are the indices of bonded atoms
        """

        pairs = np.argwhere(self.adj_mat == 1)   # extract (i,j) where a_ij=1
        sorted_pairs = np.sort(pairs, axis=1)   # sort each pair of indices
        unique_pairs = np.unique(sorted_pairs, axis=0)  # remove duplicates

        return unique_pairs

    @property
    def number_of_unique_bonds(self):
        """
        A property representing the number of unique bonds in the graph

        Returns:
            int: number of unique (considering pair-orderings) bonds in the graph
        """

        """
        A property representing the number of unique bonds in the graph

        Returns:
            int: number of unique (considering pair-orderings) bonds in the graph
        """

        return len(self.unique_bonds)

    @property
    def number_of_nodes(self):
        """
        A property representing the number of nodes in the graph

        Returns:
            int: number of nodes in the graph
        """

        """
        A property representing the number of nodes in the graph

        Returns:
            int: number of nodes in the graph
        """

        return self.adj_mat.shape[0]

def get_cns_subgraphs(graph):
    """
    Make a list of subgraphs representing common neighbors for each pair of bonded atoms

    Args:
        graph (Graph): the graph object representing the system

    Returns:
        list: a list of Graph object representing all subgraphs of common neighbors
    """

    subgraphs = [graph.make_subgraph(np.intersect1d(graph.adj_list[i], graph.adj_list[j])) for i,j in graph.unique_bonds]
    return subgraphs

def compute_signatures(graph):
    """
    Compute the CNA signatures for each pair of neighbors

    Args:
        graph (Graph): Graph object representing the system

    Returns:
        list: list of signatures for each bonds
    """

    subgraphs = get_cns_subgraphs(graph)

    signatures = []
    for subgraph in subgraphs:
        i = subgraph.number_of_nodes
        j = subgraph.number_of_unique_bonds
        k = subgraph.longest_chain
        signatures.append((i,j,k))

    return signatures

def get_occurrences(signatures):
    """
    Count occurrences for each unique signature in the system

    Args:
        signatures (list): a list of signatures. Usually the return value of compute_signatures

    Returns:
        list: a list of pairs. The first element is a signature, the second its number of occurrences
    """

    unique_sigs, counts = np.unique(signatures, axis=0, return_counts=True)

    return [(tuple(sig), count) for sig, count in zip(unique_sigs, counts)] #{tuple(sig): count for sig, count in zip(unique_sigs, counts)}

def get_CNAPs(graph):
    """
    Computes CNAPs from a given graph

    Args:
        graph (Graph): a graph representing the neighborhood relationships of the atoms

    Returns:
        list: a list of CNAPs (as lists of signatures, counts)
    """
    sigs = compute_signatures(graph)
    bonds = graph.unique_bonds

    sigs_per_node = [[] for _ in range(graph.number_of_nodes)]

    for sig, bond in zip(sigs, bonds):
        n1, n2 = bond
        sigs_per_node[n1].append(sig)
        sigs_per_node[n2].append(sig)

    out = [get_occurrences(i) for i in sigs_per_node]

    return out

def unique_CNAPs(CNAPs):
    """
    Creates a list of unique CNAPs from a list

    Args:
        CNAPs (list): a list of CNAPs

    Returns:
        list: a list of unique CNAPs
    """
    out = []
    for i in CNAPs:
        if (i not in out):
            out.append(i)

    return out

def group_CNAPs(CNAPs):
    """
    Groups CNAPs based on their type

    Args:
        CNAPs (list): a list of CNAPs

    Returns:
        list: a list where the first element is a set of unique CNAPs and the second is a list of feature vector's indices grouped by type
    """
    keys = unique_CNAPs(CNAPs)
    matching_indices = [np.array([i for i, comp in enumerate(CNAPs) if comp == key]) for key in keys]

    return [keys, matching_indices]

def make_neighbors_subgraph(graph, index):
    """
    A function to make a subgraph representing an atom's LAE from its index

    Args:
        graph (Graph): the graph representing the system
        index (int): the index of the atom in question

    Returns:
        Graph: a subgraph containing the given atom and its neighbours along with the neighborhood relationships
    """
    ns = list(graph.adj_list[index])
    ns.append(index)
    return graph.make_subgraph(ns)
