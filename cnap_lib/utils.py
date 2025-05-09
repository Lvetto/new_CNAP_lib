import numpy as np

def load_xyz(path):
    """
    Load an XYZ file.
    Returns:
        elements: list of element symbols (str)
        positions: numpy array of shape (N, 3)
    """
    elements = []
    positions = []

    with open(path, "r") as f:
        lines = f.readlines()

    # Skip first two lines (atom count and comment)
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 4:
            elements.append(parts[0])
            positions.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return elements, np.array(positions)

def build_adjacency_matrix(positions, cutoff):
    """
    Build adjacency matrix based on Euclidean distance and a cutoff.
    
    Args:
        positions (np.ndarray): shape (N, 3)
        cutoff (float): distance threshold
    
    Returns:
        np.ndarray: adjacency matrix (bool or int) of shape (N, N)
    """
    from scipy.spatial.distance import cdist

    dist_matrix = cdist(positions, positions)
    adj_matrix = (dist_matrix < cutoff) & (dist_matrix > 0.0)
    return adj_matrix.astype(int)
