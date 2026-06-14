import sys
import matplotlib.pyplot as plt

def parse_fasta(file_path):
    sequences = {}
    current_label = None
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_label = line[1:]
                sequences[current_label] = ""
            elif current_label:
                sequences[current_label] += line.replace(" ", "")
    return sequences

def levenshtein_distance(seq1, seq2, gap_penalty=1.0, mismatch_penalty=1.0):
    n, m = len(seq1), len(seq2)
    dp = [[0.0 for _ in range(m + 1)] for _ in range(n + 1)]
    
    for i in range(n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(m + 1):
        dp[0][j] = j * gap_penalty
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0.0 if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
            dp[i][j] = min(dp[i - 1][j] + gap_penalty,
                           dp[i][j - 1] + gap_penalty,
                           dp[i - 1][j - 1] + cost)
                           
    max_len = max(n, m)
    if max_len == 0:
        return 0.0
    return dp[n][m] / max_len

def calculate_distance_matrix(sequences, gap_penalty=1.0, mismatch_penalty=1.0):
    labels = list(sequences.keys())
    n = len(labels)
    matrix = {}
    for i in range(n):
        matrix[labels[i]] = {}
        for j in range(n):
            if i == j:
                matrix[labels[i]][labels[j]] = 0.0
            elif i < j:
                dist = levenshtein_distance(sequences[labels[i]], sequences[labels[j]], gap_penalty, mismatch_penalty)
                matrix[labels[i]][labels[j]] = dist
            else:
                matrix[labels[i]][labels[j]] = matrix[labels[j]][labels[i]]
    return matrix

class UPGMANode:
    def __init__(self, name, size=1, left=None, right=None, distance=0.0):
        self.name = name
        self.size = size
        self.left = left
        self.right = right
        self.distance = distance
        self.height = 0.0

def find_minimum_distance(matrix, current_clusters):
    min_dist = float('inf')
    pair = None
    for i in range(len(current_clusters)):
        for j in range(i + 1, len(current_clusters)):
            c1 = current_clusters[i]
            c2 = current_clusters[j]
            if matrix[c1][c2] < min_dist:
                min_dist = matrix[c1][c2]
                pair = (c1, c2)
    return pair, min_dist

def run_upgma(distance_matrix):
    clusters = list(distance_matrix.keys())
    cluster_sizes = {c: 1 for c in clusters}
    nodes = {c: UPGMANode(name=c, size=1) for c in clusters}
    matrix = {c1: {c2: distance_matrix[c1][c2] for c2 in distance_matrix[c1]} for c1 in distance_matrix}
    
    step = 1
    while len(clusters) > 1:
        (c1, c2), min_dist = find_minimum_distance(matrix, clusters)
        
        new_cluster_name = f"({c1},{c2})"
        new_cluster_size = cluster_sizes[c1] + cluster_sizes[c2]
        
        new_node = UPGMANode(name=new_cluster_name, size=new_cluster_size)
        new_node.left = nodes[c1]
        new_node.right = nodes[c2]
        
        height = min_dist / 2.0
        new_node.height = height
        nodes[c1].distance = height - nodes[c1].height
        nodes[c2].distance = height - nodes[c2].height
        
        nodes[new_cluster_name] = new_node
        
        matrix[new_cluster_name] = {}
        for c in clusters:
            if c != c1 and c != c2:
                d_c1 = matrix[c1][c]
                d_c2 = matrix[c2][c]
                s_c1 = cluster_sizes[c1]
                s_c2 = cluster_sizes[c2]
                
                new_dist = (s_c1 * d_c1 + s_c2 * d_c2) / (s_c1 + s_c2)
                matrix[new_cluster_name][c] = new_dist
                matrix[c][new_cluster_name] = new_dist
                
        matrix[new_cluster_name][new_cluster_name] = 0.0
        
        clusters.remove(c1)
        clusters.remove(c2)
        clusters.append(new_cluster_name)
        
        cluster_sizes[new_cluster_name] = new_cluster_size
        step += 1
        
    root_cluster = clusters[0]
    return nodes[root_cluster]

def draw_dendrogram(root_node, show=True):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    leaf_y_coords = {}
    current_y = 0
    
    def calculate_leaf_positions(node):
        nonlocal current_y
        if node.left is None and node.right is None:
            leaf_y_coords[node.name] = current_y
            current_y += 1
            return leaf_y_coords[node.name]
        
        y_left = calculate_leaf_positions(node.left)
        y_right = calculate_leaf_positions(node.right)
        
        return (y_left + y_right) / 2.0
        
    calculate_leaf_positions(root_node)
    
    def draw_branches(node):
        if node.left is None and node.right is None:
            return leaf_y_coords[node.name], node.height
            
        y_left, x_left = draw_branches(node.left)
        y_right, x_right = draw_branches(node.right)
        
        y_center = (y_left + y_right) / 2.0
        x_center = node.height
        
        ax.plot([x_left, x_center], [y_left, y_left], 'k-', lw=2)
        ax.plot([x_right, x_center], [y_right, y_right], 'k-', lw=2)
        ax.plot([x_center, x_center], [y_left, y_right], 'k-', lw=2)
        
        return y_center, x_center
        
    draw_branches(root_node)
    
    ax.set_yticks([leaf_y_coords[name] for name in leaf_y_coords])
    ax.set_yticklabels(list(leaf_y_coords.keys()), fontsize=12)
    
    ax.set_xlim(root_node.height * 1.1, -0.05)
    ax.set_xlabel("Evolutionary Distance (Differences)", fontsize=12, fontweight='bold')
    ax.set_title("Custom UPGMA Phylogenetic Dendrogram", fontsize=14, fontweight='bold')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.grid(axis='x', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    if show:
        plt.show()
    return fig

def main(fasta_file):
    print(f"\n--- 1. Input Module ---")
    print(f"Loading sequences from {fasta_file}...")
    sequences = parse_fasta(fasta_file)
    print(f"Found {len(sequences)} sequences: {list(sequences.keys())}")
    
    print("\n--- 2 & 3. Alignment and Distance Matrix Module ---")
    distance_matrix = calculate_distance_matrix(sequences)
    
    labels = list(sequences.keys())
    header = "      " + " ".join([f"{l[:5]:>7}" for l in labels])
    print("Initial Distance Matrix ($D_0$):")
    print(header)
    for l1 in labels:
        row = [f"{l1[:5]:>5}"]
        for l2 in labels:
            row.append(f"{distance_matrix[l1][l2]:>7.3f}")
        print(" ".join(row))
        
    print("\n--- 4. UPGMA Clustering Engine ---")
    print("Running proportional weighted clustering...")
    root_node = run_upgma(distance_matrix)
    print(f"Root cluster formed: {root_node.name}")
    print(f"Total Tree Height (Max distance): {root_node.height:.3f}")
    
    print("\n--- 5. Visualization Module ---")
    print("Rendering final dendrogram...")
    draw_dendrogram(root_node)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        print("Usage: python custom_upgma.py <fasta_file>")