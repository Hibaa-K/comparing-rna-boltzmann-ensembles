### This code calculates the relative entropy for combined shapes frequency
 

from scipy.stats import entropy

def read_shape_freqs(filepath):    # read shape frequencies from files
    shape_counts = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            shape = parts[0]       # shape of seq
            count = int(parts[1])  # frequency of seq
            shape_counts[shape] = count
    return shape_counts

def normalize(counter, keys):      # normalize to get probability distributions
    total = sum(counter.get(k, 0) for k in keys)
    return [counter.get(k, 0) / total if total > 0 else 0 for k in keys]

# MAIN
nat_file = "natural_lvl5_freq_5_combined_shapes_L20.txt"
rand_file = "rand_lvl5_freq_5_combined_shapes_L20.txt"

# Load frequency data
freq_natural = read_shape_freqs(nat_file)
freq_random = read_shape_freqs(rand_file)

all_shapes = set(freq_natural.keys()).union(set(freq_random.keys()))  # combine all unique shapes found in either distribution

# Get probability distributions aligned by shape
P = normalize(freq_natural, all_shapes)
Q = normalize(freq_random, all_shapes)

# Clip to avoid log(0) or division by zero in relative entropy
P = [max(p, 1e-10) for p in P]
Q = [max(q, 1e-10) for q in Q]

rel_entropy = entropy(P, Q, base=2)  # compute relative entropy (Kullback-Leibler divergence)
print(f"Relative Entropy: {rel_entropy:.2f} bits")