### This code computes shannon entropy of rna sequences and
### plots the entropies of natural & random and 
### calculates the p-values of the entropy values

import RNA
import math
import numpy as np
from scipy.stats import entropy
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Function to compute Shannon entropy for an rna sequence
def compute_shannon_entropy(seq, delta=500, exclude_mfe=False):
    fc = RNA.fold_compound(seq)
    mfe_structure, mfe_energy = fc.mfe()    # get the mfe structure of rna and its free energy
    subopt_structs = fc.subopt(delta, 1)    # generate its suboptimal sturctures within delta energy range of 5 kcal/mol

    energies = []
    if subopt_structs[0].structure!=mfe_structure and exclude_mfe==False: # include the mfe structure of the seq
        energies.append(mfe_energy)
    for s in subopt_structs:
        if exclude_mfe and (s.structure == mfe_structure):                # ignore mfe structure 
            continue
        energies.append(s.energy)                                         # add all suboptimal energies to array

    # Compute Boltzmann Probability
    kb = 0.001987                                               # Boltzmann constant
    temp = 310.15                                               # temperature in Kelvin
    RT = kb * temp
    boltzmann_factors = [math.exp(-E / RT) for E in energies]   # Boltz Prob = exp(E/RT) / Z
    Z = sum(boltzmann_factors)
    probabilities = [bf / Z for bf in boltzmann_factors]

    H = entropy(probabilities, base=2)                          # calculate Shannon entropy 
    return H

# MAIN 
# (do for natural then random if needed)
with open("natural_L20_seq_unique.txt",'r') as inp: lines=inp.readlines()
sequences = [line.strip() for line in lines]
op = open("shannon_entropy_L20.txt",'w') 

entropies_including = []
entropies_excluding = []

for seq in sequences:
    e_incl = compute_shannon_entropy(seq, exclude_mfe=False)   # calculate shannon entrpy including mfe structure of seq
    e_excl = compute_shannon_entropy(seq, exclude_mfe=True)    # calculate shannon entrpy excluding mfe structure of seq
    entropies_including.append(e_incl)
    entropies_excluding.append(e_excl)

    op.write(f"{seq} {e_incl:.3f} {e_excl:.3f}\n")             # write seq and the shannon entropies (including and excluding mfe) to a new txt file

# Compute mean entropies
mean_incl = np.mean(entropies_including) 
mean_excl = np.mean(entropies_excluding)
op.close()

print(f"Mean entropy (including MFE): {mean_incl:.2f}")
print(f"Mean entropy (excluding MFE): {mean_excl:.2f}")


#### ---- #### ---- ####

# 2. Plot distribution of entropies
nat_entropy = []
rand_entropy = []
# shuffled_entropy=[]               # add for shuffled plot
with open("shannon_entropy_natural_L20.txt",'r') as file:
    for line in file:
        parts = line.strip().split()
        entpy = float(parts[1])     # entropy inlcuding mfe(1)/excluding mfe(2)
        nat_entropy.append(entpy)

with open("shannon_entropy_rand_L20.txt",'r') as file:
    for line in file:
        parts = line.strip().split()
        entpy = float(parts[1])     # entropy inlcuding mfe(1)/excluding mfe(2)
        rand_entropy.append(entpy)

## add this for shuffled natural rna
'''with open("shannon_entropy_L20_shuffled.txt",'r') as file:
    for line in file:
        parts = line.strip().split()
        entpy = float(parts[1])     # entropy inlcuding(1)/excluding(2) mfe
        shuffled_entropy.append(entpy)'''

# Create combined density plot of both natural and random
plt.figure(figsize=(10,8))
sns.kdeplot(rand_entropy, color='blue', label="Random", bw_adjust=1.5, linewidth=2)
sns.kdeplot(nat_entropy, color='orange', label="Natural", bw_adjust=1.5, linewidth=2)
#sns.kdeplot(shuffled_entropy, color='black', label="Shuffled", bw_adjust=1.5, linewidth=2, linestyle='--')  # add this line for shuffled plot

plt.xlabel("Shannon Entropy", fontsize=23, labelpad=8)
plt.ylabel("Density", fontsize=23, labelpad=8)
plt.title("Natural Mean Entropy: 2.2, Random Mean Entropy: 2.1", fontsize=22, pad=10)
plt.suptitle("L=20 - Including MFE Structure", fontsize=20)   
plt.legend(fontsize=20)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.tight_layout()
plt.show()


#### ---- #### ---- ####

# 3. Generate P-value for Shannon Entropy
from scipy.stats import ks_2samp

nat_entropy = []
rand_entropy = []
with open("shannon_entropy_natural_L20.txt",'r') as file:
    for line in file:
        parts = line.strip().split()
        entpy = float(parts[1])     # entropy inlcuding mfe(1)/excluding mfe(2) 
        nat_entropy.append(entpy)
        
with open("shannon_entropy_rand_L20.txt",'r') as file:
    for line in file:
        parts = line.strip().split()
        entpy = float(parts[1])     # entropy inlcuding mfe(1)/excluding mfe(2)
        rand_entropy.append(entpy) 

stat, pval = ks_2samp(nat_entropy, rand_entropy)  # perform KS 2-sample test
print(f'p-value: {pval:.3e}, Stat: {stat:.3f}')