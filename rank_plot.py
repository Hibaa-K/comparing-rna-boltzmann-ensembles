### This code generates suboptimal structures for an rna sequences and 
### calculates their Boltzmann probabilities and then
### creates a plot to show the mean probability of each rank of a suboptimal structure of a sequence.

import numpy as np
import RNA
import math
import random
import matplotlib.pyplot as plt


NUM_OF_SUBOPTS = 8  # number of suboptimal structures we will generate for each sequence
T = 310.15          # temperature in Kelvin (37°C) (for Boltzmann probability calculation)
kb = 0.001987       # boltzmann constant

def get_subopt(sequence):
    fc = RNA.fold_compound(sequence)                    
    mfe_structure, mfe_energy = fc.mfe()                 # get the mfe structure of a seq and its free energy
    
    delta = 500                                          # 500 = 5 kcal/mol (since delta * 0.01 kcal/mol)
    subopts = fc.subopt(delta,1)                         # generate suboptimal structures within the given delta energy range
    subopt_structures = []

    if subopts[0].structure!=mfe_structure:              # check if first suboptimal structure isnt the same as mfe 
        subopt_structures.append((mfe_structure, mfe_energy))  
        
    while len(subopts)<NUM_OF_SUBOPTS and delta <= 2000: # generate more subopt structures by increasing delta until there are at least the number you want
        delta += 300
        subopts = fc.subopt(delta,1)

    for i in range(min(NUM_OF_SUBOPTS, len(subopts))):   # add all the needed subopts to the array
        subopt_structures.append((subopts[i].structure, subopts[i].energy)) 

    return subopt_structures


def calculate_boltzmann_probabilities(structures):       # expects a tuple as parameter (subopt structure, energy) 
    free_energies = [struct[1] for struct in structures] # get all the free energies of the suboptimal sturctures
    boltzmann_factors = [math.exp(-E /(kb * T)) for E in free_energies]
    Z = sum(boltzmann_factors)
    probabilities = [bf / Z for bf in boltzmann_factors] # calculate Boltzmann probabilties for each
    return probabilities


# This function reads a file of rna seqs and gets a sample of seqs from it and calculates the mean Boltzmann probabilities for each subopt rank
def get_mean(file):                                      
    with open(file,'r') as inp: lines=inp.readlines()
    sequences = [line.strip() for line in lines]
    rna_sample = random.sample(sequences, 10000)         # creates a sample of 10,000 seqs

    ranks = [[] for _ in range(NUM_OF_SUBOPTS)]          # creates 2D array for probabilties of each rank 

    for seq in rna_sample:                               # get suboptimals and their boltz prob. for each seq in the sample
        subopts=get_subopt(seq)
        boltz_prob=calculate_boltzmann_probabilities(subopts)
        for i in range(NUM_OF_SUBOPTS):                  
            if i < len(boltz_prob):
                ranks[i].append(boltz_prob[i])           # add the prob. to the array of the corresponding rank
            else:
                ranks[i].append(np.nan)                  # if no. of subopts less than NUM_OF_SUBOPTS, fill missing ranks with NaN so it doesn't bias the mean
        
    means = [np.nanmean(rank) for rank in ranks]         # calculate mean of each rank
    stds = [np.nanstd(rank) for rank in ranks]           # calculate standard deviation of each rank 
    return means, stds


# MAIN
nat_file="natural_L20_seq_unique.txt"                    # (add path where file is saved)
rand_file="rand_L20_seq_unique.txt"
nat_means, nat_stds=get_mean(nat_file)                   # get values for natural seqs
rand_means, rand_stds=get_mean(rand_file)                # get values for random seqs


# Plot Boltzmann Probability vs Rank 
r = np.arange(1, NUM_OF_SUBOPTS + 1)                     # for X-axis
plt.figure(figsize=(10,8))
plt.errorbar(r, nat_means, yerr=nat_stds, fmt='o-', capsize=3, linewidth=2.5, color='orange', label="Natural")
plt.errorbar(r, rand_means, yerr=rand_stds, fmt='o-', capsize=3, linewidth=2.5, color='blue', label="Random")

plt.xlabel("Ranks", fontsize=23, labelpad=12)
plt.ylabel("Mean Boltzmann Probability", fontsize=23, labelpad=12)
plt.title("Rank Plot of RNA of L=20", fontsize=23)

plt.ylim(bottom=0)
plt.xticks(fontsize=26)
plt.yticks(fontsize=26)
plt.legend(fontsize=22)
plt.show()