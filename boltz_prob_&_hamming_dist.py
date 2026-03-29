### This code generates the suboptimal structures of rna sequences and calculates thier Boltzmann probabilities 
### then for a given sequence in a sample, two structures are sampled from its Boltzmann ensemble to find the Hamming distance between them. 
### The mean hamming distance for the samples of each length of rna are then plotted.

import random
import numpy as np
import RNA
import math

kb = 0.001987  # boltzmann constant
T = 310.15     # temperature in Kelvin (37°C)

def get_subopt(sequence):
    fc = RNA.fold_compound(sequence)                    
    mfe_structure, mfe_energy = fc.mfe()    # get the mfe structure of a seq and its free energy
    
    delta = 500                             # 500 = 5 kcal/mol (since delta * 0.01 kcal/mol)
    subopts = fc.subopt(delta,1)            # generate suboptimal structures within the given delta energy range
    subopt_structures = []

    if subopts[0].structure!=mfe_structure: # check if first suboptimal structure isnt the same as mfe 
        subopt_structures.append((mfe_structure, mfe_energy))  

    for i in range(len(subopts)):           # add all the needed subopts to the array
        subopt_structures.append((subopts[i].structure, subopts[i].energy)) 

    return subopt_structures

def calculate_boltzmann_probabilities(structures):        # expects a tuple as parameter (subopt structure, energy) 
    free_energies = [struct[1] for struct in structures]  # get all the free energies of the suboptimal sturctures
    boltzmann_factors = [math.exp(-E /(kb * T)) for E in free_energies]
    Z = sum(boltzmann_factors)
    probabilities = [bf / Z for bf in boltzmann_factors]  # calculate Boltzmann probabilties for each

    prob=[]
    for i in range(len(structures)):
        prob.append((structures[i][0], probabilities[i])) # add the subopt structure and its corresponding probability to array
    return prob

def hamming_distance(seq1, seq2):          # calculate hamming distance between 2 sequences              
    if len(seq1)==len(seq2):
        count=0
        for i in range(len(seq1)):
            if seq1[i]!=seq2[i]:
                count+=1
        return count


# MAIN (single seq comparision) (do natural and then random seqs)
with open("natural_L20_seq_unique.txt",'r') as inp: lines1=inp.readlines()
sequences_nat = [line[:-1] for line in lines1]

sample_size=10000
rna_sample = random.sample(sequences_nat, sample_size)   # create a sample of 10,000 sequences

tot_count=0
distances=[]
for seq in rna_sample:
    tot_count+=1

    # Finding their suboptimal structures and Boltzmann probabilities
    subopts_seq = get_subopt(seq)
    boltzmann_probs_seq = calculate_boltzmann_probabilities(subopts_seq)
    
    # Randomly picking a suboptimal structure based on Boltzmann probabilities
    subopts=[]
    pValues=[]
    for struct,prob in boltzmann_probs_seq:
        subopts.append(struct)
        pValues.append(prob)
    chosen_struct1 = np.random.choice(subopts, p=pValues)
    chosen_struct2 = np.random.choice(subopts, p=pValues)
    print(chosen_struct1, chosen_struct2)

    #Computing Hamming distance
    distance = hamming_distance(chosen_struct1, chosen_struct2)
    print(tot_count, "Hamming distance", distance)
    distances.append(distance)

# Compute mean & standard deviation of hamming distance
mean_distance = np.mean(distances)
std_distance = np.std(distances)
print(f"Mean Hamming Dist:")
print(f"{mean_distance:.2f} (+/- {std_distance:.2f})")


#### ---- #### ---- ####

## Plot for mean hamming distance data 
import matplotlib.pyplot as plt
import numpy as np

lengths = np.array([20, 30, 45, 60, 100, 150])  # rna lengths used

# Mean Hamming distances (gotten from the code above)
natural_mean = np.array([3.3, 6.1, 6.7, 10.6, 16.6, 24.4])
random_mean  = np.array([3.2, 5.6, 9.3, 12.8, 21.5, 32.0])
# Standard deviations
natural_std = np.array([3.7, 5.8, 6.6, 10.0, 14.1, 18.6])
random_std  = np.array([3.7, 5.5, 8.2, 10.5, 15.7, 21.7])

# Horizontal offset for visibility 
delta = 1
lengths_nat = lengths - delta
lengths_rand = lengths + delta

plt.figure(figsize=(10,8))
plt.errorbar(lengths_nat, natural_mean, color='orange', yerr=natural_std, fmt='o-', capsize=5, linewidth=2.5, label='Natural')
plt.errorbar(lengths_rand, random_mean, color='blue', yerr=random_std, fmt='o-', capsize=5, linewidth=2.5,  label='Random')

plt.xlabel("RNA Sequence Length", fontsize=23, labelpad=10)
plt.ylabel("Mean Hamming Distance", fontsize=23, labelpad=10)
plt.title("Mean Self Hamming Distance", fontsize=22)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.tight_layout()
plt.legend(fontsize=20)
plt.show()