### This code calculates the energy gap between the mfe structure and first suboptimal structure of rna seqs and
### creates a plot and finally
### calculates the p-value for the energy gaps


#  1. Get energy gap between mfe and first subopt structure
import RNA

with open("natural_L20_seq_unique.txt",'r') as inp: lines=inp.readlines()
sequence_energy_gaps=[]

for line in lines:
    sequence = line[:-1]

    fc = RNA.fold_compound(sequence)
    mfe_structure, mfe_energy = fc.mfe()            # get mfe structure of rna seq and its free energy
    
    delta = 500  
    subopt_structures = fc.subopt(delta, 1)         # compute suboptimal structures within delta energy range
    
    while len(subopt_structures)<2 and delta<1000:  # make sure there are at least 2 suboptimal structures 
        delta+=200
        subopt_structures = fc.subopt(delta, 1)
    
    if subopt_structures[0].structure != mfe_structure or len(subopt_structures)==1: 
        energy_gap = subopt_structures[0].energy - mfe_energy  # use first subopt structure generated to calculate energy gap 
    else:
        energy_gap = subopt_structures[1].energy - mfe_energy  # if first subopt = mfe, use second subopt structure to calculate energy gap
    sequence_energy_gaps.append((sequence, energy_gap))        # add the rna seq and its energy gap to array


sequence_energy_gaps.sort(key=lambda x: x[1], reverse=True)    # sort the list based on energy gaps in descending order

with open("natural_energy_gap_L20.txt", 'w') as op:
    for seq, gap in sequence_energy_gaps:
        op.write(f"{seq} {gap}\n")                             # write seq and its energy gap to new txt file


# 2. Plot inset Graph (normal graph plus inset of log y axis) 
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

with open("rand_energy_gap_L20.txt",'r') as inp: lines1=inp.readlines()
with open("natural_energy_gap_20.txt",'r') as inp: lines2=inp.readlines()

rand_gap = [float(line.split()[1]) for line in lines1]        # get the energy gaps for random and natural seqs
natural_gap = [float(line.split()[1]) for line in lines2]

# Main Plot
fig, ax = plt.subplots(figsize=(10,8))
sns.kdeplot(rand_gap, label='Random', color='blue', bw_adjust=1.5, linewidth=2)  
sns.kdeplot(natural_gap, label='Natural', color='orange', bw_adjust=1.5, linewidth=2) 
ax.set_title('Density Distribution of Energy Gap for L = 20', fontsize=22, pad=18)
ax.set_xlabel("Energy Gap (kcal/mol)", fontsize=23, labelpad=12)
ax.set_ylabel("Density", fontsize=23, labelpad=12)
ax.tick_params(axis='both', labelsize=28)
ax.legend(fontsize=20, loc="lower right", bbox_to_anchor=(1, 0.1))

# Inset (log-y plot) 
ax_inset = inset_axes(ax, width="50%", height="50%", loc="upper right", borderpad=2.2)
sns.kdeplot(rand_gap, ax=ax_inset, color='blue', bw_adjust=1.5, linewidth=1.7, log_scale=(False, True))
sns.kdeplot(natural_gap, ax=ax_inset, color='orange', bw_adjust=1.5, linewidth=1.7, log_scale=(False, True))
ax_inset.set_xlabel("")
ax_inset.set_ylabel("")
ax_inset.tick_params(axis='both', labelsize=16)
ax_inset.set_title("Log-scale density", fontsize=14)

plt.show()


# 3. Generate P-values for Energy Gaps
from scipy.stats import ks_2samp

stat, pval = ks_2samp(natural_gap, rand_gap)               # Perform KS 2-sample test to get p-value
print(f'p-value: {pval:.3f}, Stat: {stat:.3f}')