### This code calculates the energy gap value threshold that separates the top 1% largest gaps from the bottom 99%
### and then plot graphs for them

import numpy as np

with open("rand_energy_gap_L20.txt",'r') as file:
    rand_gap = [float(line.split()[1]) for line in file]

with open("natural_energy_gap_L20.txt",'r') as file:
    natural_gap = [float(line.split()[1]) for line in file]

threshold = np.percentile(rand_gap, 99)                    # get the 99th percentile threshold of random rna 
print(f"Energy gap cutoff for top 1% in random RNA: {threshold:.2f}")

# Fraction of natural rna above this threshold 
count = 0 
for x in natural_gap:
    if x > threshold:
        count += 1
frac_above = count / len(natural_gap)

print(f"Natural RNA above threshold: {count} / {len(natural_gap)}")
print(f"Fraction of natural RNA above threshold: {frac_above:.2%}")


#### ---- #### ---- ####
# GRAPHS

# Plot for top 1% energy gap:
import matplotlib.pyplot as plt
import numpy as np

lengths = np.array([20, 30, 45, 60, 100, 150])
energy_thresholds = np.array([3.2, 2.3, 1.7, 1.3, 1.0, 0.7])  # energy thresholds for each rna length calculated from the code above

plt.figure(figsize=(10,8))
plt.plot(lengths, energy_thresholds, marker='o', linewidth=3)
plt.xlabel("RNA Sequence Length", fontsize=23, labelpad=10)
plt.ylabel("Energy Gap Threshold (kcal/mol)", fontsize=23, labelpad=10)
plt.title("Energy Gap Threshold (Top 1% of Random RNA)", fontsize=22)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.tight_layout()
plt.show()


# Plot for natural seq above this threshold:
import matplotlib.pyplot as plt
import numpy as np

lengths = np.array([20, 30, 45, 60, 100, 150])
percent_above = np.array([0.6, 0.8, 1.9, 1.8, 1.6, 1.7])      # percentage of natural rna above the threshold for each rna length

plt.figure(figsize=(10,8))
plt.plot(lengths, percent_above, marker='o', linewidth=3)
plt.axhline(y=1.0, linestyle='--', linewidth=2)               # reference line at 1%

plt.xlabel("RNA Sequence Length", fontsize=23, labelpad=10)
plt.ylabel("Percentage of RNA Above Threshold", fontsize=23, labelpad=10)
plt.title("Natural RNA Above Random Energy Gap Threshold", fontsize=22)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.tight_layout()
plt.show()