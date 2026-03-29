### This code has 3 main sections:
### 1. Generate suboptimal structures for an rna seq  
### 2. Combine the shapes of those suboptimals and find the frequency of each combined shape for every length of rna
### 3. Plot the frequency of combined shapes

### Note: to generate shapes of RNA please check the txt file "Steps to download RNAshapes on Linux" as it cant be done here in python


import RNA

# 1. Getting subopt structures of rna seqs
with open("natural_L20_seq_unique.txt",'r') as inp: lines=inp.readlines() # file of rna seqs whos subopts will be generated
op = open("natural_L20_mfe_struct.txt",'w')
op1= open("natural_L20_subopt1_struct.txt",'w')
op2= open("natural_L20_subopt2_struct.txt",'w')
op3= open("natural_L20_subopt3_struct.txt",'w')
op4= open("natural_L20_subopt4_struct.txt",'w')

for line in lines:
    sequence=line[:-1]
    fc = RNA.fold_compound(sequence)
    mfe = fc.mfe()                                    # get mfe structure of an rna seq
    op.write(mfe[0]+'\n')                             # write mfe seq to seperate file
    
    delta = 500  
    subopt_structures = fc.subopt(delta,1)            # compute suboptimal structures within delta energy range
    
    while len(subopt_structures) < 5 and delta<=2000: # make sure to generate at least 5 suboptimal structures
        delta += 200  
        subopt_structures = fc.subopt(delta,1)        # if not, increase delta and generate subopts again
    
    if len(subopt_structures)>4:
        if subopt_structures[0].structure == mfe[0]:  # ignore first subopt structure if its the same as mfe
            op1.write(subopt_structures[1].structure+'\n')
            op2.write(subopt_structures[2].structure+'\n')
            op3.write(subopt_structures[3].structure+'\n')
            op4.write(subopt_structures[4].structure+'\n')
        else:
            op1.write(subopt_structures[0].structure+'\n')
            op2.write(subopt_structures[1].structure+'\n')
            op3.write(subopt_structures[2].structure+'\n')
            op4.write(subopt_structures[3].structure+'\n')
    else:                                             # if no. of subopts < 5 
        op1.write(subopt_structures[0].structure+'\n')
        op2.write(subopt_structures[0].structure+'\n')
        op3.write(subopt_structures[0].structure+'\n')
        op4.write(subopt_structures[0].structure+'\n')
op.close()
op1.close()
op2.close()
op3.close()
op4.close() 


#### ---- #### ---- ####

# 2. Combining mfe and 4 subopts for each rna seq
with open("Shape_Lvl5_mfe_natural_L20.txt",'r') as inp: lines1=inp.readlines()
with open("Shape_Lvl5_subopt1_natural_L20.txt",'r') as inp: lines2=inp.readlines()
with open("Shape_Lvl5_subopt2_natural_L20.txt",'r') as inp: lines3=inp.readlines()
with open("Shape_Lvl5_subopt3_natural_L20.txt",'r') as inp: lines4=inp.readlines()
with open("Shape_Lvl5_subopt4_natural_L20.txt",'r') as inp: lines5=inp.readlines()

mfe=[]
for line in lines1:
    shape=line[:-1]
    mfe.append(shape)

subopt1=[]
for line in lines2:
    shape=line[:-1]
    subopt1.append(shape)

subopt2=[]
for line in lines3:
    shape=line[:-1]
    subopt2.append(shape)

subopt3=[]
for line in lines4:
    shape=line[:-1]
    subopt3.append(shape)

subopt4=[]
for line in lines5:
    shape=line[:-1]
    subopt4.append(shape)

op= open("natural_lvl5_combined_shapes_5_L20.txt",'w') 
for i in range(len(mfe)):                              # write combines subopt shapes of all rna seqs to a new txt file
    op.write(mfe[i]+"*"+subopt1[i]+"*"+subopt2[i]+"*"+subopt3[i]+"*"+subopt4[i]+'\n') 
op.close()


# Finding freq of combined phenotypes for each (do natural first & then random)
from collections import Counter
with open("natural_lvl5_combined_shapes_5_L20.txt",'r') as inp: lines=inp.readlines()

freq = Counter()
for line in lines:
    sh=line[:-1]
    freq[sh]+=1                                        # increment count for each combined shape

op = open("natural_lvl5_freq_5_combined_shapes_L20.txt",'w') 
for x in freq.most_common():           
    op.write(x[0]+" "+ str(x[1])+'\n')                 # write the combined shape with its frequency in a new txt file
op.close()


#### ---- #### ---- ####

# 3. Plotting the frequency of combined shapes
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt

with open("natural_lvl5_freq_5_combined_shapes_L20.txt",'r') as inp: lines=inp.readlines()
nat_shapes={}
for line in lines:
    parts = list(line.split())
    shape = parts[0]                                  # get the combined shape  
    frequency = int(parts[1])                         # get the frequency of the combined shape
    nat_shapes[shape] = frequency                     # add the freq to the shape for natural rna
   
with open("rand_lvl5_freq_5_combined_shapes_L20.txt",'r') as inp: lines=inp.readlines()
rand_shapes={}
for line in lines:
    parts = list(line.split())
    shape = parts[0]                                  # get the combined shape
    frequency = int(parts[1])                         # get the frequency of the combined shape
    rand_shapes[shape] = frequency                    # add the freq to the shape for random rna 

common_shapes = set(nat_shapes.keys()).intersection(set(rand_shapes.keys())) # find shapes in common between natural and random
natural = [nat_shapes[shape] for shape in common_shapes]                     # get frequency of the common shapes
rand = [rand_shapes[shape] for shape in common_shapes]
print(len(common_shapes))

# Get the Pearson correlation using log10 values
log_nat = np.log10(natural)                           
log_rand = np.log10(rand)
corr, pval = pearsonr(log_rand, log_nat)
print(f"p-value: {pval}, Pearson correlation: {corr:.3f}")


# Create scatter plot
plt.figure(figsize=(10,8))
plt.plot([1e-1, 1e5], [1e-1, 1e5], color='black', linewidth=3)
plt.scatter(x=rand, y=natural, color='blue', s=100)
plt.xlabel('Random', fontsize=23, labelpad=10)
plt.ylabel('Natural', fontsize=23, labelpad=10)
plt.title('Combined phenotype shape frequency L=20 (5)', fontsize=22)

plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.tight_layout()
plt.show()