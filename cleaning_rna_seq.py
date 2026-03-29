### This file contains code to clean and get unique rna sequences for different lengths,
### removing the Schistosoma rna for the special case of samplig bias in L=45, 
### and lastly shuffling natural rna sequences.


#### ---- #### ---- ####
# 1. Converting fasta file to txt file (Getting only the RNA seq)
with open("natural_L20.fasta", 'r') as infile, open("natural_L20_seq.txt", 'w') as outfile: # get fasta file from RNAcentral database
    rna_sequence = ""
    for line in infile:
        if line.startswith(">"):               # skip this line in file since it does not contain a seq
            outfile.write(rna_sequence + "\n") # write only the seq to a new txt file
            rna_sequence = ""
        else:
            rna_sequence += line.strip()  
    outfile.write(rna_sequence + "\n")         # write down the last seq to new file


# 2. Getting only the unique sequences from the txt file
with open("natural_L20_seq.txt", 'r') as inp: lines=inp.readlines()
op= open("natural_L20_seq_uni.txt", 'w') 

unique_seq=set()               # create a set to store unique rna sequences 
for line in lines:
    if line not in unique_seq: 
        op.write(line)         # write seq to new file if unique
        unique_seq.add(line)   # add seq to set 
op.close()


# 3. Making sure all sequences are only ATCG nucleotides
with open("natural_L20_seq_uni.txt",'r') as inp: lines=inp.readlines()
nat_seq = [line.strip() for line in lines]
op = open("natural_L20_seq_unique.txt",'w') 

letters=['A', 'T', 'C', 'G']
s1= set(letters)           # create a set of expected nucleotides of rna
for seq in nat_seq:        # go through each seq to check the nuleotide composition
    nuc=set(list(seq))     # seperate the seq into a set of individual nucleotides
    if (nuc==s1):          
        op.write(seq+'\n') # write seqs comprised of only ATCG to a new txt file
    else:
        print(nuc)
op.close()


# 4. Generating unique random sequences
import random
op = open("rand_L20_seq_unique.txt",'w') 
length=20
count=12438                                    # same number as the count of unique natural sequences 
unique_seq=set()                               # create set to store unique seqs
while len(unique_seq) < count:
    random_rna_sequence = ''.join(random.choice('ATCG') for _ in range(length))  # create a random seq of required length
    if random_rna_sequence not in unique_seq:  # check if seq is unique
        unique_seq.add(random_rna_sequence)    # if unique, add to the set and write seq into a new txt file
        op.write(random_rna_sequence + "\n")   
op.close()


#### ---- #### ---- ####
# For the special case of sampling bias in L=45:
# Getting L=45 seq without Schistosoma Type RNA
with open("natural_L45.fasta", 'r') as infile, open("natural_L45_seq(No_Schistosoma).txt", 'w') as outfile:
    rna_sequence = ""
    skip_current = False
    count=0

    for line in infile:
        if line.startswith(">"):
            if rna_sequence and not skip_current: # write previous sequence if it wasn't skipped and isn't empty into new txt file
                outfile.write(rna_sequence + "\n")

            rna_sequence = ""                     # reset for new seq
            header = line.strip()
            skip_current = ("Schistosoma" in header and "Hammerhead ribozyme" in header) # skip these seqs
            if skip_current:
                count+=1
        else:
            if not skip_current:
                rna_sequence += line.strip()

    if rna_sequence and not skip_current:         # write last seq
        outfile.write(rna_sequence + "\n")
    print(count)


#### ---- #### ---- ####
# SHUFFLING NATURAL RNA (For same GC content)
import random

with open("natural_L20_seq_unique.txt", "r") as inp: lines=inp.readlines()
op = open("natural_L20_seq_shuffled.txt",'w') 

for seq in lines:
    seq = seq.strip()           
    nucleotides = list(seq)              # get all the nucleotides of the seq
    random.shuffle(nucleotides)          
    scrambled_seq = "".join(nucleotides) # join the shuffled nucleotides to create a new seq
    #print(scrambled_seq)
    op.write(scrambled_seq + "\n")       # write the shuffled seq to a new txt file
op.close()