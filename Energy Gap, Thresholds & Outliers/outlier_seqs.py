### This code has 3 sections:
### 1. Get the types of the rna seq from the original fasta file 
### 2. Calcultes the IQR for the outliers sequences
### 3. Find the RNA types of the outlier sequences
 

# 1. Get RNA types of all seqs of a certain length from fasta file
import numpy as np

seq_to_type = {}
header = None
sequence = ''                                   # will combine multiple sequence lines from fasta file for longer seqs
with open("natural_L20.fasta", "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):                # if we were building a sequence before this header, save it first
            if header is not None and sequence:
                full_seq = sequence
                if full_seq not in seq_to_type: # making sure seq is unique
                    seq_to_type[full_seq] = header

            header = line                       # new entry
            sequence = ''
        else:
            sequence+=line                      # add rna seq

    if header is not None and sequence:         # save last sequence after file ends
        full_seq = sequence
        if full_seq not in seq_to_type:
            seq_to_type[full_seq] = header


#### ---- #### ---- ####

# 2. Calculating IQR and upper cutoff for outlier seqs
def compute_upper_cutoff(data):
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    upper = Q3 + 1.5 * IQR
    return upper

with open("natural_energy_gap_L20.txt",'r') as inp: lines1=inp.readlines()  # energy gap files that were already created
natural_seq = [line.split()[0] for line in lines1]                          # get seq and its energy gap from the file
natural_gap = [float(line.split()[1]) for line in lines1]

upper_nat = compute_upper_cutoff(natural_gap)
num_outliers=np.sum(natural_gap > upper_nat)                                # number of outliers

op = open("natural_outliers(energy_gap)_L20.txt",'w') 
natural_outliers = []
for seq, gap in zip(natural_seq, natural_gap):
    if gap > upper_nat:                                                    # check if seq is an outlier
        seq_type = seq_to_type.get(seq)                                    # get the type of the outlier seq
        natural_outliers.append((seq, gap, seq_type))
        op.write(f'{seq} {gap}\n')                                         # write outlier seq to a new txt file
        op.write(seq_type+'\n')
op.close()

print(f'Natural: Upper Cutoff={upper_nat:.3} and no. of outliers={num_outliers}')


#### ---- #### ---- ####

# 3. Finding RNA Types (Composition) of Outlier Sequences 
import re
from collections import Counter

type_counts = Counter()
unknown_headers = []

with open('natural_outliers(energy_gap)_L20.txt', "r") as f:
    for line in f:
        line = line.strip()
        if not line.startswith(">"):   # only looking for headers (RNA Types)
            continue
        header = line
        h = header.lower()

        rna_type = "unknown"           # classification of rna seq

        # piRNA 
        if re.search(r"\bpir[-_]?\d", h) or "_pir_" in h or "pirna" in h:
            rna_type = "piRNA"

        # microRNA (miRNA)
        elif re.search(r"\bmir[-]?[a-z0-9]+", h) or "microrna" in h or "mirna" in h:
            rna_type = "miRNA"

        # ribosomal RNA
        elif ("ribosomal rna" in h) or re.search(r"\brrna\b", h) or re.search(r"\b(5\.8s|16s|18s|23s|28s)\b", h):
            rna_type = "rRNA"

        # tRNA
        elif "trna" in h or "transfer rna" in h:
            rna_type = "tRNA"

        # transfer-messenger RNA:
        elif "transfer-messenger rna" in h or "tmrna" in h or "partial ssra" in h:
            rna_type = "tmRNA"
        
        # signal recognition particle RNA
        elif "signal recognition particle" in h or "srp rna" in h or "small srp" in h:
            rna_type = "SRP_RNA"

        # crispr RNA
        elif "crispr rna" in h or "crispr" in h:
            rna_type = "crRNA"
        
        # ribosomal leader
        elif "ribosomal protein leader" in h:
            rna_type = "r_leader_RNA"
        
        # spliced leader
        elif "spliced leader" in h or re.search(r"\bsl\b", h):
            rna_type = "SL_RNA"

        # small nucleolar RNA
        elif "snorna" in h or "small nucleolar rna" in h:
            rna_type = "snoRNA"

        # small nuclear RNA 
        elif "snrna" in h or "small nuclear rna" in h or "spliceosomal rna" in h or re.search(r"\bu(1|2|4|5|6|11|12)\b", h):
            rna_type = "snRNA"

        # clustered small RNA (specific)
        elif "oocyte_clustered_small_rna" in h:
            rna_type = "oocyte_clustered_small_RNA"

        # small RNA 
        elif re.search(r"\bsrna\b", h) or "small rna" in h:
            rna_type = "small_RNA"

        # ribozyme
        elif "ribozyme" in h:
            rna_type = "ribozyme"

        # riboswitch
        elif "riboswitch" in h or "fmn sequence" in h or "amp sequence" in h or "adocbl sequence" in h:
            rna_type = "riboswitch"

        # RNase_P
        elif "rnase_p" in h:
            rna_type = "RNase_P_RNA"
        
        # regulatory
        elif "pemk rna" in h:
            rna_type = "regulatory_RNA"
        
        # Internal Transcribed Spacer (ITS)
        elif "internal transcribed spacer" in h or re.search(r"\bits[12]\b", h):
            rna_type = "ITS"
        
        # Group II catalytic intron
        elif "group ii catalytic intron" in h or "group-ii" in h or re.search(r"group\s*ii", h):
            rna_type = "Group_II_intron"

        # Y RNA
        elif "y rna" in h or re.search(r"\bhy\d\b", h):
            rna_type = "Y_RNA"
        
        # Vault RNA
        elif "vault rna" in h:
            rna_type = "Vault_RNA"

        # miscellaneous RNA
        elif "miscellaneous rna" in h or "misc rna" in h:
            rna_type = "miscRNA"

        # long non-coding RNA
        elif "lncrna" in h or "long non" in h:
            rna_type = "lncRNA"

        # non-coding RNA
        elif "ncrna" in h or "npcr" in h or "noncoding" in h or "non-coding" in h:
            rna_type = "other_ncRNA"

        type_counts[rna_type] += 1
        if rna_type == "unknown":
            unknown_headers.append(header)

total = sum(type_counts.values())
print(f"\nTotal types found: {total}")

print("\nRNA Type Composition:")
for t, c in type_counts.most_common():
    pct = (c / total * 100) if total > 0 else 0    # calculate percentage of each compostion 
    print(f"{t:12s} : {c:5d} ({pct:.1f}%)")


# Write unknown headers to a new txt file 
with open("outliers(energy_gap)_Unknown_Types_L20.txt", "w") as out:
    out.write("Unknown headers\n")
    for h in unknown_headers:
        out.write(h + "\n")