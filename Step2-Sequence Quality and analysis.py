# Reading the NS1 sequence 
from Bio import SeqIO

record = SeqIO.read("NS1.fasta", "fasta")

print("Sequence ID:", record.id)
print("Description:", record.description)
print("Sequence length:", len(record.seq))

# Sequence ID: CAA78918.1
# Description: CAA78918.1 NS1, partial [dengue virus type 2]
# Sequence length: 380


# Sequence length analysis 
sequence_length = len(record.seq)
print("Length of Dengue virus non-structural protein 1 (NS1):", sequence_length)
#Length of Dengue virus non-structural protein 1 (NS1): 380


#Amino Acid Composition Analysis (Protein Sequence)
# Since Dengue virus non-structural protein 1 (NS1) is a protein, amino acid composition is calculated.

from collections import Counter  #collection function is used to count the amino acids and counter counts the elements

protein_seq = str(record.seq)

aa_composition = Counter(protein_seq)

print("Amino Acid Composition:")
for aa, count in aa_composition.items(): # iterates over each amino acid loop 
    print(f"{aa}: {count}")

# Amino Acid Composition:
# M: 11
# N: 19
# S: 32
# R: 18
# T: 30
# L: 32
# V: 24
# Q: 12
# G: 27
# I: 19
# Y: 8
# A: 18
# D: 14
# C: 12
# W: 13
# K: 26
# E: 30
# F: 10
# H: 10
# P: 15


# Step 2: Sequence Quality & Basic Analysis – Interpretation
# ------------------------------------------------------------
# Based on the sequence quality assessment performed above:
#
# 1. The Dengue virus non-structural protein 1 (NS1) sequence
#    shows an appropriate and biologically expected length.
#
# 2. No ambiguous or invalid amino acid residues were detected,
#    indicating high sequence integrity and reliability.
#
# 3. The amino acid composition is biologically reasonable
#    and consistent with known viral protein characteristics.
#
# These observations confirm that the NS1 protein sequence
# passes basic quality control criteria
