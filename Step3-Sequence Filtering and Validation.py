# Sequence filtering & Validation 

from Bio import SeqIO

# Biological thresholds are defined to ensure that the Dengue virus
# non-structural protein 1 (NS1) sequence is complete and suitable
# for reliable functional and homology-based analysis.
MIN_PROTEIN_LENGTH = 100

# Ambiguous amino acids (X, B, Z, J) represent uncertain or mixed
# residues that can reduce the accuracy of BLAST searches and
# functional annotation, and are therefore filtered out.
AMBIGUOUS_AA = {"X", "B", "Z", "J"}

# The NS1 protein sequence is read from a FASTA file using Biopython
# to ensure proper parsing of biological sequence data and metadata.
record = SeqIO.read("NS1.fasta", "fasta")
sequence = str(record.seq).upper()

# Printing the sequence identifier allows traceability and confirms
# that the correct viral protein sequence is being evaluated.
print("Evaluating sequence:", record.id)

#  Criterion 1: Length check 
# Very short protein sequences may represent incomplete fragments
# or sequencing artifacts. NS1 is a functional viral protein with
# a known length range, so a minimum length threshold is applied.
if len(sequence) < MIN_PROTEIN_LENGTH:
    print("Sequence rejected: length below minimum threshold")
    decision = False
else:
    print("Length check passed")

#  Criterion 2: Ambiguous residue check 
# Ambiguous amino acids introduce uncertainty in sequence alignment
# and homology detection, which can compromise downstream analyses
# such as BLAST and functional annotation.
if any(aa in AMBIGUOUS_AA for aa in sequence):
    print("Sequence rejected: ambiguous amino acids detected")
    decision = False
else:
    print("No ambiguous amino acids detected")
    decision = True

# ---- Final decision ----
# Only sequences that pass all biological validation criteria are
# retained, ensuring high-confidence input for downstream analyses
# including homology search, structural prediction, and annotation.
if decision:
    print("Sequence retained for downstream analysis")
else:
    print("Sequence discarded")

# output - 
# Evaluating sequence: CAA78918.1
# Length check passed
# No ambiguous amino acids detected
# Sequence retained for downstream analysis
