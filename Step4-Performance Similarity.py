from Bio.Blast import NCBIWWW
from Bio import SeqIO

# Read the validated NS1 protein sequence
record = SeqIO.read("NS1.fasta", "fasta")

# BLASTP is used to compare the NS1 protein against
# the non-redundant protein database to identify homologous proteins.
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=record.seq
)

# Save BLAST results locally in XML format for downstream parsing
with open("NS1_blast_results.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

print("BLAST search completed and results saved.")
