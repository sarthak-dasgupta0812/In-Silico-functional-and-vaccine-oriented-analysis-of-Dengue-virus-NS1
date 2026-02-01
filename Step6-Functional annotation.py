# Functional Annotation 
# Retrieving Uniprot Annotation 

from Bio import ExPASy
from Bio import SwissProt

# UniProt accession obtained from the top BLAST hit
# (Example accession for Dengue virus NS1)
uniprot_accession = "P29990"

# Retrieve UniProt record from ExPASy
handle = ExPASy.get_sprot_raw(uniprot_accession)
record = SwissProt.read(handle)

# Basic annotation fields
print("Protein name:", record.description)
print("Organism:", record.organism)
print("Sequence length:", record.sequence_length)


# extracting functional Information 
print("\nFunctional comments:")
for comment in record.comments:
    if comment.startswith("FUNCTION"):
        print(comment)


# Identifying Biological roleand Pathogenic relevance 
print("\nKeywords associated with NS1:")
for keyword in record.keywords:
    print(keyword)



# ------------------------------------------------------------
# Step 5: Functional Annotation – Interpretation
# ------------------------------------------------------------
# Functional annotation of the Dengue virus non-structural
# protein 1 (NS1) was performed using UniProt data derived
# from top BLAST homologs.
#
# The annotation indicates that NS1 plays a crucial role in
# viral replication and immune evasion, interacting directly
# with host immune pathways.
#
# Its conserved nature across Dengue virus serotypes and
# secretion into host circulation highlight its importance
# in viral pathogenicity and disease progression.
#
# This step provides biological meaning to the sequence and
# prepares it for structural analysis or vaccine-related studies.
# ------------------------------------------------------------
