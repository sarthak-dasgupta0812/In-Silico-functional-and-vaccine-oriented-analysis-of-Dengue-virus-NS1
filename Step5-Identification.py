# Parsing Blast results with identifying homologs 
from Bio.Blast import NCBIXML

# Load the BLAST XML output
with open("NS1_blast_results.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)

# Number of homologous hits found
print("Total alignments found:", len(blast_record.alignments))

#Total alignments found: 50

# Extracting the top homolog - 
# The top alignment usually represents the closest homolog
top_alignment = blast_record.alignments[0]

print("Top Hit Description:")
print(top_alignment.title)

print("Alignment Length:", top_alignment.length)
# Total alignments found: 50
# Top Hit Description:
# emb|CAA78918.1| NS1, partial [dengue virus type 2]
# Alignment Length: 380

#High-Scoring Segment Pair (HSP) Analysis
# HSP represents the most conserved matching region
top_hsp = top_alignment.hsps[0]

print("Alignment Score:", top_hsp.score)
print("E-value:", top_hsp.expect)

print("Query sequence segment:")
print(top_hsp.query)

print("Matched sequence segment:")
print(top_hsp.sbjct)

print("Alignment pattern:")
print(top_hsp.match)


# Total alignments found: 50
# Top Hit Description:
# emb|CAA78918.1| NS1, partial [dengue virus type 2]
# Alignment Length: 380
# Alignment Score: 2049.0
# E-value: 0.0
# Query sequence segment:
# MNSRSTSLSVSQVLVGIVTLYLGVMVQADSGCVVSWKNKELKCGSGIFVTDNVHTRTEQYKFQPESPSKLASAIQKAHEEGICGIRSVTRLENLMWKQITSELNHILSENEVKLTIMTGDIKGIMQVGKRSLRPQPTELRYSWKTWGKAKMLSTELHNQTFLIDGPETAECPNTNRAWNSLEVEDYGFGVFTTNIWLRLREKQDAFCDSKLMSAAIKDNRAVHADMGYWIESALNDTWKIEKASFIEVKSCHWPKSHTLWSNGVLESEMVIPKNFAGPKSQHNNRPGYHTQTAGPWHLGKLEMDFDFCEGTTVVVTEDCGNRGPSLRTTTASGKLITEWCCRSCTLPPLRYRGEDGCWYGMEIRPLKEKEENLVSSLVTA
# Matched sequence segment:
# MNSRSTSLSVSQVLVGIVTLYLGVMVQADSGCVVSWKNKELKCGSGIFVTDNVHTRTEQYKFQPESPSKLASAIQKAHEEGICGIRSVTRLENLMWKQITSELNHILSENEVKLTIMTGDIKGIMQVGKRSLRPQPTELRYSWKTWGKAKMLSTELHNQTFLIDGPETAECPNTNRAWNSLEVEDYGFGVFTTNIWLRLREKQDAFCDSKLMSAAIKDNRAVHADMGYWIESALNDTWKIEKASFIEVKSCHWPKSHTLWSNGVLESEMVIPKNFAGPKSQHNNRPGYHTQTAGPWHLGKLEMDFDFCEGTTVVVTEDCGNRGPSLRTTTASGKLITEWCCRSCTLPPLRYRGEDGCWYGMEIRPLKEKEENLVSSLVTA
# Alignment pattern:
# MNSRSTSLSVSQVLVGIVTLYLGVMVQADSGCVVSWKNKELKCGSGIFVTDNVHTRTEQYKFQPESPSKLASAIQKAHEEGICGIRSVTRLENLMWKQITSELNHILSENEVKLTIMTGDIKGIMQVGKRSLRPQPTELRYSWKTWGKAKMLSTELHNQTFLIDGPETAECPNTNRAWNSLEVEDYGFGVFTTNIWLRLREKQDAFCDSKLMSAAIKDNRAVHADMGYWIESALNDTWKIEKASFIEVKSCHWPKSHTLWSNGVLESEMVIPKNFAGPKSQHNNRPGYHTQTAGPWHLGKLEMDFDFCEGTTVVVTEDCGNRGPSLRTTTASGKLITEWCCRSCTLPPLRYRGEDGCWYGMEIRPLKEKEENLVSSLVTA
# PS C:\Users\ADMIN\OneDrive\Desktop\Dengue virus polyprotein> 



# ------------------------------------------------------------
# Step 4: Homology Search (BLAST) – Interpretation
# ------------------------------------------------------------
# BLASTP was performed to identify homologous proteins for the
# Dengue virus non-structural protein 1 (NS1).
#
# The analysis revealed significant similarity to NS1 proteins
# from related Dengue virus strains and Flaviviruses, supported
# by high alignment scores and low E-values.
#
# Conserved regions identified through HSP analysis suggest
# essential functional and structural roles for NS1.
#
# These results provide strong evolutionary and functional
# evidence, enabling confident functional annotation in
# subsequent steps.


