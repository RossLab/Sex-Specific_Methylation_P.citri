#------------------------------------------------------
# Calling Transposable Elements in Planococcus citri
#------------------------------------------------------

# I followed the pipeline from https://github.com/SamuelHLewis/TEAnnotator
# where RepeatModeler and RepeatMasker are used to identify de novo TEs.

# The pipeline also calls for interproscan to be run which annotates
# TEs with protein domains from the Pfam database, you can then filter
# the de novo TEs and keep only ones with known TE domains. I didn't
# run this extra bit as Peter Sarkis said it wasn't needed because most
# of the de novo TEs do have the domains when he checked. 

# The final output is a .gff file with all of the TE annotations,
# you can then use bedtools or something to get a fasta of the TEs if needed.