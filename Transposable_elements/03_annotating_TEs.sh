# ----------------------------------------------------------------
# Annotating the TEs with some useful name
# ----------------------------------------------------------------

# RepBase is now a paid for service but through Edinburgh can get the database I need
# (ref Kamil for command)
/ceph/software/repeatmasker/RepeatMasker-4.1.0/util/queryRepeatDatabase.pl \
-species Arthropoda > Arthropoda.repeat.lib.RepBase.2020.fa

# ----------------------------------------------------------------
# Make a fasta of the TEs to annotate

cut -f1,4,5,9 PCITRI_v0_TEs.gff > PCITRI_v0_TEs.bed
bedtools getfasta -fi ../citri/PCITRI.assembly.v0.fa -bed PCITRI_v0_TEs.bed > PCITRI_v0_TEs.fa

# Make blast databases (urgh checky duplicated TE that needs removing)
sed  '/PCITRI_00251:78779-78800/,+1d' PCITRI_v0_TEs.fa > PCITRI_v0_TEs_dedup.fa

/ceph/software/blast/ncbi-blast-2.8.0+/bin/makeblastdb \
-in PCITRI_v0_TEs_dedup.fa \
-dbtype nucl \
-parse_seqids \
-out PCITRI_v0_TEs


# XXX STUCK HERE this is apparently a dedup seq id but I can't seem to find it or remove it :S
sed  '/NMAR\-15\_HMEL\#DNA/,+1d' Arthropoda.repeat.lib.RepBase.2020.fa > Arthropoda.repeat.lib.RepBase.2020_dedup.fa

/ceph/software/blast/ncbi-blast-2.8.0+/bin/makeblastdb \
-in Arthropoda.repeat.lib.RepBase.2020_dedup.fa \
-dbtype nucl \
-parse_seqids \
-out Arthropoda_TEs

# NOTE: need to do this for all TEs actually as need to check for enrichment of diff meth type!