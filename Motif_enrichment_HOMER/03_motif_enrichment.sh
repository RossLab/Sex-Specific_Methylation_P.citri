# -------------------------------------------------------------------
# HOMER motif enrichment analysis: methylated vs unmethylated
# -------------------------------------------------------------------

#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/motif_enrichment_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/motif_analysis/02_fasta_inputs/*.fa ./

echo "finding enriched motifs"
# p is for cores, b is for background set, -len is length of motif to look for (ran for standard 10 and also 5)

homer2 denovo -len 5 -p 8 -i all_exons_1-3.fa -b all_exons_3plus.fa > all_exon_1-3_vs_all_exon3plus_enriched_motifs.txt

homer2 denovo -len 5 -p 8 -i methylated_exons_1-3.fa -b unmethylated_exons_1-3.fa > meth_exons1-3_vs_unmeth_exons1-3_enriched_motifs.txt
homer2 denovo -len 5 -p 8 -i methylated_exons_3plus.fa -b unmethylated_exons_3plus.fa > meth_exons3plus_vs_unmeth_exons3plus_enriched_motifs.txt
homer2 denovo -len 5 -p 8 -i methylated_promotors.fa -b unmethylated_promotors.fa > meth_promotors_vs_unmeth_promotors_enriched_motifs.txt

echo "moving outputs"
mv ./*txt /data/ross/mealybugs/analyses/hollie/sex-specific/motif_analysis

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"


# -------------------------------------------------------------------
# Making HTML results, couldn't include in above
# -------------------------------------------------------------------

# HOMER doesn't like writing to a proper output, of course, so need to run each of 
# the below and then rename the output directory and .html that is produced, urgh
compareMotifs.pl all_exon_1-3_vs_all_exon3plus_enriched_motifs.txt ./ -pvalue 0.05

compareMotifs.pl meth_exons1-3_vs_unmeth_exons1-3_enriched_motifs.txt ./ -pvalue 0.05
compareMotifs.pl meth_exons3plus_vs_unmeth_exons3plus_enriched_motifs.txt ./ -pvalue 0.05
compareMotifs.pl meth_promotors_vs_unmeth_promotors_enriched_motifs.txt ./ -pvalue 0.05