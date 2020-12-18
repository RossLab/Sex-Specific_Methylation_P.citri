#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/blast_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------
echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/genome/citri_v0/Planococcus_citri_genes.fa ./
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome/trinity_annotations_sex_specific.fasta ./

echo "make blast database"
/ceph/software/blast/ncbi-blast-2.8.0+/bin/makeblastdb \
-in Planococcus_citri_genes.fa \
-dbtype nucl \
-parse_seqids \
-out P_citri_genome_db

echo "blast trinity transcripts to genome annotations"
blastn -query trinity_annotations_sex_specific.fasta \
-db P_citri_genome_db \
-num_threads 16 \
-max_target_seqs 1 \
-outfmt 6 \
-evalue 1e-3 > trinity_to_genome_blast.txt

#---------------------------------------------

echo "moving outputs"
mv ./trinity_to_genome_blast.txt /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
