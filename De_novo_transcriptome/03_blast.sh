#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/blast_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

# NOTE: conda activate holl_trinity

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome/trinity_annotations_sex_specific.fasta ./
rsync -r /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome/trinotate_stuff/ ./
#---------------------------------------------

echo "homology searching"
blastx -query trinity_annotations_sex_specific.fasta -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
blastp -query trinity_annotations_sex_specific.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
hmmscan --cpu 16 --domtblout TrinotatePFAM.out Pfam-A.hmm trinity_annotations_sex_specific.fasta.transdecoder.pep > pfam.log

#---------------------------------------------

echo "moving outputs"
mv *.outfmt6 /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome
mv *pfam.log /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome
echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"