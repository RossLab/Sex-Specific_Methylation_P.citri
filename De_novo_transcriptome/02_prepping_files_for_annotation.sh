#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/annotation_prep_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

# NOTE: conda activate holl_trinity

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome/trinity_annotations_sex_specific.fasta ./

#---------------------------------------------

echo "build sqlite database and prepare sequences"
Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

echo "predict ORFs and protein coding"
# With homology
TransDecoder.LongOrfs -t trinity_annotations_sex_specific.fasta
# Without homology
TransDecoder.Predict -t trinity_annotations_sex_specific.fasta


#---------------------------------------------

echo "moving outputs"
mv * /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome/trinotate_stuff

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"