#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/meth_extraction_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/alignment/*.bam ./
rsync -r /data/ross/mealybugs/analyses/hollie/genome ./

echo "doing the shiz"
mkdir methylation_extraction

for file in $(ls *.bam)
do
    bismark_methylation_extractor -p \
    --comprehensive \
    --multicore 10 \
    --bedgraph \
    --cytosine_report \
    --genome_folder $SCRATCH/genome \
    --scaffolds \
    --output ./methylation_extraction \
    ${file}
done


echo "moving outputs"
mv ./methylation_extraction /data/ross/mealybugs/analyses/hollie/sex-specific

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
