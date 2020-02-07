#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/deduplication_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/alignment/original_bams/*.bam ./

echo "doing the shiz"
for file in $(ls *bam)
do
	deduplicate_bismark -p --bam ${file}
done


echo "moving outputs"
mv ./*dedup* /data/ross/mealybugs/analyses/hollie/sex-specific/alignment

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
