#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/sorting_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/alignment/*.bam ./

echo "doing the shiz"

for file in $(ls *.bam)
do
    base=$(basename $file "_1_bismark_bt2_pe.deduplicated.bam")
    samtools sort -o ${base}_deduplicated_sorted.bam -@ 10 ${file}
done

wait

for file in $(ls *sorted.bam)
do
    samtools index ${file}
done

echo "moving outputs"
mv ./*sorted* /data/ross/mealybugs/analyses/hollie/sex-specific

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
