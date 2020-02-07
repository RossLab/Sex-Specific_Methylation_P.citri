#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/methylkit_inputs_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/alignment/sorted_final_bams/*bam* ./
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/scripts/making_files.R ./

echo "doing the shiz"
R --save -q -f  making_files.R

echo "moving outputs"
mv ./*.txt /data/ross/mealybugs/analyses/hollie/sex-specific

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
