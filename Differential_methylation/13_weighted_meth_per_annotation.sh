#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/weighted_meth_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/methylation//02_methylation_extraction/coverage/*.txt ./
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/methylation//02_methylation_extraction/coverage/weighted_meth_per_sample.R ./

echo "doing the shiz"
R --save -q -f  weighted_meth_per_sample.R

echo "moving outputs"
mv ./*.txt /data/ross/mealybugs/analyses/hollie/sex-specific/methylation

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
