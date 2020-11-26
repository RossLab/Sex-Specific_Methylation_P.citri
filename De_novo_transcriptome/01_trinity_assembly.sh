#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/XXX_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

# NOTE: conda activate afilia_trinity
#---------------------------------------------

echo "copying data in"

#---------------------------------------------

echo "doing the shiz"

# TRIM READS IF NEEDED BEFORE THIS

Trinity --seqType fq \ 
--left <read1> \
--right <read2> \
--SS_lib_type RF --max_memory 100G --CPU 32 --full_cleanup

#---------------------------------------------

echo "moving outputs"
mv ./*.txt /data/ross/mealybugs/analyses/hollie/sex-specific

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
