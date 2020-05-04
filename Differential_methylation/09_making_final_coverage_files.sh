#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/making_cov_files_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/methylation_extraction/coverage/*.cov.gz ./
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/methylation_extraction/coverage/making_final_coverage_files.R ./

echo "doing the shiz"
gunzip *cov.gz

for file in $(ls *cov)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bismark.cov")
    cut -f1,2,5,6 ${file} > ${base}_coverage.txt
done

R --save -q -f  making_final_coverage_files.R

echo "moving outputs"
mv ./*final* /data/ross/mealybugs/analyses/hollie/sex-specific

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
