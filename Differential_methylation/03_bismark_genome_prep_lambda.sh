#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/genome_prep_lambda_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
#cp /data/ross/mealybugs/analyses/hollie/sex-specific/trimmed_raw/*.fq.gz ./
cp /data/ross/mealybugs/analyses/hollie/lambda_genome/lambda.fa ./

echo "doing the shiz"
#for file in $(ls *1.fq.gz)
#do
#	base=$(basename $file "1.fq.gz")
#done

bismark_genome_preparation ./

echo "moving outputs"
mv * /data/ross/mealybugs/analyses/hollie/lambda_genome

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
