#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/alignment_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/trimmed_raw/*.fq.gz ./
rsync -r /data/ross/mealybugs/analyses/hollie/genome ./

echo "doing the shiz"
for file in $(ls *1.fq.gz)
do
	base=$(basename ${file} "1.fq.gz")
	bismark --multicore 8 -o ./alignment \
	--genome $SCRATCH/genome \
	-1 ${base}1.fq.gz \
	-2 ${base}2.fq.gz
done


echo "moving outputs"
mv ./alignment /data/ross/mealybugs/analyses/hollie/sex-specific

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
