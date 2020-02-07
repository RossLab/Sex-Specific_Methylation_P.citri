#$ -V
#$ -cwd
#$ -j y
#$ -o /ceph/users/hmarshall/logs/fastqc_bisulfite_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
cp /data/ross/mealybugs/raw/bgi_bisulfite/0_reads/*.fq.gz ./

echo "doing the shiz"
for file in $(ls *.gz)
do
	fastqc -t 9 ${file}
done

echo "moving outputs"
mv *html /data/ross/mealybugs/analyses/hollie/sex-specific/01_fastqc_bisulfite

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
