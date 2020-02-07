#$ -V
#$ -cwd
#$ -j y
#$ -o /ceph/users/hmarshall/logs/trimming_bisulfite_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
cp /data/ross/mealybugs/raw/bgi_bisulfite/0_reads/*.fq.gz ./

echo "doing the shiz"
for file in $(ls *1.fq.gz)
do
	base=$(basename $file "1.fq.gz")
	cutadapt -j 0 \
	-u 4 -U 4 \
	-o trim_${base}1.fq.gz \
	-p trim_${base}2.fq.gz \
	${base}1.fq.gz \
	${base}2.fq.gz
done

echo "moving outputs"
mv trim_* /data/ross/mealybugs/analyses/hollie/sex-specific/trimmed_raw

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
