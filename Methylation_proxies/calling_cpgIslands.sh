#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/cpg_islands_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/genome/PCITRI.assembly.v0.fa ./

echo "doing the shiz"
/data/ross/mealybugs/analyses/hollie/bin/cpgiscan --length 200 --gap 100 \
-G cpg_island_annotation.gff PCITRI.assembly.v0.fa


echo "moving outputs"
mv ./*gff /data/ross/mealybugs/analyses/hollie/genome

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
