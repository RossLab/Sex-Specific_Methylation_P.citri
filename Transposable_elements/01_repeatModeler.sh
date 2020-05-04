# ----------------------------------------------------------------
# Running RepeatModeler to find de novo TEs in the mealybug genome
# ----------------------------------------------------------------

# Following protocol of Lewis et al. (2019) pre-print from Laura's lab
# Downloaded and configured RepeatModeler and RepeatMasker along with
# all dependencies according to their respective websites

# Make a database of the fasta file you want to annotate TEs in 
/data/ross/mealybugs/analyses/hollie/bin/RepeatModeler-2.0/BuildDatabase \
--name database_PCITRI_repeatmodeler PCITRI.assembly.v0.fa

# The RMBlast uses 4 cores at a time, so specify cores in multiples of 4
/data/ross/mealybugs/analyses/hollie/bin/RepeatModeler-2.0/RepeatModeler \
-database database_PCITRI_repeatmodeler -pa 32

# ----------------------------------------------------------------
# Script:

#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/repeat_modeler_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`


echo "copying data in"
cp -r /data/ross/mealybugs/analyses/hollie/genome ./

echo "doing the shiz"
cd genome
/data/ross/mealybugs/analyses/hollie/bin/RepeatModeler-2.0/RepeatModeler \
-database database_PCITRI_repeatmodeler -pa 32

echo "moving outputs"
cd ..
mv ./genome /data/ross/mealybugs/analyses/hollie/genome
echo "a clean directory is a happy directory"
rm -r $SCRATCH/*


end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"


# Takes a couple of days to run, the output file is interest is within a subfolder created
# by repeatModeler and ends with <file>.classififed (this is the library needed for repeatMasker)