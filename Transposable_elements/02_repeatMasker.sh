# ----------------------------------------------------------------
# RepeatMasker to annotate de novo TEs in the mealybug genome
# ----------------------------------------------------------------

# Use output from repeatModeler <file>.classified

# ----------------------------------------------------------------
# Script:

#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/repear_masker_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`


echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/genome/PCITRI.assembly.v0.fa ./
rsync /data/ross/mealybugs/analyses/hollie/genome/consensi.fa.classified ./

echo "doing the shiz"
# NOTE: repeat masker likes it on one line, due to extra spaces introduced across multiple lines
/data/ross/mealybugs/analyses/hollie/bin/RepeatMasker/RepeatMasker -pa 32 -norna -cutoff 250 -no_is -lib consensi.fa.classified -x -gff -a PCITRI.assembly.v0.fa

#-norna \ #doesn't mask small RNA / pseudo genes
#-cutoff 250 \ #cutoff when using lib (value from Sam Lewis Github)
#-no_is \ #skips bacterial insertion element check
#-lib consensi.fa.classified \
#-x \ #returns repetitive regions marked with X's not N's
#-gff \ #creates gff output
#-a \ #writes alignments out 

echo "moving outputs"
cd ..
mv ./* /data/ross/mealybugs/analyses/hollie/genome/repeatMasker

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*


end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"