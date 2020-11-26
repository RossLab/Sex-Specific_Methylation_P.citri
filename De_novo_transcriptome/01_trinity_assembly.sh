#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/trinity_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

# NOTE: conda activate holl_trinity
#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/02_trimmed/*.fq.gz ./

#---------------------------------------------

echo "gettin them annotations"
Trinity --seqType fq \
--left trim_PC_F1_1.fq.gz,trim_PC_F2_1.fq.gz,trim_PC_F3_1.fq.gz,trim_PC_M1_1.fq.gz,trim_PC_M2_1.fq.gz,trim_PC_M3_1.fq.gz \
--right trim_PC_F1_2.fq.gz,trim_PC_F2_2.fq.gz,trim_PC_F3_2.fq.gz,trim_PC_M1_2.fq.gz,trim_PC_M2_2.fq.gz,trim_PC_M3_2.fq.gz \
--SS_lib_type RF --max_memory 100G --CPU 32 --full_cleanup

#---------------------------------------------

echo "moving outputs"
mv *.fasta /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
