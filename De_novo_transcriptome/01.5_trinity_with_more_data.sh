#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/trinity_with_extra_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

# NOTE: conda activate holl_trinity
#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/02_trimmed/*.fq.gz ./
rsync /data/ross/mealybugs/analyses/ase_in_mealybugs/intraspecific/1_reads/CW* ./
rsync /data/ross/mealybugs/analyses/ase_in_mealybugs/intraspecific/1_reads/WC* ./
#---------------------------------------------

echo "gettin them annotations"
Trinity --seqType fq \
--left trim_PC_F1_1.fq.gz,trim_PC_F2_1.fq.gz,trim_PC_F3_1.fq.gz,trim_PC_M1_1.fq.gz,trim_PC_M2_1.fq.gz,trim_PC_M3_1.fq.gz,CW10_1.trim.fastq.gz,CW11_1.trim.fastq.gz,CW1_1.trim.fastq.gz,CW4_1.trim.fastq.gz,WC1_1.trim.fastq.gz,WC19_1.trim.fastq.gz,WC4_1.trim.fastq.gz,WC6_1.trim.fastq.gz \
--right trim_PC_F1_2.fq.gz,trim_PC_F2_2.fq.gz,trim_PC_F3_2.fq.gz,trim_PC_M1_2.fq.gz,trim_PC_M2_2.fq.gz,trim_PC_M3_2.fq.gz,CW10_2.trim.fastq.gz,CW11_2.trim.fastq.gz,CW1_2.trim.fastq.gz,CW4_2.trim.fastq.gz,WC1_2.trim.fastq.gz,WC19_2.trim.fastq.gz,WC4_2.trim.fastq.gz,WC6_2.trim.fastq.gz \
--SS_lib_type RF --max_memory 150G --CPU 32 --full_cleanup

#---------------------------------------------

echo "moving outputs"
mv *.fasta /data/ross/mealybugs/analyses/hollie/sex-specific/expression/de_novo_transcriptome/with_extra_data

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
