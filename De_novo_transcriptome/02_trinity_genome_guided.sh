#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/trinity_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

# NOTE: conda activate holl_trinity
#---------------------------------------------

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/03_bams/*.bam ./

#---------------------------------------------

echo "making one bam"
samtools merge -@ 10 all_bams.bam trim_PC_F1_Aligned.sortedByCoord.out.bam \
trim_PC_M1_Aligned.sortedByCoord.out.bam trim_PC_F2_Aligned.sortedByCoord.out.bam \
trim_PC_M2_Aligned.sortedByCoord.out.bam trim_PC_F3_Aligned.sortedByCoord.out.bam \
trim_PC_M3_Aligned.sortedByCoord.out.bam

echo "sorting bam"
samtools sort -@ 10 -o all_bams_sorted.bam all_bams.bam

echo "gettin them annotations"
 Trinity --genome_guided_bam all_bams_sorted.bam \
         --genome_guided_max_intron 10000 \
         --SS_lib_type RF --max_memory 100G --CPU 32 
#---------------------------------------------

echo "moving outputs"
mv *.fasta /data/ross/mealybugs/analyses/hollie/sex-specific/expression

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"