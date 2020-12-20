#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/quantification_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------
# NOTE removed invorrect format lines from original pasa .gff3:
# grep -v "#" pcitri_v0_updated_with_pasa.gff3 > pcitri_v0_updated_with_pasa_no_hashes.gff3
# sed -i '/^[[:space:]]*$/d' pcitri_v0_updated_with_pasa_no_hashes.gff3

echo "copying data in"
rsync /data/ross/mealybugs/analyses/hollie/genome/citri_v0/pcitri_v0_updated_with_pasa_no_hashes.gff3 ./
rsync /data/ross/mealybugs/analyses/hollie/genome/citri_v0/PCITRI.assembly.v0.fa ./
#rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/03_bams/*bam ./
rsync /data/ross/mealybugs/analyses/hollie/sex-specific/expression/02_trimmed/*.gz ./
#---------------------------------------------
# https://github.com/deweylab/RSEM

echo "prepare rsem reference"
rsem-prepare-reference --gff3 pcitri_v0_updated_with_pasa_no_hashes.gff3 --star -p 32 PCITRI.assembly.v0.fa pcitri

#echo "make valid bams"
#for file in $(ls *.bam)
#do
#	base=$(basename ${file} "_Aligned.sortedByCoord.out.bam")
#    convert-sam-for-rsem -p 10 ${file} ${base}_converted.bam
#done

#echo "calc RSEM"
#for file in $(ls *_converted.bam)
#do
#	base=$(basename ${file} "_converted.bam")
#    rsem-calculate-expression --alignments --paired-end -p 32 ${file} pcitri ${base}
#done

echo "run RSEM fresh"
for file in $(ls *1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    rsem-calculate-expression --star -p 20 --paired-end ${base}_1.fq.gz ${base}_2.fq.gz pcitri ${base}
done

#---------------------------------------------

echo "moving outputs"
mv * /data/ross/mealybugs/analyses/hollie/sex-specific/expression/04_RSEM_new

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
