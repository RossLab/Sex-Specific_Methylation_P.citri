############################################################################
### Making genome information files
############################################################################

# Use one genome-wide cytosine report from methylation_extraction from bismark
# to create a file which consists of 
# just the scaffold name and the CpG position in a text file
# Only take + stranf coorsinate otherwise get each CpG counted twice
grep "+" trim_F1_1_bismark_bt2_pe.deduplicated.CpG_report.txt > plus_only.txt
cut -f1,2 plus_only.txt > total_cpgs_in_genome.txt

#---------------------------------------------------

# Make gene name and start and end positions
grep "gene" PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3 > genes.txt
cut -f1,4,5,7,9 genes.txt > genes_next.txt
sed 's/ID=//g' genes_next.txt > new.txt
echo -e "scaffold\tstart\tend\tstrand\tgene_id" | cat - new.txt > genes_with_start_and_end.txt

#---------------------------------------------------
# From coverage files need files to look like:
# chr, position, total coverage, count cystosines
gunzip *cov.gz

for file in $(ls *cov)
do
    base=$(basename ${file} "_1_bismark_bt2_pe.deduplicated.bismark.cov")
    cut -f1,2,5,6 ${file} > ${base}_coverage.txt
done







