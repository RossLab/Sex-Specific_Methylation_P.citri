#GFF file doesn't contain intron regions so need to define them and add to a bed file with all annotation so can plot methylation one features.

#Found link which explains how to define introns from a gff:https://davetang.org/muse/2013/01/18/defining-genomic-regions/

#In brief:
#Use subtractBed command (part of bedtools) to subtract the exon regions from the overall gene region, leaving introns behind.

#First define all exon without overlap then take from gene:
cat PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3 |
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' |
sortBed |
mergeBed -i - | gzip > exon_merged.bed.gz

gunzip exon_merged.bed.gz

cat PCITRI.assembly.v0.braker.planococcus_citri.gt.gff3 |
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
sortBed |
subtractBed -a stdin -b exon_merged.bed > intron.bed

#Then can put this file into R along with the main annotation file and use SQL to add in intron info with gene name.