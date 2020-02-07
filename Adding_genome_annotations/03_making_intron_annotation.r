#------------------------------------------------------------------------------------------------
### Take the intron file from script 02_making_intron_annotation.sh and add gene ID
#------------------------------------------------------------------------------------------------

# See https://wiki.lamp.le.ac.uk/mallonlab/index.php/Define_Intronic_Regions_in_GFF 

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files")
library(readr)
library(dplyr)
library(sqldf)

genes <-read.csv.sql("genes_with_start_and_end.txt",
                                      sql = "select * from file", sep = "\t", header = TRUE) 

intron_annot<-read.csv.sql("intron.bed",
                                     sql = "select * from file", sep = "\t", header = FALSE)
colnames(intron_annot) <- c("scaffold", "start", "end")

# Add intron information
output <- sqldf("SELECT intron.scaffold,
                intron.start,
                intron.end,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.strand,
                genes.gene_id
                FROM intron_annot AS intron
                LEFT JOIN genes AS genes
                ON intron.scaffold = genes.scaffold
                AND (intron.start >= genes.start AND intron.start <= genes.end)")

output<-subset(output, !output$gene_id=="NA")
output<-output[!duplicated(output),]
output <- output[,-c(4,5,6)]

output$feature<-"intron"
write.table(output, file="PCITRI_v0_introns.csv", 
            row.names = F, col.names = T, sep = '\t', quote = F)

