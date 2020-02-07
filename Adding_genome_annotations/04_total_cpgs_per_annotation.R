## -------------------------------------------------------------------------
## Weighted Methylation per annotation: make file with cpg count per annoation
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files")
library(sqldf)
library(readr)
library(doBy)

## -------------------------------------------------------------------------

# Making the file which has the total CpGs per gene information

cpgs <- read.delim("total_cpgs_in_genome.txt", header=F)
colnames(cpgs) <- c("scaffold", "cpg_position")
cpgs$cpg_position <- as.numeric(cpgs$cpg_position)

# --------------------------------------------------------------------

# GENES

genes <- read.delim("genes_with_start_and_end.txt", header=T)
genes$start <- as.numeric(genes$start)
genes$end <- as.numeric(genes$end)

output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.gene_id
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$gene_id),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ gene_id+scaffold+start+end, data = output, FUN=sum)

write.table(final, file="PCITRI_v0_genes_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)

# --------------------------------------------------------------------

# EXONS (FIRST 3 ONLY)

genes <- read.delim("first_3_exons_annotation.txt", header=F)
colnames(genes) <- c("scaffold","annotation_source","feature","start","end","strand","info")

output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.info
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$info),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ info+scaffold+start+end, data = output, FUN=sum)

write.table(final, file="PCITRI_v0_first3exons_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)

# --------------------------------------------------------------------

# PROMOTORS (HOLL DEFINED 1000bp)

genes <- read.delim("PCITRI_v0_gene_promotors_1000bp.csv", header=T, sep="\t")

output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.gene_id
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$gene_id),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ gene_id+scaffold+start+end, data = output, FUN=sum)

write.table(final, file="PCITRI_v0_1000bp_promotors_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)

# --------------------------------------------------------------------

# PROMOTORS (STEVIE DEFINED 2000bp)

genes <- read.delim("PCITRI.v0.PROMOTORS_STEVIE.bed", header=F)
colnames(genes) <- c("scaffold","start","end","gene_id","something", "strand","info","feature","info2")


output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.gene_id
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$gene_id),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ gene_id+scaffold+start+end, data = output, FUN=sum)

write.table(final, file="PCITRI_v0_2000bp_promotors_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)

# --------------------------------------------------------------------

# INTRONS

genes <- read.delim("PCITRI_v0_introns.csv", header=T)

output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.gene_id
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$gene_id),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ gene_id+scaffold+start+end, data = output, FUN=sum)

write.table(final, file="PCITRI_v0_introns_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)


# --------------------------------------------------------------------

# EXONS (WITHOUT FIRST 3)

genes <- read.delim("exons_without_first3_annotation.txt", header=F)
colnames(genes) <- c("scaffold","something","feature","start","end", "strand","info")


output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.info
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$info),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ info+scaffold+start+end, data = output, FUN=sum)

write.table(final, file="PCITRI_v0_exons_withoutFirst3_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)

# --------------------------------------------------------------------

# TE's

genes <- read.delim("PCITRI_v0_TEs.bed", header=F)
colnames(genes) <- c("scaffold","annotation_source","something","start","end","something1","strand",
                     "something2","info")

output <- sqldf("SELECT cpgs.scaffold,
                cpgs.cpg_position,
                genes.scaffold,
                genes.start,
                genes.end,
                genes.info
                FROM cpgs AS cpgs
                LEFT JOIN genes AS genes 
                ON cpgs.scaffold = genes.scaffold
                AND (cpgs.cpg_position >= genes.start AND cpgs.cpg_position <= genes.end)")
output <- output[!is.na(output$info),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ info+scaffold+start+end, data = output, FUN=sum)

colnames(final) <- NULL
write.table(final, file="PCITRI_v0_TEs_with_total_cpgs.txt", col.names=F, row.names=F, quote=F)



# --------------------------------------------------------------------
# Stuff it, put what I currently have together

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files/annotation_files_with_total_cpgs")

promotors_1000bp <- read_table2("PCITRI_v0_1000bp_promotors_with_total_cpgs.txt")
promotors_2000bp <- read_table2("PCITRI_v0_2000bp_promotors_with_total_cpgs.txt")
exons_withoutFirst3 <- read_table2("PCITRI_v0_exons_withoutFirst3_with_total_cpgs.txt")
exons_first3 <- read_table2("PCITRI_v0_first3exons_with_total_cpgs.txt")
genes <- read_table2("PCITRI_v0_genes_with_total_cpgs.txt")
introns <- read_table2("PCITRI_v0_introns_with_total_cpgs.txt")
TEs <- read_table2("PCITRI_v0_TEs_with_total_cpgs.txt", col_names = F)

# Re-name exon infor columns to match others so can bind
colnames(exons_withoutFirst3)[1] <- "gene_id"
colnames(exons_first3)[1] <- "gene_id"

# Fix the TE file
TEs$info <- paste0(TEs$X1, TEs$X2, TEs$X3, TEs$X4)
TEs <- TEs[,c(5,6,7,8,9)]
colnames(TEs) <- c("scaffold","start","end","cpg_counter.sum","gene_id")

# Add identifying column
promotors_1000bp$feature <- "promotors_1000bp"
promotors_2000bp$feature <- "promotors_2000bp"
exons_withoutFirst3$feature <- "exon_notFirst3"
exons_first3$feature <- "exon_first3"
genes$feature <- "whole_gene"
introns$feature <- "intron"
TEs$feature <- "TE"

all_data <- rbind(promotors_1000bp, promotors_2000bp, exons_withoutFirst3, exons_first3, genes, introns, TEs)
colnames(all_data)[5] <- "cpg_count"

write.table(all_data, file="merged_annotations.txt", sep="\t", quote =F, col.names = T, row.names = F)
