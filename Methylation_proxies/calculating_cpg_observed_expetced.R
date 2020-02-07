#---------------------------------------------
# Calculating CpG observed / expected 
#---------------------------------------------

library(readr)
library(seqinr)
library(stringr)
library(sqldf)
library(ggplot2)
library(ggpmisc)

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation")

#---------------------------------------------
# Get gene lengths for all genes
#---------------------------------------------
genes_with_start_and_end <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files/annotation_files/genes_with_start_and_end.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)

genes_with_start_and_end$length <- genes_with_start_and_end$end - genes_with_start_and_end$start
range(genes_with_start_and_end$length) # 72-207868

#---------------------------------------------
# Count number of C's and G's per gene
#---------------------------------------------
gene_sequence <- read.fasta("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files/PCITRI_genes.fa")

seq_names <- character(length = 39801)
c_count <- numeric(length = 39801)
g_count <- numeric(length = 39801)

for(i in seq_along(gene_sequence)){
  seq_names[i] <- names(gene_sequence)[i]
  count <- as.data.frame(count(gene_sequence[[i]],1))
  c_count[i] <- count[2,"Freq"]
  g_count[i] <- count[3,"Freq"]
}

counts_c_g_per_gene <- data.frame(seq_names,c_count,g_count)

#---------------------------------------------
# Match up gene length with counts
#---------------------------------------------

# Need to split the seq_names into columns for matching
head(counts_c_g_per_gene)
counts_c_g_per_gene$scaffold <- str_extract(counts_c_g_per_gene$seq_names, "[^:]+")
counts_c_g_per_gene$end <- str_extract(counts_c_g_per_gene$seq_names, "[^-]+$")
counts_c_g_per_gene$start <- str_extract(counts_c_g_per_gene$seq_names, "[^:]+$")
counts_c_g_per_gene$start <- gsub("-.*","",counts_c_g_per_gene$start)

# Match on scaffold and gene start and end to get gene id 
output <- sqldf("SELECT genes.scaffold,
                    genes.start,
                    genes.end,
                    genes.gene_id,
                    genes.length,
                    counts.scaffold,
                    counts.start,
                    counts.end,
                    counts.c_count,
                    counts.g_count
                    FROM genes_with_start_and_end AS genes
                    LEFT JOIN counts_c_g_per_gene AS counts
                    ON genes.scaffold = counts.scaffold
                    AND (genes.start = counts.start AND genes.end = counts.end)")
output <- output[,-c(6,7,8)]

#---------------------------------------------
# Add in the CpG count per gene
#---------------------------------------------

cpg_per_gene <- read_table2("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files/annotation_files_with_total_cpgs/PCITRI_v0_genes_with_total_cpgs.txt")

all_data <- merge(output, cpg_per_gene)

#---------------------------------------------
# Calculate the o/e for each gene
#---------------------------------------------
# CpG o/e = (length2/length)*(CpG count/(C count * G count))
# From Simola et al. (2013) doi:10.1101/gr.155408.113
# Althought this makes no sense as (length*length) / length is just length ???
# From Liu et al (2019) doi:10.3390/genes10020137  (L*#CpG)/(#C*#G)
# Other sources to normalise by length you do: length^2 / length -1
# Just gives the same as Liu et al. so will stick with that

all_data$cpg_ob_ex <- all_data$length * (all_data$cpg_counter.sum/(all_data$c_count*all_data$g_count))

head(all_data)
range(all_data$cpg_ob_ex)
plot(density(all_data$cpg_ob_ex))

# Weird ... a few genes Cpg o/e > 1, up to 9 :S biological or mistake?
all_data_subset <- subset(all_data, cpg_ob_ex < 2)
plot(density(all_data_subset$cpg_ob_ex)) #39228 / 39752 have cpg o/e < 1.5

#---------------------------------------------
# Make a fancy CpG o/e plot
#---------------------------------------------

# Find the value of the biggest peak
density(all_data_subset$cpg_ob_ex)
which.max(density(all_data_subset$cpg_ob_ex)$y) # top peak 281
density(all_data_subset$cpg_ob_ex)$x[281] # top peak to plot on X = 1.113309

# Find the trough so we can find the second peak
DensityY <- density(all_data_subset$cpg_ob_ex)$y
DensityX <- density(all_data_subset$cpg_ob_ex)$x
# so only looking between the two peaks to avoid the tails
MinYDensity<- min(DensityY[DensityX < 1.1 & DensityX > 0.5]) 
which(DensityY == MinYDensity) #149
DensityX[149] #0.5476148

# Find the value of the second biggest peak
MaxY<- max(density(all_data_subset$cpg_ob_ex)$y[density(all_data_subset$cpg_ob_ex)$x < 0.5476148])
which(density(all_data_subset$cpg_ob_ex)$y == MaxY) #138
density(all_data_subset$cpg_ob_ex)$x[138]#0.5004735

ggplot(all_data_subset, aes(x=cpg_ob_ex))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2, colour = "black", bins = 40 )+
  geom_line(stat="density",size=1.5, colour="black")+
  geom_vline(aes(xintercept=0.5004735),linetype="dashed", size=1.5,colour="black")+
  geom_vline(aes(xintercept=1.113309),linetype="dashed",size=1.5,colour="black")+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  xlab("CpG Observed/Expected")+
  ylab("Density")
  #annotate("text", x = 1.75, y = 1.2, label = "P. citri", size=15)


#-----------------------------------------------------------------------------------------------------

# I wonder if methylated promtors (according to BS data) have a noticable difference in CpG o/e value

#---------------------------------------------
# Good check for the proxy holding true
#---------------------------------------------

proms_with_total_cpgs <- read_table2("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files/annotation_files_with_total_cpgs/PCITRI_v0_1000bp_promotors_with_total_cpgs.txt")

# ISSUE: these .fa are for 2000bp promotors not 1000bp promotors urgh
meth_proms_fasta <- read.fasta("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/methylated_promotors.fa")
unmeth_proms_fasta <- read.fasta("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/unmethylated_promotors.fa")

#---------------------------------------------
# For meth proms
seq_names <- character(length = 6750)
c_count <- numeric(length = 6750)
g_count <- numeric(length = 6750)

for(i in seq_along(meth_proms_fasta)){
  seq_names[i] <- names(meth_proms_fasta)[i]
  count <- as.data.frame(count(meth_proms_fasta[[i]],1))
  c_count[i] <- count[2,"Freq"]
  g_count[i] <- count[3,"Freq"]
}

counts_c_g_per_meth_prom <- data.frame(seq_names,c_count,g_count)


# For meth proms
seq_names <- character(length = 37425)
c_count <- numeric(length = 37425)
g_count <- numeric(length = 37425)

for(i in seq_along(unmeth_proms_fasta)){
  seq_names[i] <- names(unmeth_proms_fasta)[i]
  count <- as.data.frame(count(unmeth_proms_fasta[[i]],1))
  c_count[i] <- count[2,"Freq"]
  g_count[i] <- count[3,"Freq"]
}

counts_c_g_per_unmeth_prom <- data.frame(seq_names,c_count,g_count)

#---------------------------------------------
# Split seq name
head(counts_c_g_per_unmeth_prom)
counts_c_g_per_unmeth_prom$scaffold <- str_extract(counts_c_g_per_unmeth_prom$seq_names, "[^:]+")
counts_c_g_per_unmeth_prom$end <- str_extract(counts_c_g_per_unmeth_prom$seq_names, "[^-]+$")
counts_c_g_per_unmeth_prom$start <- str_extract(counts_c_g_per_unmeth_prom$seq_names, "[^:]+$")
counts_c_g_per_unmeth_prom$start <- gsub("-.*","",counts_c_g_per_unmeth_prom$start)

head(counts_c_g_per_meth_prom)
counts_c_g_per_meth_prom$scaffold <- str_extract(counts_c_g_per_meth_prom$seq_names, "[^:]+")
counts_c_g_per_meth_prom$end <- str_extract(counts_c_g_per_meth_prom$seq_names, "[^-]+$")
counts_c_g_per_meth_prom$start <- str_extract(counts_c_g_per_meth_prom$seq_names, "[^:]+$")
counts_c_g_per_meth_prom$start <- gsub("-.*","",counts_c_g_per_meth_prom$start)

#---------------------------------------------
# annotate with gene id and number cpgs

head(proms_with_total_cpgs)
head(counts_c_g_per_meth_prom)

output_meth <- merge(proms_with_total_cpgs, counts_c_g_per_meth_prom)
output_unmeth <- merge(proms_with_total_cpgs, counts_c_g_per_unmeth_prom)

output_meth$length <- 1000
output_unmeth$length <- 1000

#---------------------------------------------
# calculate cpg o/e

output_meth$cpg_ob_ex <- output_meth$length * (output_meth$cpg_counter.sum/(output_meth$c_count*output_meth$g_count))
head(output_meth)
range(output_meth$cpg_ob_ex) #0.03881988 111.11111111
plot(density(output_meth$cpg_ob_ex))

output_unmeth$cpg_ob_ex <- output_unmeth$length * (output_unmeth$cpg_counter.sum/(output_unmeth$c_count*output_unmeth$g_count))
head(output_unmeth)
range(output_unmeth$cpg_ob_ex) #0.03496503 55.55555556
plot(density(output_unmeth$cpg_ob_ex))

#---------------------------------------------
# make one data frame for plotting

output_meth$prom_status <- "methylated"
output_unmeth$prom_status <-"unmethylated"

output_meth_subset <- subset(output_meth, cpg_ob_ex < 2) #6612/6750
output_unmeth_subset <- subset(output_unmeth, cpg_ob_ex < 2) #22857/23376

all_prom_data <- rbind(output_meth_subset, output_unmeth_subset)

#---------------------------------------------
# make plot

# Find the value of the biggest peak for methylated
plot(density(output_meth_subset$cpg_ob_ex))
which.max(density(output_meth_subset$cpg_ob_ex)$y) # top peak 222
density(output_meth_subset$cpg_ob_ex)$x[222] # top peak to plot on X = 0.8621213

# Find the value of the biggest peak for unmethylated
plot(density(output_unmeth_subset$cpg_ob_ex))
which.max(density(output_unmeth_subset$cpg_ob_ex)$y) # top peak 274
density(output_unmeth_subset$cpg_ob_ex)$x[274] # top peak to plot on X = 1.09268


ggplot(all_prom_data, aes(x=cpg_ob_ex, colour=prom_status))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2, bins = 40 )+
  geom_line(stat="density",size=1.5)+
  geom_vline(aes(xintercept=0.8621213),linetype="dashed", size=1.5,colour="black")+
  geom_vline(aes(xintercept=1.09268),linetype="dashed",size=1.5,colour="black")+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        title = element_text(size=20))+
  xlab("CpG Observed/Expected")+
  ylab("Density")+
  ggtitle("Promotors")


ggplot(all_prom_data, aes(x=cpg_ob_ex))+
  geom_histogram(aes(y=..density..), position="identity", alpha=0.2, bins = 40,colour = "black" )+
  geom_line(stat="density",size=1.5)+
  geom_vline(aes(xintercept=0.8621213),linetype="dashed", size=1.5,colour="black")+
  geom_vline(aes(xintercept=1.09268),linetype="dashed",size=1.5,colour="black")+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        title = element_text(size=20))+
  xlab("CpG Observed/Expected")+
  ylab("Density")+
  ggtitle("Promotors")
