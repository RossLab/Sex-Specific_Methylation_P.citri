## -------------------------------------------------------------------------
# Filtering diff meth CpGs from methylkit with weighted meth of features
## -------------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation")
library(sqldf)
library(readr)
library(doBy)
library(ggplot2)
library(dplyr)

annotation <- read_delim("PCITRI_weighted_meth_annotation_common_both_sexes.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

diff_meth_sites <- read_csv("differential_meth_methylkit/F_vs_M__DMRs_min15percentDiff_qval0.01_MSCfilter.csv")
diff_meth_sites <- diff_meth_sites[,-c(1,2)]
colnames(diff_meth_sites)[7] <- "meth_diff"
colnames(diff_meth_sites)[2] <- "cpg_position"

## -------------------------------------------------------------------------

output <- sqldf("SELECT sample.chr,
                    sample.cpg_position,
                    sample.meth_diff,
                    sample.qvalue,
                    annot.scaffold,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.cpg_count,
                    annot.feature,
                    annot.male_mean_weightedMeth,
                    annot.female_mean_weightedMeth,
                    annot.methylated
                    FROM diff_meth_sites AS sample
                    LEFT JOIN annotation AS annot
                    ON sample.chr = annot.scaffold
                    AND (sample.cpg_position >= annot.start AND sample.cpg_position <= annot.end)")
output <- output[,-1]
output <- output[!(output$feature=="promotors_1000bp"),]

## -------------------------------------------------------------------------
# Where are these CpGs
## -------------------------------------------------------------------------

not_in_feature <- output[is.na(output$gene_id),] #21,820 (12%) Hmmm due to poor annotation???

output$feature[is.na(output$feature)] <- "Not_annotated" 
output <- subset(output, !feature == "whole_gene")

output$feature <- as.factor(output$feature)

ggplot(output, aes(x=feature))+
  geom_bar()+
  guides(fill=FALSE)+
  xlab("Feature")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  scale_x_discrete(breaks = c("exon_first3", "exon_notFirst3", "intron","Not_annotated","promotors_2000bp","TE"),
                   labels = c("Exons 1-3", "Exons 4+", "Introns", "No Annotation", "Promotors","TEs"),
                   limits = c("TE","exon_first3","promotors_2000bp", "Not_annotated", "intron","exon_notFirst3"))


## -------------------------------------------------------------------------
# Which genes are the promotor meth and which are the exons 1-3 meth
## -------------------------------------------------------------------------

promotor_genes <- output[output$feature == "promotors_2000bp",]
length(unique(promotor_genes$gene_id)) #8400 genes with diff cpg in promotor


exon_genes <- output[output$feature == "exon_first3",]
exon_genes$gene_id <- gsub("^.*Parent=", "", exon_genes$gene_id)
exon_genes$gene_id <- gsub("\\..*", "", exon_genes$gene_id)
length(unique(exon_genes$gene_id)) #6177 genes with diff cpg in exon


promotor_genes_only <- unique(promotor_genes$gene_id)
exon_genes_only <- unique(exon_genes$gene_id)

both <- Reduce(intersect, list(promotor_genes_only,exon_genes_only)) #3075 genes with both

## -------------------------------------------------------------------------
# How many CpGs per promotor are diff meth 
## -------------------------------------------------------------------------

library(data.table)

# Here we are saying a min of 3 diff meth cpgs per gene to count
number_diff_cpgs_per_promotor <- count(promotor_genes, gene_id)
median(number_diff_cpgs_per_promotor$n) #1-142, mean = 7.48, median = 2
hist(number_diff_cpgs_per_promotor$n)
nrow(number_diff_cpgs_per_promotor[number_diff_cpgs_per_promotor$n > 3,]) #2776 DROPS BIG TIME
promotors_with_3_cpgs <- subset(number_diff_cpgs_per_promotor, n >3)
ggplot(number_diff_cpgs_per_promotor, aes(x=n)) + 
  geom_histogram(color="black", fill="white", binwidth=3)+#bin every 3 (so second column = 3cpgs)
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  xlab("Number of Significant CpGs per Promotor")+
  ylab("Number of Promotors")

number_diff_cpgs_per_exon<- count(exon_genes, gene_id)
median(number_diff_cpgs_per_exon$n) # 1-230, mean=12.60, median=3
hist(number_diff_cpgs_per_exon$n)
nrow(number_diff_cpgs_per_exon[number_diff_cpgs_per_exon$n > 3,]) #2727 DOESN'T DROP AS MUCH
exons_with_3_cpgs <- subset(number_diff_cpgs_per_exon, n >3)
ggplot(number_diff_cpgs_per_exon, aes(x=n)) + 
  geom_histogram(color="black", fill="white", binwidth=3)+#bin every 3 (so second column = 3cpgs)
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  xlab("Number of Significant CpGs per Exon")+
  ylab("Number of Exons")

## -------------------------------------------------------------------------
# After have lists, filter on weighted meth difference of annotation
## -------------------------------------------------------------------------

output_filter1 <- output

# Add column for weighted methylation difference between males and females as a %
# Here a negative value means the male has less methylation
output_filter1$percent_meth_difference_of_feature <- ((output_filter1$male_mean_weightedMeth -
  output_filter1$female_mean_weightedMeth) / output_filter1$male_mean_weightedMeth)*100

# Most of the annotations by this point already have >15% diff
output_filter2 <- output_filter1[(output_filter1$percent_meth_difference_of_feature > 15 |
                                   output_filter1$percent_meth_difference_of_feature < -15),]

# Lets see how many genes now are left
promotor_genes <- subset(output_filter2, feature =="promotors_2000bp")
length(unique(promotor_genes$gene_id)) #4132 genes with diff cpg in promotor and 15% diff

exon_genes <- subset(output_filter2, feature =="exon_first3")
exon_genes$gene_id <- gsub("^.*Parent=", "", exon_genes$gene_id)
exon_genes$gene_id <- gsub("\\..*", "", exon_genes$gene_id)
length(unique(exon_genes$gene_id)) #4118 genes with diff cpg in exon and 15% diff

# What if we just keep genes with 3+ diff meth CpGs from methylkit
promotor_genes <- promotor_genes[promotor_genes$gene_id %in% promotors_with_3_cpgs$gene_id,]
length(unique(promotor_genes$gene_id)) #2709
exon_genes <- exon_genes[exon_genes$gene_id %in% exons_with_3_cpgs$gene_id,]
length(unique(exon_genes$gene_id)) #2736

promotor_genes_only <- unique(promotor_genes$gene_id)
exon_genes_only <- unique(exon_genes$gene_id)

both <- Reduce(intersect, list(promotor_genes_only,exon_genes_only)) #1522 genes with both

## -------------------------------------------------------------------------
# Before we get into genes for promotors and exons lets just look at TEs
## -------------------------------------------------------------------------
TEs <- subset(output_filter2, feature =="TE")
head(TEs)
TEs$uniquie_identified <- paste0(TEs$start, TEs$gene_id)
length(unique(TEs$uniquie_identified)) #28,056 TEs!

unique_TEs_only <- unique(TEs$uniquie_identified)
number_diff_cpgs_per_TE<- count(TEs,uniquie_identified )
TEs_with_3_cpgs <- subset(number_diff_cpgs_per_TE, n >3) #4605
# How many diff CpGs per TE annotation, distribution
hist(number_diff_cpgs_per_TE$n)
range(number_diff_cpgs_per_TE$n) # 1-257 diff CpGs per TE

ggplot(number_diff_cpgs_per_TE, aes(x=n)) + 
  geom_histogram(color="black", fill="white", binwidth=3)+#bin every 3 (so second column = 3cpgs)
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  xlab("Number of Significant CpGs per TE")+
  ylab("Number of TEs")
  

TEs <- TEs[TEs$uniquie_identified %in% TEs_with_3_cpgs$uniquie_identified,]
length(unique(TEs$uniquie_identified)) #4605

hyper_in_female_TEs<- subset(TEs, percent_meth_difference_of_feature < 0)
length(unique(hyper_in_female_TEs$uniquie_identified)) #4458/4605

hyper_in_male_TEs <- subset(TEs, percent_meth_difference_of_feature > 0)
length(unique(hyper_in_male_TEs$uniquie_identified)) #147/4605
# Hyper in both TEs? None


par(mfrow=c(2,2))
boxplot(hyper_in_female_TEs$male_mean_weightedMeth, main="female hyper TEs", 
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_female_TEs$female_mean_weightedMeth, main="female hyper TEs", 
        xlab="female", ylim=c(0,1))
boxplot(hyper_in_male_TEs$male_mean_weightedMeth, main="male hyper TEs",
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_male_TEs$female_mean_weightedMeth,  main="male hyper TEs",
        xlab="female", ylim=c(0,1))


## -------------------------------------------------------------------------
# Are these GENES hypermethylated in females or males, any with both? for promotos and exons
## -------------------------------------------------------------------------

head(promotor_genes)

hyper_in_female_promotor_genes <- subset(promotor_genes, percent_meth_difference_of_feature < 0)
length(unique(hyper_in_female_promotor_genes$gene_id)) #2645/2709

hyper_in_male_promotor_genes <- subset(promotor_genes, percent_meth_difference_of_feature > 0)
length(unique(hyper_in_male_promotor_genes$gene_id)) #64/2709

# Hyper in both promtors? None

par(mfrow=c(2,2))
boxplot(hyper_in_female_promotor_genes$male_mean_weightedMeth, main="female hyper promotors", 
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_female_promotor_genes$female_mean_weightedMeth, main="female hyper promotors", 
        xlab="female", ylim=c(0,1))
boxplot(hyper_in_male_promotor_genes$male_mean_weightedMeth, main="male hyper promotors",
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_male_promotor_genes$female_mean_weightedMeth,  main="male hyper promotors",
        xlab="female", ylim=c(0,1))



head(exon_genes)

hyper_in_female_exon_genes <- subset(exon_genes, percent_meth_difference_of_feature < 0)
length(unique(hyper_in_female_exon_genes$gene_id)) #2709/2736

hyper_in_male_exon_genes <- subset(exon_genes, percent_meth_difference_of_feature > 0)
length(unique(hyper_in_male_exon_genes$gene_id)) #33/2736

# Hyper in both sexes in exons? 6 genes!



boxplot(hyper_in_female_exon_genes$male_mean_weightedMeth, main="female hyper exons", 
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_female_exon_genes$female_mean_weightedMeth, main="female hyper exons", 
        xlab="female", ylim=c(0,1))
boxplot(hyper_in_male_exon_genes$male_mean_weightedMeth, main="male hyper exons",
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_male_exon_genes$female_mean_weightedMeth,  main="male hyper exons",
        xlab="female", ylim=c(0,1))


# Are genes with hypermeth in promotors the same as with hyper meth in exons?
hyper_female_both <- merge(hyper_in_female_promotor_genes, hyper_in_female_exon_genes, by = "gene_id")
length(unique(hyper_female_both$gene_id)) #1521

hyper_male_both <- merge(hyper_in_male_promotor_genes, hyper_in_male_exon_genes, by = "gene_id")
length(unique(hyper_male_both$gene_id)) #0


hyper_female_exon_hyper_male_promotor <-  merge(hyper_in_female_exon_genes, hyper_in_male_promotor_genes, by = "gene_id")
length(unique(hyper_female_exon_hyper_male_promotor$gene_id)) #0

hyper_male_exon_hyper_female_promotor <-  merge(hyper_in_male_exon_genes, hyper_in_female_promotor_genes, by = "gene_id")
length(unique(hyper_male_exon_hyper_female_promotor$gene_id)) #4

# Make a list of genes found in multiple datasets so can filter and get the unique ones
hyper_female_both_genes <- unique(hyper_female_both$gene_id)
hyper_male_both_genes <- unique(hyper_male_both$gene_id)
hyper_female_exon_hyper_male_promotor_genes <- unique(hyper_female_exon_hyper_male_promotor$gene_id)
hyper_male_exon_hyper_female_promotor_genes <- unique(hyper_male_exon_hyper_female_promotor$gene_id)

all_common_genes <- c(hyper_female_both_genes,hyper_male_both_genes,hyper_female_exon_hyper_male_promotor_genes,
                          hyper_male_exon_hyper_female_promotor_genes) #1525 correct


# Which genes then are unique to either exon or promotor
unique_female_hyper_promotor <- subset(hyper_in_female_promotor_genes, !(gene_id %in% all_common_genes))
length(unique(unique_female_hyper_promotor$gene_id)) #1123

unique_male_hyper_promotor <- subset(hyper_in_male_promotor_genes, !(gene_id %in% all_common_genes))
length(unique(unique_male_hyper_promotor$gene_id)) #64

unique_female_hyper_exon <- subset(hyper_in_female_exon_genes, !(gene_id %in% all_common_genes))
length(unique(unique_female_hyper_exon$gene_id)) #1188

unique_male_hyper_exon <- subset(hyper_in_male_exon_genes, !(gene_id %in% all_common_genes))
length(unique(unique_male_hyper_exon$gene_id)) #29


## -------------------------------------------------------------------------
# Make an Upset plot to show these overlaps better 
## -------------------------------------------------------------------------

library(UpSetR)

main_data <- as.data.frame(unique(c(promotor_genes_only, exon_genes_only)))
colnames(main_data) <- "gene_id"

main_data$`Hypermethylated in Female Promotor` <- 0
main_data$`Hypermethylated in Female Promotor`[main_data$gene_id %in% hyper_in_female_promotor_genes$gene_id] <- 1
main_data$`Hypermethylated in Male Promotor` <- 0
main_data$`Hypermethylated in Male Promotor`[main_data$gene_id %in% hyper_in_male_promotor_genes$gene_id] <- 1

main_data$`Hypermethylated in Female Exon` <- 0
main_data$`Hypermethylated in Female Exon`[main_data$gene_id %in% hyper_in_female_exon_genes$gene_id] <- 1
main_data$`Hypermethylated in Male Exon` <- 0
main_data$`Hypermethylated in Male Exon`[main_data$gene_id %in% hyper_in_male_exon_genes$gene_id] <- 1

upset(main_data, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity")



