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
  guides()+
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

promotor_genes$length <- promotor_genes$end - promotor_genes$start
promotor_genes_normalised <- merge(number_diff_cpgs_per_promotor, promotor_genes, by="gene_id")
promotor_genes_normalised$normalised_cpg_count <- promotor_genes_normalised$n / promotor_genes_normalised$length

ggplot(promotor_genes_normalised, aes(x=normalised_cpg_count)) + 
  geom_histogram(color="black", fill="white")+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  xlab("Number of Significant CpGs per Promotor")+
  ylab("Number of Promotors")+
  xlim(0,0.15)

number_diff_cpgs_per_exon<- count(exon_genes, gene_id)
median(number_diff_cpgs_per_exon$n) # 1-230, mean=12.60, median=3
hist(number_diff_cpgs_per_exon$n)
nrow(number_diff_cpgs_per_exon[number_diff_cpgs_per_exon$n > 3,]) #2761 DOESN'T DROP AS MUCH
exons_with_3_cpgs <- subset(number_diff_cpgs_per_exon, n >3)

exon_genes$length <- exon_genes$end - exon_genes$start
exon_genes_normalised <- merge(number_diff_cpgs_per_exon, exon_genes, by="gene_id")
exon_genes_normalised$normalised_cpg_count <- exon_genes_normalised$n / exon_genes_normalised$length

ggplot(exon_genes_normalised, aes(x=normalised_cpg_count)) + 
  geom_histogram(color="black", fill="white")+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  xlab("Number of Significant CpGs per Exon")+
  ylab("Number of Exons")+
  xlim(0,0.15)

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
#length(unique(promotor_genes$gene_id)) #7816 genes with diff cpg in promotor and 15% diff
exon_genes <- subset(output_filter2, feature =="exon_first3")
exon_genes$gene_id <- gsub("^.*Parent=", "", exon_genes$gene_id)
exon_genes$gene_id <- gsub("\\..*", "", exon_genes$gene_id)
#length(unique(exon_genes$gene_id)) #6019 genes with diff cpg in exon and 15% diff

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

# Filter out all TEs overlapping a significant promotor or exon
# First make a list of CpGs which are diff meth in TEs
cpgs_in_TEs <- unique(TEs$cpg_position)
cpgs_in_proms <- unique(output_filter2$cpg_position[output_filter2$feature=="promotors_2000bp"])
cpgs_in_exons <- unique(output_filter2$cpg_position[output_filter2$feature=="exon_first3"])

cpgs_in_both <- cpgs_in_exons[cpgs_in_exons %in% cpgs_in_proms] #25808 cpgs in both proms and exons ... urgh

cpgs_unique_to_TEs <- cpgs_in_TEs[!cpgs_in_TEs %in% cpgs_in_proms]
cpgs_unique_to_TEs <- cpgs_unique_to_TEs[!cpgs_in_TEs %in% cpgs_in_exons]

TEs_subset <- TEs[TEs$cpg_position %in% cpgs_unique_to_TEs,]

# Filter by 3cpgs per TE
unique_TEs_only <- unique(TEs_subset$uniquie_identified)
number_diff_cpgs_per_TE<- count(TEs_subset,uniquie_identified )
TEs_with_3_cpgs <- subset(number_diff_cpgs_per_TE, n >3) #1061

# How many diff CpGs per TE annotation, distribution
hist(number_diff_cpgs_per_TE$n)
range(number_diff_cpgs_per_TE$n) # 1-115 diff CpGs per TE

# Normalise by length (but what does this actually tell us ...)
# Trying to figure out if there are more CpGs diff meth is that becasue it's longer?
TEs$length <- TEs$end - TEs$start
number_diff_cpgs_per_TE_normalised <- merge(number_diff_cpgs_per_TE, TEs, by="uniquie_identified")
number_diff_cpgs_per_TE_normalised$normalised_cpg_count <- number_diff_cpgs_per_TE_normalised$n / number_diff_cpgs_per_TE_normalised$length

ggplot(number_diff_cpgs_per_TE_normalised, aes(x=normalised_cpg_count)) + 
  geom_histogram(color="black", fill="white")+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))+
  xlab("Number of Significant CpGs per TE")+
  ylab("Number of TEs")+
  xlim(0,0.15)

# Side question, can we put all three distributions together for num cpgs per feature
head(number_diff_cpgs_per_TE_normalised)
head(exon_genes_normalised)
head(promotor_genes_normalised)

TE_data <- number_diff_cpgs_per_TE_normalised[,c(11,17)]
exon_data <-exon_genes_normalised[,c(10,15)]
prom_data <-promotor_genes_normalised[,c(10,15)]
all_data <- rbind(TE_data,exon_data,prom_data)

ggplot(all_data, aes(x=normalised_cpg_count, fill = feature)) + 
  geom_histogram(colour="black")+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.title = element_blank())+
  xlab("Normalised Significant CpG Count")+
  ylab("Fequency")+
  xlim(0,0.15)+
  scale_fill_discrete(breaks = c("exon_first3","promotors_2000bp","TE"),
                   labels = c("Exons 1-3", "Promotors","TEs"))

TEs_subset <- TEs_subset[TEs_subset$uniquie_identified %in% TEs_with_3_cpgs$uniquie_identified,]
length(unique(TEs_subset$uniquie_identified)) #1061

hyper_in_female_TEs<- subset(TEs_subset, percent_meth_difference_of_feature < 0)
length(unique(hyper_in_female_TEs$uniquie_identified)) #1053/1061

hyper_in_male_TEs <- subset(TEs_subset, percent_meth_difference_of_feature > 0)
length(unique(hyper_in_male_TEs$uniquie_identified)) #8/1061
# Hyper in both TEs? None


hyper_in_male_TEs$origin <- "hypermeth_male"
hyper_in_female_TEs$origin <- "hypermeth_female"
male_data <- hyper_in_male_TEs[,c(4:7,15)]
female_data <- hyper_in_female_TEs[,c(4:7,15)]
TE_for_writing_out <- rbind(male_data, female_data)
TE_for_writing_out <- TE_for_writing_out[!duplicated(TE_for_writing_out),]
write.table(TE_for_writing_out, file="diff_meth_TEs.bed", col.names = F, row.names = F,
            sep="\t", quote = F)


par(mfrow=c(2,2))
boxplot(hyper_in_female_TEs$male_mean_weightedMeth, main="female hyper TEs", 
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_female_TEs$female_mean_weightedMeth, main="female hyper TEs", 
        xlab="female", ylim=c(0,1))
boxplot(hyper_in_male_TEs$male_mean_weightedMeth, main="male hyper TEs",
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_male_TEs$female_mean_weightedMeth,  main="male hyper TEs",
        xlab="female", ylim=c(0,1))

# Re-make where are the cpgs taking into account TEs overlapping genes
cpgs_in_exons_or_proms <- unique(cpgs_in_exons, cpgs_in_proms)
cpgs_not_unique_to_TEs <- cpgs_in_TEs[cpgs_in_TEs %in% cpgs_in_exons_or_proms]

output_1 <- output[!(output$cpg_position %in% cpgs_not_unique_to_TEs & output$feature == "TE"),]
output_1 <- output_1[!(output_1$cpg_position %in% cpgs_in_exons_or_proms 
                       & output_1$feature == "intergenic"),]

output_1$hyper_sex <- ifelse(output_1$male_mean_weightedMeth > output_1$female_mean_weightedMeth,
                             "male","female")
unique(output_1$hyper_sex)

ggplot(output_1, aes(x=feature,fill=hyper_sex))+
  geom_bar()+
  xlab("Feature")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  scale_x_discrete(breaks = c("exon_first3", "exon_notFirst3", "intron",
                              "promotors_2000bp","TE","intergenic"),
                   labels = c("Exons 1-3", "Exons 4+", "Introns",
                              "Promotors","TEs","Intergenic"),
                   limits = c("exon_first3","promotors_2000bp","TE","intergenic", 
                               "intron","exon_notFirst3"))+
  labs(fill='Hypermethylated\nSex')+
  theme(axis.title = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.text.x=element_text(angle=45,hjust=1,size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=22))+
  scale_fill_manual(limits=c("male","female"),
                    labels=c("Male","Female"),
                    values=c("steelblue1","lightpink"))
# This has worked it just doesn't make a big difference to the numbers


#Out of curiosity just plot CpGs unique to a single feature
#output_1 <- output[!duplicated(output$cpg_count),]
#ggplot(output_1, aes(x=feature))+
 # geom_bar()+
#  guides(fill=FALSE)+
 # xlab("Feature")+
#  ylab("Number of Significant CpGs")+
 # theme_bw()+
#  scale_x_discrete(breaks = c("exon_first3", "exon_notFirst3", "intron","Not_annotated","promotors_2000bp","TE"),
 #                  labels = c("Exons 1-3", "Exons 4+", "Introns", "No Annotation", "Promotors","TEs"),
  #                 limits = c("exon_first3","promotors_2000bp","TE", "Not_annotated", "intron","exon_notFirst3"))
# This has worked it just doesn't make a big difference to the numbers

## -------------------------------------------------------------------------
# Are these GENES hypermethylated in females or males, any with both? for promotos and exons
## -------------------------------------------------------------------------

head(promotor_genes)

hyper_in_female_promotor_genes <- subset(promotor_genes, percent_meth_difference_of_feature < 0)
length(unique(hyper_in_female_promotor_genes$gene_id)) #2645/2709

hyper_in_male_promotor_genes <- subset(promotor_genes, percent_meth_difference_of_feature > 0)
length(unique(hyper_in_male_promotor_genes$gene_id)) #64/2709

# Goodness of fit
observed = c(2645, 64)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) 

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

# Make a nicer boxplot for supplementary
female_hyper <- hyper_in_female_promotor_genes[,c(10,11)]
colnames(female_hyper) <- c("Male","Female")
female_hyper$category <- "female_hyper"
female_hyper_melt <- melt(female_hyper)
colnames(female_hyper_melt) <- c("Category","Sex","Weighted Methylation")
range(female_hyper$Female) #0.08219821 0.93270224
mean(female_hyper$Female)#0.6984731
sd(female_hyper$Female)#0.1838038

male_hyper <- hyper_in_male_promotor_genes[,c(10,11)]
colnames(male_hyper) <- c("Male","Female")
male_hyper$category <- "male_hyper"
male_hyper_melt <- melt(male_hyper)
colnames(male_hyper_melt) <- c("Category","Sex","Weighted Methylation")
range(male_hyper$Male) #0.04793464 0.34374614
mean(male_hyper$Male) # 0.1171903
sd(male_hyper$Male)#0.07344032

both <- rbind(female_hyper_melt, male_hyper_melt)

ggplot(both, aes(x=Category, y=`Weighted Methylation`, fill=Sex))+
  geom_boxplot()+
  theme_bw()+
  ggtitle("Promotors")+
  xlab("Hypermethylated Sex")+
  ylab("Mean Weighted Methylation")+
  scale_fill_manual("",breaks=c("Male", "Female"),
                    values = c("steelblue1","pink1"))+
  scale_x_discrete(limits=c("female_hyper", "male_hyper"),
                   labels=c("Female", "Male"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22),
        plot.title = element_text(size=22))


head(exon_genes)

hyper_in_female_exon_genes <- subset(exon_genes, percent_meth_difference_of_feature < 0)
length(unique(hyper_in_female_exon_genes$gene_id)) #2709/2736

hyper_in_male_exon_genes <- subset(exon_genes, percent_meth_difference_of_feature > 0)
length(unique(hyper_in_male_exon_genes$gene_id)) #33/2736

# Goodness of fit
observed = c(2709, 33)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) 

# Hyper in both sexes in exons? 6 genes!

boxplot(hyper_in_female_exon_genes$male_mean_weightedMeth, main="female hyper exons", 
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_female_exon_genes$female_mean_weightedMeth, main="female hyper exons", 
        xlab="female", ylim=c(0,1))
boxplot(hyper_in_male_exon_genes$male_mean_weightedMeth, main="male hyper exons",
        xlab="male", ylim=c(0,1))
boxplot(hyper_in_male_exon_genes$female_mean_weightedMeth,  main="male hyper exons",
        xlab="female", ylim=c(0,1))

# Make a nicer boxplot for supplementary
female_hyper <- hyper_in_female_exon_genes[,c(10,11)]
colnames(female_hyper) <- c("Male","Female")
female_hyper$category <- "female_hyper"
female_hyper_melt <- melt(female_hyper)
colnames(female_hyper_melt) <- c("Category","Sex","Weighted Methylation")
range(female_hyper$Female) #0.07474839 0.97086837
mean(female_hyper$Female)#0.7330598
sd(female_hyper$Female)#0.1521472

male_hyper <- hyper_in_male_exon_genes[,c(10,11)]
colnames(male_hyper) <- c("Male","Female")
male_hyper$category <- "male_hyper"
male_hyper_melt <- melt(male_hyper)
colnames(male_hyper_melt) <- c("Category","Sex","Weighted Methylation")
range(male_hyper$Male) #0.03368477 0.45143344
mean(male_hyper$Male) # 0.1437782
sd(male_hyper$Male)#0.1228486

both <- rbind(female_hyper_melt, male_hyper_melt)

ggplot(both, aes(x=Category, y=`Weighted Methylation`, fill=Sex))+
  geom_boxplot()+
  theme_bw()+
  ggtitle("Exons 1-3")+
  xlab("Hypermethylated Sex")+
  ylab("Mean Weighted Methylation")+
  scale_fill_manual("",breaks=c("Male", "Female"),
                    values = c("steelblue1","pink1"))+
  scale_x_discrete(limits=c("female_hyper", "male_hyper"),
                   labels=c("Female", "Male"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22),
        plot.title = element_text(size=22))


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

library(grid)
upset(main_data, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity",
      mainbar.y.label =NULL)
#grid.text("Exons 1-3",x = 0.15, y=0.95, gp=gpar(fontsize=36))
grid.text("Intersection Size",x = 0.45, y=0.60, gp=gpar(fontsize=20), rot = 90)

# overlap of female hypermeth, genes in genome 32590, prom:2709, exon:2736, overlap:1521
phyper(1521, 2736, 29881, 2709, lower.tail = F)

## -------------------------------------------------------------------------
# Write out all the gene lists for later use
## -------------------------------------------------------------------------

head(promotor_genes_only) # all promotor genes 2709
write.table(as.data.frame(promotor_genes_only), file="./diff_meth_gene_lists/diff_meth_promotor_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)

head(exon_genes_only) # all exon genes 2736
write.table(as.data.frame(exon_genes_only), file="./diff_meth_gene_lists/diff_meth_exon_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)

length(unique(TEs_subset$uniquie_identified)) # all TEs 824
write.table(as.data.frame(TEs_subset$uniquie_identified), file="./diff_meth_gene_lists/diff_meth_TEs_IDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)


# Also forthe below output the weighted meth scores for correlations
both # genes with both exon and promotor 1522
write.table(as.data.frame(both), file="./diff_meth_gene_lists/diff_meth_common_promotor_exon_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)

prom_only <- promotor_genes_only[!promotor_genes_only %in% both] # 1187
write.table(as.data.frame(prom_only), file="./diff_meth_gene_lists/diff_meth_unique_promotor_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)

exon_only <- exon_genes_only[!exon_genes_only %in% both] # 1214
write.table(as.data.frame(exon_only), file="./diff_meth_gene_lists/diff_meth_unique_exon_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)


both_with_weighted_meth <- promotor_genes[promotor_genes$gene_id %in% both,]
write.table(as.data.frame(both_with_weighted_meth), file="./diff_meth_gene_lists/diff_meth_common_promotor_exon_all_info.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

prom_with_weighted_meth <- promotor_genes[promotor_genes$gene_id %in% prom_only,]
write.table(as.data.frame(prom_with_weighted_meth), file="./diff_meth_gene_lists/diff_meth_unique_promotor_all_info.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

exon_with_weighted_meth <- exon_genes[exon_genes$gene_id %in% exon_only,]
write.table(as.data.frame(exon_with_weighted_meth), file="./diff_meth_gene_lists/diff_meth_unique_exon_all_info.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

# Make file of all proms and if they're diff or not for Peter
prom_with_weighted_meth <- promotor_genes[promotor_genes$gene_id %in% promotor_genes_only,]
head(prom_with_weighted_meth)
for_peter <- prom_with_weighted_meth[,c(4,5,6,7,9,10,11)]
head(for_peter2)

for_peter2 <- for_peter[!duplicated(for_peter),]
for_peter2$sig_diff_methylated <- "yes"

head(annotation)
annotation_proms <- subset(annotation, feature=="promotors_2000bp")
annotation_proms <- annotation_proms[,c(1,2,3,4,5,7,8)]
nrow(annotation_proms)

not_sig_proms <- annotation_proms[!annotation_proms$gene_id %in% for_peter2$gene_id,]
head(not_sig_proms)
not_sig_proms$sig_diff_methylated <- "no"

all_data_proms <- rbind(for_peter2, not_sig_proms)
head(all_data_proms)
nrow(all_data_proms)
write.table(as.data.frame(all_data_proms), file="pcitri_diff_meth_promotors.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

FPKM_values_by_sex <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/FPKM_values_by_sex.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
data_wide <- dcast(FPKM_values_by_sex, gene_id ~ origin, value.var="FPKM")

head(all_data_proms)
all_data_proms$hyper_status <- ifelse(all_data_proms$male_mean_weightedMeth >
                                        all_data_proms$female_mean_weightedMeth,
                                  "male", "female")

all_data_proms$hyper_status[all_data_proms$sig_diff_methylated == "no"] <- "not_sig"

both_data <- merge(all_data_proms, data_wide, by="gene_id")
head(both_data)
both_data <- both_data[,c(9,10,11)]
both_data_melt <- melt(both_data)
head(both_data_melt)

#### Define summary function (ref:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 

summary_all<-summarySE(both_data_melt, measurevar = "value", 
                       groupvars = c("hyper_status","variable"))
head(summary_all)

ggplot(summary_all, aes(x=hyper_status, y=value, fill=variable))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Gene Set")+
  ylab("FPKM")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("female","male","not_sig"),
                   labels = c("Female Hypermethylated","Male Hypermethylated",
                              "Not Significant"))

both_data <- merge(all_data_proms, data_wide, by="gene_id")
head(both_data)
both_data$hyper_status[(both_data$male_mean_weightedMeth < 0.01 &
                                both_data$female_mean_weightedMeth < 0.01)]<-"not_meth"
both_data <- both_data[,c(9,10,11)]
both_data_melt <- melt(both_data)
head(both_data_melt)

summary_all<-summarySE(both_data_melt, measurevar = "value", 
                       groupvars = c("hyper_status","variable"))
head(summary_all)

ggplot(summary_all, aes(x=hyper_status, y=value, fill=variable))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  ggtitle("Promotors")+
  xlab("Gene Set")+
  ylab("FPKM")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("female","male","not_sig","not_meth"),
                   labels = c("Female\nHypermethylated","Male\nHypermethylated",
                              "Not Significant","Unmethylated"),
                   limits = c("not_meth","not_sig","female","male"))
