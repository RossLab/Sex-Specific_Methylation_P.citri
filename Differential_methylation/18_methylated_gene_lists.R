## -------------------------------------------------------------------------
# Getting methylated gene lists
## -------------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation")
library(readr)
library(doBy)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Hmisc)
library(scales)
## -------------------------------------------------------------------------
annotation <- read_delim("PCITRI_weighted_meth_annotation_by_sex.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
annotation$feature <- as.factor(annotation$feature)
annotation_males <- subset(annotation, origin =="male")
annotation_females <- subset(annotation, origin =="female")

both <- merge(annotation_males, annotation_females, by=c("gene_id","start","end"))
both <- both[,-c(7,9,10,11,12)]
colnames(both) <- c("gene_id","start","end","scaffold","feature","cpg_count","male_mean_weightedMeth","female_mean_weightedMeth")

head(both)

## -------------------------------------------------------------------------
#Genes with meth proms
## -------------------------------------------------------------------------
proms <- both[both$feature == "promotors_2000bp",]

proms$male_category[proms$male_mean_weightedMeth == 0] <- "none"
proms$male_category[proms$male_mean_weightedMeth > 0 &
                      proms$male_mean_weightedMeth < 0.3] <- "low"
proms$male_category[proms$male_mean_weightedMeth > 0.3 &
                      proms$male_mean_weightedMeth < 0.7] <- "medium"
proms$male_category[proms$male_mean_weightedMeth > 0.7] <- "high"
table(proms$male_category)
#   high    low medium   none 
# 303  28949   2283    220 


proms$female_category[proms$female_mean_weightedMeth == 0] <- "none"
proms$female_category[proms$female_mean_weightedMeth > 0 &
                      proms$female_mean_weightedMeth < 0.3] <- "low"
proms$female_category[proms$female_mean_weightedMeth > 0.3 &
                      proms$female_mean_weightedMeth < 0.7] <- "medium"
proms$female_category[proms$female_mean_weightedMeth > 0.7] <- "high"
table(proms$female_category)
#high    low medium   none 
# 1454  27973   1991    340 

## -------------------------------------------------------------------------
# Stats

# Need to compare proportions but chi-squared goodness of fit needs expected
# proportions which we don't have. Looked up and I should use the two-proportions z-test
nrow(proms)

# High= sig more female high
high <- prop.test(x = c(303, 1454), n = c(31758, 31758))
a <- high$p.value

# Medium = sig more male medium
medium <- prop.test(x = c(2283, 1991), n = c(31758, 31758))
b <- medium$p.value

# Low = sig more male low
low <- prop.test(x = c(28949, 27973), n = c(31758, 31758))
c <- low$p.value

# Low = sig more female none
none <- prop.test(x = c(220, 340), n = c(31758, 31758))
d <- none$p.value

## -------------------------------------------------------------------------
#Genes with meth exons
## -------------------------------------------------------------------------
exons <- both[both$feature == "exon_first3",]
exons$gene_id <- gsub("^.*Parent=", "", exons$gene_id)
exons$gene_id <- gsub("\\..*", "", exons$gene_id)
exons_males <- aggregate(male_mean_weightedMeth ~ gene_id, data=exons,
          FUN = mean)
exons_females <- aggregate(female_mean_weightedMeth ~ gene_id, data=exons,
                         FUN = mean)


exons_males$male_category[exons_males$male_mean_weightedMeth == 0] <- "none"
exons_males$male_category[exons_males$male_mean_weightedMeth > 0 &
                            exons_males$male_mean_weightedMeth < 0.3] <- "low"
exons_males$male_category[exons_males$male_mean_weightedMeth > 0.3 &
                            exons_males$male_mean_weightedMeth < 0.7] <- "medium"
exons_males$male_category[exons_males$male_mean_weightedMeth > 0.7] <- "high"
table(exons_males$male_category)
#high    low medium   none 
#366  29861   2851    176

exons_females$female_category[exons_females$female_mean_weightedMeth == 0] <- "none"
exons_females$female_category[exons_females$female_mean_weightedMeth > 0 &
                                exons_females$female_mean_weightedMeth < 0.3] <- "low"
exons_females$female_category[exons_females$female_mean_weightedMeth > 0.3 &
                                exons_females$female_mean_weightedMeth < 0.7] <- "medium"
exons_females$female_category[exons_females$female_mean_weightedMeth > 0.7] <- "high"
table(exons_females$female_category)
#high    low medium   none 
#1815  28862   2221    357 


## -------------------------------------------------------------------------
# Stats
nrow(exons_females)
nrow(exons_males)

# High= sig more female high
high <- prop.test(x = c(366, 1815), n = c(33255, 33255))
e <- high$p.value

# Medium = sig more male medium
medium <- prop.test(x = c(2221, 2851), n = c(33255, 33255))
f <- medium$p.value

# Low = sig more male low
low <- prop.test(x = c(28862, 29861), n = c(31758, 31758))
g <- low$p.value

# Low = sig more female none
none <- prop.test(x = c(357, 176), n = c(31758, 31758))
h <- none$p.value


# correction for multiple testing all of the above:
ps <- c(a,b,c,d,e,f,g,h)
p.adjust(ps, method = "hochberg", n = length(ps))

## -------------------------------------------------------------------------
# Need to pull out the gene list for each and the total gene list as background for the 
# GO enrichment analysis
## -------------------------------------------------------------------------
head(exons_males)
head(exons_females)
head(proms)

write.table(unique(proms$gene_id[proms$male_category=="none"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            no_male_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(proms$gene_id[proms$male_category=="low"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            low_male_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(proms$gene_id[proms$male_category=="medium"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            medium_male_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(proms$gene_id[proms$male_category=="high"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            high_male_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)

write.table(unique(proms$gene_id[proms$female_category=="none"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            no_female_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(proms$gene_id[proms$female_category=="low"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            low_female_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(proms$gene_id[proms$female_category=="medium"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            medium_female_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(proms$gene_id[proms$female_category=="high"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            high_female_prom_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)


write.table(unique(exons_males$gene_id[exons_males$male_category=="none"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            no_male_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(exons_males$gene_id[exons_males$male_category=="low"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            low_male_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(exons_males$gene_id[exons_males$male_category=="medium"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            medium_male_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(exons_males$gene_id[exons_males$male_category=="high"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            high_male_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)

write.table(unique(exons_males$gene_id[exons_females$female_category=="none"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            no_female_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(exons_males$gene_id[exons_females$female_category=="low"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            low_female_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(exons_males$gene_id[exons_females$female_category=="medium"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            medium_female_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(exons_males$gene_id[exons_females$female_category=="high"]), 
            file="~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses/gene_lists/methylaetd_genes/
            high_female_exon_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)


