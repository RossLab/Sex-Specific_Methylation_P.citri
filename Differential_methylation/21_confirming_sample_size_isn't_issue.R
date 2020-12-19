## -------------------------------------------------------------------------
# Comparing male meth vs 4x female samples pooled
## -------------------------------------------------------------------------
# Doing to to show the large amount of male meth that randomly distributed isn't
# due to the number of individuals used during the pooling set for DNA extraction
# pooling 4x female samples = 60indiv, one male sample = 60indiv

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation/differential_meth_methylkit")

library(readr)
library(reshape2)

## -------------------------------------------------------------------------
# Raw methylation data from methylkit
F_vs_M_objectmethbase <- read_delim("F_vs_M_objectmethbase.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
head(F_vs_M_objectmethbase)

F_vs_M_objectmethbase <- F_vs_M_objectmethbase[,-c(3,4,7,10,13,16,19,22,25,28,31)]
head(F_vs_M_objectmethbase)

# Make a female data set with 4 samples
F_vs_M_objectmethbase$female_coverage_4_samples <- F_vs_M_objectmethbase$coverage1 + F_vs_M_objectmethbase$coverage2 +
  F_vs_M_objectmethbase$coverage3 + F_vs_M_objectmethbase$coverage4

F_vs_M_objectmethbase$female_cs <- F_vs_M_objectmethbase$numCs1 + F_vs_M_objectmethbase$numCs2 + F_vs_M_objectmethbase$numCs3 +
  F_vs_M_objectmethbase$numCs4 

# make a male data set with 1 sample: test with a couple of them
F_vs_M_objectmethbase$male_coverage_1_sample <- F_vs_M_objectmethbase$coverage6 

F_vs_M_objectmethbase$male_cs <- F_vs_M_objectmethbase$numCs6 

F_vs_M_objectmethbase <- F_vs_M_objectmethbase[, c(1,2,21,22,23,24)]

F_vs_M_objectmethbase$female_proportion <- 1- (F_vs_M_objectmethbase$female_coverage_4_samples - F_vs_M_objectmethbase$female_cs)/ F_vs_M_objectmethbase$female_coverage_4_samples
F_vs_M_objectmethbase$male_proportion <- 1- (F_vs_M_objectmethbase$male_coverage_1_sample - F_vs_M_objectmethbase$male_cs)/ F_vs_M_objectmethbase$male_coverage_1_sample

head(F_vs_M_objectmethbase)

female <- F_vs_M_objectmethbase[,c(1,2,3,7)]
male <- F_vs_M_objectmethbase[,c(1,2,5,8)]
head(female)
head(male)


# ID, chrom, loc.start, loc.end, num.mark (coverage), seg.mean(proportion C)
colnames(female) <- c("chrom","loc.start","num.mark","seg.mean")
female$ID <- "Female"
female$loc.end <- female$loc.start + 1
female <- female[,c(5,1,2,6,3,4)]
write.table(female, file="female_4_samples_meth.seg", sep = "\t",
            quote = F, col.names = T, row.names = F)

colnames(male) <- c("chrom","loc.start","num.mark","seg.mean")
male$ID <- "Male"
male$loc.end <- male$loc.start + 1
male <- male[,c(5,1,2,6,3,4)]
write.table(male, file="male_1_sample_meth.seg", sep = "\t",
            quote = F, col.names = T, row.names = F)

# Put these into IGV and look