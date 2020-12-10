## -------------------------------------------------------------------------
# Genome browser for methylation - lollipops
## -------------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation/differential_meth_methylkit")

library(readr)
library(reshape2)

## -------------------------------------------------------------------------
# Raw methylation data from methylkit
F_vs_M_objectmethbase <- read_delim("F_vs_M_objectmethbase.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
head(F_vs_M_objectmethbase)

## -------------------------------------------------------------------------
# Make the file more useable: average female / male proportion meth C per positon
F_vs_M_objectmethbase <- F_vs_M_objectmethbase[,-c(3,4,7,10,13,16,19,22,25,28,31)]
head(F_vs_M_objectmethbase)

F_vs_M_objectmethbase$female_coverage <- F_vs_M_objectmethbase$coverage1 + F_vs_M_objectmethbase$coverage2 +
                                             F_vs_M_objectmethbase$coverage3 + F_vs_M_objectmethbase$coverage4 +
                                             F_vs_M_objectmethbase$coverage5

F_vs_M_objectmethbase$male_coverage <- F_vs_M_objectmethbase$coverage6 + F_vs_M_objectmethbase$coverage7 +
                                             F_vs_M_objectmethbase$coverage8 + F_vs_M_objectmethbase$coverage9


F_vs_M_objectmethbase$female_cs <- F_vs_M_objectmethbase$numCs1 + F_vs_M_objectmethbase$numCs2 + F_vs_M_objectmethbase$numCs3 +
                                       F_vs_M_objectmethbase$numCs4 + F_vs_M_objectmethbase$numCs5

F_vs_M_objectmethbase$male_cs <- F_vs_M_objectmethbase$numCs6 + F_vs_M_objectmethbase$numCs7 +
                                       F_vs_M_objectmethbase$numCs8 + F_vs_M_objectmethbase$numCs9

F_vs_M_objectmethbase <- F_vs_M_objectmethbase[, c(1,2,21,22,23,24)]

F_vs_M_objectmethbase$female_proportion <- 1- (F_vs_M_objectmethbase$female_coverage - F_vs_M_objectmethbase$female_cs)/ F_vs_M_objectmethbase$female_coverage
F_vs_M_objectmethbase$male_proportion <- 1- (F_vs_M_objectmethbase$male_coverage - F_vs_M_objectmethbase$male_cs)/ F_vs_M_objectmethbase$male_coverage

head(F_vs_M_objectmethbase)
F_vs_M_objectmethbase <- F_vs_M_objectmethbase[, c(1,2,7,8)]
colnames(F_vs_M_objectmethbase) <- c("chromosome","position","female_proportion_methylated","male_proportion_methylated")

write.table(F_vs_M_objectmethbase, file = "proportion_methylation_per_C_whole_genome_by_sex.txt", sep='\t',
            quote = F, col.names = T, row.names = F)

## -------------------------------------------------------------------------
# Now need to change the file to fit the requirments of: http://maplab.imppc.org/methylation_plotter/

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation/regional_genome_browser")

data1 <- F_vs_M_objectmethbase
head(data1)

# Try with a subset
subset1 <- head(data1, n = 100)
head(subset1)
subset1 <- subset1[,-1]
colnames(subset1) <- c("position","Female","Male")

# Flip columns and rows
df_transpose = as.data.frame(t(subset1))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="positions_genome_browser_methylation.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)

# Just uploaded the test to the website above and it works !!!!!
# Max num of CpGs is 100, so I'll make a few (jusr change region and output file name)
region1 <- data1[c(900:1000),]
head(region1)
region1 <- region1[,-1]
colnames(region1) <- c("position","Female","Male")

# Flip columns and rows
df_transpose = as.data.frame(t(region1))
df_transpose <- df_transpose[c(2,3,1),]
write.table(df_transpose, file="region10_genome_browser.txt", sep = "\t",
            quote = F, col.names = F, row.names = T)
