## -------------------------------------------------------------------------
# Genome browser for methylation 
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




## -------------------------------------------------------------------------
# New attempt make a .igv file to use with the online IGV
# Should be: Chromosome, Start, End, Feature, Female_Level, Male_Level

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation/regional_genome_browser")

meth_info <- read_delim("proportion_methylation_per_C_whole_genome_by_sex.txt", 
                                                               "\t", escape_double = FALSE, trim_ws = TRUE)
head(meth_info)
colnames(meth_info)<- c("Chromosome","Start","Female","Male")
meth_info$End <- meth_info$Start + 1
meth_info$Name <- "CpG"
meth_info <- meth_info[,c(1,2,5,6,3,4)]
write.table(meth_info, file="meth_levels_by_sex.igv", sep = "\t",
            quote = F, col.names = T, row.names = F) # Not working

# Try bed format
meth_info <- meth_info[,-4]
meth_info1 <- melt(meth_info, id.vars = c("Chromosome","Start","End"))
head(meth_info1)

female <- subset(meth_info1, variable =="Female")
male <- subset(meth_info1, variable =="Male")

write.table(female, file="female_meth.bed", sep = "\t",
            quote = F, col.names = F, row.names = F) # works but doens't give me a level visualised by shade

# Try seg format
# ID, chrom, loc.start, loc.end, num.mark (coverage), seg.mean(proportion C)
head(meth_info1)

F_vs_M_objectmethbase <- read_delim("F_vs_M_objectmethbase.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
head(F_vs_M_objectmethbase)
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
F_vs_M_objectmethbase <- F_vs_M_objectmethbase[,-c(5,6)]

female <- F_vs_M_objectmethbase[,-c(4,6)]
male <- F_vs_M_objectmethbase[,-c(3,5)]
head(female)
head(male)

# ID, chrom, loc.start, loc.end, num.mark (coverage), seg.mean(proportion C)
colnames(female) <- c("chrom","loc.start","num.mark","seg.mean")
female$ID <- "Female"
female$loc.end <- female$loc.start + 1
female <- female[,c(5,1,2,6,3,4)]
write.table(female, file="female_meth.seg", sep = "\t",
            quote = F, col.names = T, row.names = F)

colnames(male) <- c("chrom","loc.start","num.mark","seg.mean")
male$ID <- "Male"
male$loc.end <- male$loc.start + 1
male <- male[,c(5,1,2,6,3,4)]
write.table(male, file="male_meth.seg", sep = "\t",
            quote = F, col.names = T, row.names = F)

# These work, wack them into IGV along with the genome .fa and .fa.fai and the .gff3
# Try just chr1 for both
female_subset <- subset(female, chrom =="PCITRI_00005")
female_subset$seg.mean<-as.numeric(female_subset$seg.mean)
write.table(female_subset, file="female_meth_subset_chr5.seg", sep = "\t",
            quote = F, col.names = T, row.names = F)
# ' #type=DNA_METHYLATION ' between the ' needs to be included in the second line

male_subset <- subset(male, chrom =="PCITRI_00005")
male_subset$seg.mean<-as.numeric(male_subset$seg.mean)
write.table(male_subset, file="male_meth_subset_chr5.seg", sep = "\t",
            quote = F, col.names = T, row.names = F)

# FINALLY this works on the downloaded version of IGV