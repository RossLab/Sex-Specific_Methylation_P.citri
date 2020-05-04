#--------------------------------------------------------------------
# Making intergenic regions or unannotated regions to quantify background methylation
#--------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation")

library(readr)
library(dplyr)

merged_annotations <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files/merged_annotations.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
head(merged_annotations)

# We can remove some annotations which are within others to decrease the file size
annotations_cleaner <- merged_annotations[merged_annotations$feature == "promotors_2000bp"|
                                            merged_annotations$feature =="whole_gene" |
                                            merged_annotations$feature == "TE",]
annotations_cleaner <- annotations_cleaner[,-c(1,5)]
              
# Order the file and make in a format the python script wants
ordered <- annotations_cleaner %>% arrange(scaffold, start)
head(ordered)
colnames(ordered) <- c("chr","start","end","gene")
ordered <- ordered[,c(4,1,2,3)]

write.table(ordered, file="annotations_for_getting_intergenic_regions.gtf",
            col.names = T, row.names = F, quote = F, sep = '\t')

# Run getting_intergenic_regions.py as:
# python getting_intergenic_regions.py annotations_for_getting_intergenic_regions.gtf > intergenic.bed

intergenic <- read_delim("intergenic.bed", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

# Remove rows which are incorrect becasue they wer created when a TE is found
# within gene conordinates
head(intergenic)

intergenic$remove <- ifelse(intergenic$X2 > intergenic$X3, "yes","no")
intergenic_for_sure <- intergenic[intergenic$remove == "no",]

# Neaten the file up for later use
intergenic_for_sure$X4 <- "intergenic"
intergenic_for_sure <- intergenic_for_sure[,-5]
colnames(intergenic_for_sure) <- c("scaffold","start","end","feature")

write.table(intergenic_for_sure, file="PCITRI.assembly.v0.intergenic.txt",
            col.names = T, row.names = F, quote = F, sep = '\t')