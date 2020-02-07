#---------------------------------------------------------------
# Bed files for motif analysis
#---------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/motif_analysis")
library(readr)

annotation <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/PCITRI_weighted_meth_annotation_common_both_sexes.txt", 
                                                                "\t", escape_double = FALSE, trim_ws = TRUE)
annotation$feature <- as.factor(annotation$feature)
levels(annotation$feature)

#---------------------------------------------------------------
# General exons 1-3 vs general exons 3+
#---------------------------------------------------------------

exons_first3 <- subset(annotation, feature == "exon_first3")
exons_3plus <- subset(annotation, feature == "exon_notFirst3")

colnames(exons_first3) <- NULL
colnames(exons_3plus) <- NULL

write.table(exons_first3[,c(4,2,3)], file = "all_exons_1-3.bed",
            sep="\t", quote = F, col.names = F, row.names = F)
write.table(exons_3plus[,c(4,2,3)], file = "all_exons_3plus.bed",
            sep="\t", quote = F, col.names = F, row.names = F)

#---------------------------------------------------------------
# Methylated and unmethylated annotations for comparisons
#---------------------------------------------------------------

# Gettgin exons 1-3 meth and un-meth
exons_first3 <- subset(annotation, feature == "exon_first3")

exons_first3_methylations <- subset(exons_first3, methylated == "yes")
exons_first3_unmeth <- subset(exons_first3, methylated == "no")

colnames(exons_first3_methylations) <- NULL
colnames(exons_first3_unmeth) <- NULL

write.table(exons_first3_methylations[,c(4,2,3)], file = "methylated_exons_1-3.bed",
            sep="\t", quote = F, col.names = F, row.names = F)
write.table(exons_first3_unmeth[,c(4,2,3)], file = "unmethylated_exons_1-3.bed",
            sep="\t", quote = F, col.names = F, row.names = F)

# All other exons meth and un-meth
exons_3plus <- subset(annotation, feature == "exon_notFirst3")

exons_3plus_methylations <- subset(exons_3plus, methylated == "yes")
exons_3plus_unmeth <- subset(exons_3plus, methylated == "no")

colnames(exons_3plus_methylations) <- NULL
colnames(exons_3plus_unmeth) <- NULL

write.table(exons_3plus_methylations[,c(4,2,3)], file = "methylated_exons_3plus.bed",
            sep="\t", quote = F, col.names = F, row.names = F)
write.table(exons_3plus_unmeth[,c(4,2,3)], file = "unmethylated_exons_3plus.bed",
            sep="\t", quote = F, col.names = F, row.names = F)


# Promotors meth and unmeth
proms <- subset(annotation, feature == "promotors_1000bp")

proms_methylations <- subset(proms, methylated == "yes")
proms_unmeth <- subset(proms, methylated == "no")

colnames(proms_methylations) <- NULL
colnames(proms_unmeth) <- NULL

write.table(proms_methylations[,c(4,2,3)], file = "methylated_promotors.bed",
            sep="\t", quote = F, col.names = F, row.names = F)
write.table(proms_unmeth[,c(4,2,3)], file = "unmethylated_promotors.bed",
            sep="\t", quote = F, col.names = F, row.names = F)
