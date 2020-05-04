# ------------------------------------------------------
# Making base genome-wide GO file and background sets
# ------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses")
library(readr)
library(tidyr)

# ------------------------------------------------------
# Make a compatible GO annotation file
Pcitri_genes_with_GO <- read_table2("Pcitri_genes_with_GO.txt", 
                                    col_names = FALSE)
colnames(Pcitri_genes_with_GO) <- c("geneID","goID")
new_annotations <- separate_rows(Pcitri_genes_with_GO, goID, sep =';')
write.table(new_annotations, file="P_citri_GO_terms.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

# ------------------------------------------------------
# Make the expression gene lists
logFC_DEgenes <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/logFC_DEgenes.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
logFC_DEgenes$X9[logFC_DEgenes$X8 > 10] <- "Ex_Fbias"
logFC_DEgenes$X9[logFC_DEgenes$X8 < -10] <- "Ex_Mbias"

all_genes_in_RNAseq <- as.data.frame(unique(logFC_DEgenes$X1)) #15570
write.table(all_genes_in_RNAseq, file="./gene_lists/all_genes_in_RNAseq.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

unbias_genes <- subset(logFC_DEgenes, X9 =="Unbias") #5022
write.table(as.data.frame(unbias_genes$X1), file="./gene_lists/unibas_exp_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)
all_bias_genes <- subset(logFC_DEgenes, !X9 =="Unbias") #10548
write.table(as.data.frame(all_bias_genes$X1), file="./gene_lists/all_bias_exp_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

mbias_genes <- subset(logFC_DEgenes, X9 =="Mbias") #4756
write.table(as.data.frame(mbias_genes$X1), file="./gene_lists/mbias_exp_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)
fbias_genes <- subset(logFC_DEgenes, X9 =="Fbias")#5270
write.table(as.data.frame(fbias_genes$X1), file="./gene_lists/fbias_exp_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

exmbias_genes <- subset(logFC_DEgenes, X9 =="Ex_Mbias")#344
write.table(as.data.frame(exmbias_genes$X1), file="./gene_lists/ex_mbias_exp_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)
exfbias_genes <- subset(logFC_DEgenes, X9 =="Ex_Fbias")#178
write.table(as.data.frame(exfbias_genes$X1), file="./gene_lists/ex_fbias_exp_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

# Make general methylated gene lists
weightedMeth_exons_promotors_only <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/weightedMeth_exons_promotors_only.txt", 
                                                "\t", escape_double = FALSE, trim_ws = TRUE)

proms <- subset(weightedMeth_exons_promotors_only, feature =="promotors_2000bp")
exons <- subset(weightedMeth_exons_promotors_only, feature =="exon_first3")

methylated_proms <- subset(proms, weightedMeth.mean > 0.1)
methylated_proms <- as.data.frame(unique(methylated_proms$gene_id)) #7425
write.table(as.data.frame(methylated_proms), file="./gene_lists/genes_with_methylated_proms_in_either_sex.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

methylated_exons <- subset(exons, weightedMeth.mean > 0.1) # 7455
methylated_exons <- as.data.frame(unique(methylated_exons$gene_id))
write.table(as.data.frame(methylated_exons), file="./gene_lists/genes_with_methylated_exons_in_either_sex.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

# Make diff meth lists
diff_meth <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/meth_paired_with_exp/diff_meth_genes_labeled_with_hyper_by_sex.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

all_genes_in_meth <- as.data.frame(unique(diff_meth$gene_id)) #35782
write.table(as.data.frame(all_genes_in_meth), file="./gene_lists/all_genes_in_meth_data.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

not_diff_meth <- subset(diff_meth, diff_meth_feature =="not_diff_meth")
not_diff_meth <- as.data.frame(unique(not_diff_meth$gene_id)) #31898
write.table(as.data.frame(not_diff_meth), file="./gene_lists/not_diff_meth_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

diff_meth <- subset(diff_meth, !diff_meth_feature =="not_diff_meth")
diff_meth <- as.data.frame(unique(diff_meth$gene_id)) #3923
write.table(as.data.frame(diff_meth), file="./gene_lists/diff_meth_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

prom_female_meth <- subset(diff_meth, diff_meth_feature =="prom_female_meth")
prom_female_meth <- as.data.frame(unique(prom_female_meth$gene_id)) #2650
write.table(as.data.frame(prom_female_meth), file="./gene_lists/prom_female_hypermeth_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)
prom_male_meth <- subset(diff_meth, diff_meth_feature =="prom_male_meth")
prom_male_meth <- as.data.frame(unique(prom_male_meth$gene_id)) #236
write.table(as.data.frame(prom_male_meth), file="./gene_lists/prom_male_hypermeth_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)


exon_female_meth <- subset(diff_meth, diff_meth_feature =="exon_female_meth")
exon_female_meth <- as.data.frame(unique(exon_female_meth$gene_id)) #2706
write.table(as.data.frame(exon_female_meth), file="./gene_lists/exon_female_hypermeth_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)
exon_male_meth <- subset(diff_meth, diff_meth_feature =="exon_male_meth")
exon_male_meth <- as.data.frame(unique(exon_male_meth$gene_id)) #211
write.table(as.data.frame(exon_male_meth), file="./gene_lists/exon_male_hypermeth_genes.txt", sep="\t", quote = F,
            col.names = F, row.names = F)

# ------------------------------------------------------
# Read in my new annotation file to check it's all good

P_citri_GO_terms <- read_delim("P_citri_GO_terms.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

# ------------------------------------------------------
# Make background go sets

colnames(all_genes_in_meth) <- "geneID"
all_meth_background <- merge(P_citri_GO_terms, all_genes_in_meth)
write.table(as.data.frame(all_meth_background), file="./background_go_sets/all_meth_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

colnames(all_genes_in_RNAseq) <- "geneID"
background <- merge(P_citri_GO_terms, all_genes_in_RNAseq)
write.table(as.data.frame(background), file="./background_go_sets/all_RNAseq_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

all_bias_genes <- as.data.frame(unique(all_bias_genes$X1))
colnames(all_bias_genes) <- "geneID"
background <- merge(P_citri_GO_terms, all_bias_genes)
write.table(as.data.frame(background), file="./background_go_sets/all_biased_exp_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

colnames(diff_meth) <- "geneID"
background <- merge(P_citri_GO_terms, diff_meth)
write.table(as.data.frame(background), file="./background_go_sets/all_diff_meth_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

# background sets for comparing the extreme exp bias against just exp bias
fbias_genes <- as.data.frame(unique(fbias_genes$X1))
colnames(fbias_genes) <- "geneID"
exfbias_genes <- as.data.frame(unique(exfbias_genes$X1))
colnames(exfbias_genes) <- "geneID"
all_female_bias <- rbind(fbias_genes, exfbias_genes)#5448
background <- merge(P_citri_GO_terms, all_female_bias)
write.table(as.data.frame(background), file="./background_go_sets/all_female_bias_exp_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

mbias_genes <- as.data.frame(unique(mbias_genes$X1))
colnames(mbias_genes) <- "geneID"
exmbias_genes <- as.data.frame(unique(exmbias_genes$X1))
colnames(exmbias_genes) <- "geneID"
all_male_bias <- rbind(mbias_genes, exmbias_genes)#5100
background <- merge(P_citri_GO_terms, all_male_bias)
write.table(as.data.frame(background), file="./background_go_sets/all_male_bias_exp_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

###
list_diff_iso_genes <- read_csv("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/list_diff_iso_genes.txt")
colnames(list_diff_iso_genes) <- "geneID"
background <- merge(P_citri_GO_terms, list_diff_iso_genes)
write.table(as.data.frame(background), file="./background_go_sets/all_alt_spliced_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

FPKMs_logFC_bias_catergory <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/FPKMs_logFC_bias_catergory.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
FPKMs_logFC_bias_catergory<-as.data.frame(FPKMs_logFC_bias_catergory$gene_id[!FPKMs_logFC_bias_catergory$bias=="Unbias"])
colnames(FPKMs_logFC_bias_catergory) <- "geneID"
background <- merge(P_citri_GO_terms, FPKMs_logFC_bias_catergory)
write.table(as.data.frame(background), file="./background_go_sets/all_diff_exp_inc_sex_limited.txt", sep="\t", quote = F,
            col.names = T, row.names = F)
