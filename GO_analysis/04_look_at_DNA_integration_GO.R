# Check with diff meth genes are related to the GO term realted to TEs

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs")

library(readr)

Pcitri_genes_with_GO <- read_delim("GO_analyses/Pcitri_genes_with_GO.txt", 
                                   "\t", escape_double = FALSE, col_names = FALSE, 
                                   trim_ws = TRUE)



exon_female_hypermeth_genes <- read_csv("GO_analyses/gene_lists/exon_female_hypermeth_genes.txt", 
                                        col_names = FALSE)

hyper_genes_with_GO <- merge(exon_female_hypermeth_genes, Pcitri_genes_with_GO, by="X1")
# 355 / 585 of these genes have the GO terms DNA integration GO:0015074


prom_female_hypermeth_genes <- read_csv("GO_analyses/gene_lists/prom_female_hypermeth_genes.txt", 
                                        col_names = FALSE)
hyper_genes_with_GO2 <- merge(prom_female_hypermeth_genes, Pcitri_genes_with_GO, by="X1")
# 215 / 464 of these genes have the GO terms DNA integration GO:0015074

in_both <- merge(hyper_genes_with_GO, hyper_genes_with_GO2, by ="X1")
with_go_term <- in_both[grepl("GO:0015074", in_both), ] #214
just_genes <- as.data.frame(with_go_term$X1)
colnames(just_genes) <- "gene_id"

# Have a look at the expression of those 214 genes 
FPKMs_logFC_bias_catergory <- read_delim("transcription/FPKMs_logFC_bias_catergory.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)

genes_with_exp <- merge(just_genes, FPKMs_logFC_bias_catergory, by = "gene_id")
# Lol only 23 of them are expressed
table(genes_with_exp$bias) # 13/23 are male expressed
