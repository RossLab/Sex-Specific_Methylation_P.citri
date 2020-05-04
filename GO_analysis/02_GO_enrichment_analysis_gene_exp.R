#-----------------------------------------------
# Script created by Alun Jones, see paper Bebane et al. (2019) Neonics and bumblebees...
#-----------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses")
library(GOstats)
library(GSEABase)
library(treemap)
library(readr)

#-----------------------------------------------
# Read in background GO set and make compatible with GOstats
#GO_annotations <- read.table("./background_go_sets/all_RNAseq_genes_background.txt",header = T) # for sex-biased
GO_annotations <- read.table("./background_go_sets/all_diff_exp_inc_sex_limited.txt",header = T) #for extreme and limited
GO_annotations[,3] <- paste("IEA")
names(GO_annotations) <- c("genes","GOIds","evi")
GO_annotations[,3] <- paste("IEA")
GO_annotations <- GO_annotations[c(2,3,1)]


#-----------------------------------------------
# Create necessary objects

GO_frame <- GOFrame(GO_annotations,organism = "Planococcus citri")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))


#-----------------------------------------------
# Read in gene's of interest 
FPKMs_logFC_bias_catergory <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/FPKMs_logFC_bias_catergory.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
table(FPKMs_logFC_bias_catergory$bias)

#my_genes <- as.data.frame(FPKMs_logFC_bias_catergory$gene_id[FPKMs_logFC_bias_catergory$bias=="Mbias"])
#my_genes <- as.data.frame(FPKMs_logFC_bias_catergory$gene_id[FPKMs_logFC_bias_catergory$bias=="Fbias"])

#my_genes <- as.data.frame(FPKMs_logFC_bias_catergory$gene_id[FPKMs_logFC_bias_catergory$bias=="extreme_male_biased"])
#my_genes <- as.data.frame(FPKMs_logFC_bias_catergory$gene_id[FPKMs_logFC_bias_catergory$bias=="extreme_female_biased"])
#my_genes <- as.data.frame(FPKMs_logFC_bias_catergory$gene_id[FPKMs_logFC_bias_catergory$bias=="male_limited_exp"])
my_genes <- as.data.frame(FPKMs_logFC_bias_catergory$gene_id[FPKMs_logFC_bias_catergory$bias=="female_limited_exp"])

my_genes <- as.data.frame(na.omit(my_genes[,1]))
colnames(my_genes) <- "genes"
my_genes <- as.vector(my_genes[,1])

#-----------------------------------------------
# Keep only genes with annotated GOs

my_genes <- my_genes[my_genes %in% universe]
# male biased genes: 1975/4756
# female biased genes: 1911/5270
# extreme male biased genes: 24/140
# extreme female biased genes: 4/28
# limited male genes: 29/204
# limited female genes: 43/150

 

#-----------------------------------------------
# Set up paramters for hypergeometric test

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = my_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      param_list <- c(param_list,parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.05)


#-----------------------------------------------
# Hypergeometric test

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}


GO_enrichment <- Hyper_G_test(param_list = param_list)


#-----------------------------------------------
# Choose what you want to look for, here the choice is biological process over-represented 

Result <- summary(GO_enrichment[["BP_over"]])

#-----------------------------------------------
# FDR correction

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#-----------------------------------------------
# Make an output for REVIGO and write out

REVIGO <- Result_FDR[,1:2]

#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_male_biased_exp.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_female_biased_exp.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_ex_male_biased_exp.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_ex_female_biased_exp.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_male_limited_exp.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_female_limited_exp.txt",row.names = F,sep = "\t",quote = F)


