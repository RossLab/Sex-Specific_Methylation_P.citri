#-----------------------------------------------
# Script created by Alun Jones, see paper Bebane et al. (2019) Neonics and bumblebees...
#-----------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses")
library(GOstats)
library(GSEABase)
library(treemap)

#-----------------------------------------------
# Read in background GO set and make compatible with GOstats

GO_annotations <- read.table("./P_citri_GO_terms.txt")
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

#my_genes <- read.csv("./gene_lists/exon_female_hypermeth_genes.txt", header=F)
#my_genes <- read.csv("./gene_lists/exon_male_hypermeth_genes.txt", header=F)
#my_genes <- read.csv("./gene_lists/prom_female_hypermeth_genes.txt", header=F)
my_genes <- read.csv("./gene_lists/prom_male_hypermeth_genes.txt", header=F)


my_genes <- as.data.frame(na.omit(my_genes$V1))
colnames(my_genes) <- "genes"
my_genes <- as.vector(my_genes[,1])

#-----------------------------------------------
# Keep only genes with annotated GOs

my_genes <- my_genes[my_genes %in% universe]
# exon female hypermeth: 585/2706
# exon male hypermeth: 62/211
# prom female hypermeth: 464/2650
# prom male hypermeth: 69/236

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

#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_exon_female_hypermeth.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_exon_male_hypermeth.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_prom_female_hypermeth.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"./enriched_GO_terms/enriched_GOs_prom_male_hypermeth.txt",row.names = F,sep = "\t",quote = F)
