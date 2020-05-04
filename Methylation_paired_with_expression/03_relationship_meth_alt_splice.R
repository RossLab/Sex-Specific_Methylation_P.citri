# -----------------------------------------------
# Relationship methylation and alternative splicing
# -----------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/meth_paired_with_exp")

library(readr)
library(ggplot2)
library(FSA)
library(reshape2)
library(Hmisc)
library(tidyr)
library(multcomp)
library(UpSetR)
library(grid)

# -----------------------------------------------
# Read in all data
# -----------------------------------------------

# Methylation
methylation_all_data <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/weightedMeth_exons_promotors_only.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)

methylation_all_data_1 <- dcast(methylation_all_data, feature + gene_id ~ origin, value.var= "weightedMeth.mean")
#methylation_all_data_1$meth_diff <- ((methylation_all_data_1$male -
#                                       methylation_all_data_1$female) / methylation_all_data_1$male)*100
methylation_all_data_1$meth_diff <- methylation_all_data_1$female - methylation_all_data_1$male


methylation_all_data <- melt(methylation_all_data_1, id.vars=c("feature", "gene_id","meth_diff"))
colnames(methylation_all_data) <- c("feature","gene_id","meth_diff","origin","weighted_meth")

diff_meth_exon_geneIDs <- read_csv("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/diff_meth_gene_lists/diff_meth_exon_geneIDs.txt", 
                                   col_names = FALSE)
diff_meth_promotor_geneIDs <- read_csv("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/diff_meth_gene_lists/diff_meth_promotor_geneIDs.txt", 
                                       col_names = FALSE)

colnames(diff_meth_exon_geneIDs) <- "gene_id"
colnames(diff_meth_promotor_geneIDs) <- "gene_id"

diff_meth_exon_geneIDs$diff_meth_feature <- "exon"
diff_meth_promotor_geneIDs$diff_meth_feature <- "promotor"

diff_meth_genes <- rbind(diff_meth_exon_geneIDs,diff_meth_promotor_geneIDs)

diff_meth_genes_all_info <- merge(methylation_all_data, diff_meth_genes, all=T)
diff_meth_genes_all_info$diff_meth_feature[is.na(diff_meth_genes_all_info$diff_meth_feature)] <- "not_diff_meth"


diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "exon" & diff_meth_genes_all_info$meth_diff < 0)] <- "exon_male_meth"
diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "exon" & diff_meth_genes_all_info$meth_diff > 0)] <- "exon_female_meth"
diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "promotor" & diff_meth_genes_all_info$meth_diff < 0)] <- "prom_male_meth"
diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "promotor" & diff_meth_genes_all_info$meth_diff > 0)] <- "prom_female_meth"

diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "promotor" & diff_meth_genes_all_info$meth_diff == 0)] <-"not_diff_meth"
diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "exon" & diff_meth_genes_all_info$meth_diff == 0)] <-"not_diff_meth"
diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "promotor" & is.na(diff_meth_genes_all_info$meth_diff))] <-"not_diff_meth"
diff_meth_genes_all_info$diff_meth_feature[(diff_meth_genes_all_info$diff_meth_feature ==
                                              "exon" & is.na(diff_meth_genes_all_info$meth_diff))] <-"not_diff_meth"


# Alternative splicing data 209 genes which use diff isoforms between sexes
alt_splice_with_exp_bias <- read_delim("../transcription/alt_splice_with_exp_bias.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
# do it's fair, non-alt spliced genes must be present in the RNA-seq data
FPKMs_logFC_bias_catergory <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/FPKMs_logFC_bias_catergory.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)

head(FPKMs_logFC_bias_catergory)
FPKMs_logFC_bias_catergory$alt_spliced <- "no"
FPKMs_logFC_bias_catergory$alt_spliced[FPKMs_logFC_bias_catergory$gene_id %in% alt_splice_with_exp_bias$gene_id] <- "yes"
table(FPKMs_logFC_bias_catergory$alt_spliced)
FPKMs_logFC_bias_catergory <- FPKMs_logFC_bias_catergory[,c(1,5,6)]


# Q) what are the methylation levels of these genes in males and females compared to non-switching genes?
all_data <- merge(diff_meth_genes_all_info, FPKMs_logFC_bias_catergory, by=c("gene_id"), all=T)
# Remove alt spliced genes not in the meth data
all_data <- all_data[!is.na(all_data$weighted_meth),] 
all_data <- all_data[!is.na(all_data$alt_spliced),] 

head(all_data)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Relationship with exon methylation
exon_data <- all_data[all_data$feature=="exon_first3",]
#range(exon_data$weighted_meth[exon_data$alt_spliced=="yes"]) #0.001076106 0.153879087
#all_data_low <- all_data[all_data$weighted_meth < 0.2,]

#exon_data_low <- all_data_low[all_data_low$feature=="exon_first3",]

summary_exon<-summarySE(exon_data, measurevar = "weighted_meth", 
                        groupvars = c("alt_spliced","origin"))
unique(exon_data$alt_spliced)

ggplot(summary_exon, aes(x=alt_spliced, y=weighted_meth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-se, ymax=weighted_meth+se),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Exons 1-3")+
  xlab("Gene Group")+
  ylab("Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("no", "yes"),
                   labels=c("Not Alternatively\nSpliced","Alternatively\nSpliced"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=22),
        plot.title = element_text(size=22))

# Stats
model1<-lm(weighted_meth ~ origin * alt_spliced , data=exon_data)
model2<-lm(weighted_meth ~ origin + alt_spliced , data=exon_data)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Sex sig, alt splice status is sig, no interaction effect

exon_data$SHD<-interaction(exon_data$origin,exon_data$alt_spliced)
model1_new<-lm(weighted_meth~-1+SHD, data=exon_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) 
# Everything is different to everything except :
# male.yes - male.no == 0     -0.014502   0.006861  -2.114  0.12106  
# female.yes - female.no == 0 -0.023972   0.006861  -3.494  0.00197 **  



# Relationship with promotor methylation
prom_data <- all_data[all_data$feature=="promotors_2000bp",]
#range(prom_data$weighted_meth[prom_data$alt_spliced=="yes"]) #0.000000 0.568094

#all_data_low <- all_data[all_data$weighted_meth < 0.6,]

#prom_data_low <- all_data_low[all_data_low$feature=="promotors_2000bp",]

summary_prom<-summarySE(prom_data, measurevar = "weighted_meth", 
                        groupvars = c("alt_spliced","origin"))

ggplot(summary_prom, aes(x=alt_spliced, y=weighted_meth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-se, ymax=weighted_meth+se),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Promotors")+
  xlab("Gene Group")+
  ylab("Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("no", "yes"),
                   labels=c("Not Alternatively\nSpliced","Alternatively\nSpliced"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=22),
        plot.title = element_text(size=22))

# Stats
model1<-lm(weighted_meth ~ origin * alt_spliced , data=prom_data)
model2<-lm(weighted_meth ~ origin + alt_spliced , data=prom_data)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Sex sig, alt splice status is sig

prom_data$SHD<-interaction(prom_data$origin,prom_data$alt_spliced)
model1_new<-lm(weighted_meth~-1+SHD, data=prom_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) 
# as with exons:
# female.yes - female.no == 0 -0.024217   0.006993  -3.463  0.00218 **  
# male.yes - male.no == 0     -0.013108   0.007012  -1.869  0.20319  







exon_data <- all_data[all_data$feature=="exon_first3",]
exon_data$bias[exon_data$alt_spliced=="no"]<-"not_spliced"

summary_exon<-summarySE(exon_data, measurevar = "weighted_meth", 
                        groupvars = c("bias","origin"))
unique(exon_data$bias)

ggplot(summary_exon, aes(x=bias, y=weighted_meth, fill=origin))+
geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-se, ymax=weighted_meth+se),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Exons 1-3")+
  xlab("Gene Group")+
  ylab("Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("not_spliced","Unbias",
                            "Mbias","Fbias"),
                   labels=c("Not AS",
                            "AS\nUnbiased","AS\nMale Biased","AS\nFemale Biased"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=20,hjust=-0.01),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=22),
        plot.title = element_text(size=22))

# Stats
model1<-lm(weighted_meth ~ origin * bias , data=exon_data)
model2<-lm(weighted_meth ~ origin + bias , data=exon_data)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Sex sig, alt splice status is sig, no interaction effect

exon_data$SHD<-interaction(exon_data$origin,exon_data$bias)
model1_new<-lm(weighted_meth~-1+SHD, data=exon_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) 



# Relationship with promotor methylation
prom_data <- all_data[all_data$feature=="promotors_2000bp",]
prom_data$bias[prom_data$alt_spliced=="no"]<-"not_spliced"

summary_prom<-summarySE(prom_data, measurevar = "weighted_meth", 
                        groupvars = c("bias","origin"))

ggplot(summary_prom, aes(x=bias, y=weighted_meth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-se, ymax=weighted_meth+se),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Promotors")+
  xlab("Gene Group")+
  ylab("Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("not_spliced","Unbias",
                            "Mbias","Fbias"),
                   labels=c("Not AS",
                            "AS\nUnbiased","AS\nMale Biased","AS\nFemale Biased"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=20,hjust=-0.01),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=22),
        plot.title = element_text(size=22))

# Stats
model1<-lm(weighted_meth ~ origin * bias , data=prom_data)
model2<-lm(weighted_meth ~ origin + bias , data=prom_data)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Sex not sig, alt splice status is sig

prom_data$SHD<-interaction(prom_data$origin,prom_data$bias)
model1_new<-lm(weighted_meth~-1+SHD, data=prom_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) 
# Everything is different to everything except female/male alt splice is not diff




#Significant overlap with diff meth???

head(prom_data)
upset_data <- prom_data[,c(1,6,7)]

prom_female_hyper <- subset(upset_data, diff_meth_feature == "prom_female_meth")
prom_male_hyper <- subset(upset_data, diff_meth_feature == "prom_male_meth")

# These are all alt spliced
Mbias <- subset(upset_data, bias == "Mbias")
Fbias <- subset(upset_data, bias == "Fbias")
Unbias <- subset(upset_data, bias == "Unbias")


upset_data1 <- as.data.frame(unique(prom_data$gene_id))
colnames(upset_data1) <- "gene_id"

upset_data1$prom_female_hyper <- 0
upset_data1$prom_female_hyper[upset_data1$gene_id %in% prom_female_hyper$gene_id] <- 1
upset_data1$prom_male_hyper <- 0
upset_data1$prom_male_hyper[upset_data1$gene_id %in% prom_male_hyper$gene_id] <- 1
upset_data1$Mbias <- 0
upset_data1$Mbias[upset_data1$gene_id %in% Mbias$gene_id] <- 1
upset_data1$Fbias <- 0
upset_data1$Fbias[upset_data1$gene_id %in% Fbias$gene_id] <- 1
upset_data1$Unbias <- 0
upset_data1$Unbias[upset_data1$gene_id %in% Unbias$gene_id] <- 1

row.names(upset_data1) <- upset_data1$gene_id
upset_data1 <- upset_data1[,-1]

colnames(upset_data1) <- c("Female Hypermethylated","Male Hypermethylated",
                         "Male Biased","Female Biased","Unbiased")

upset(upset_data1, nsets =5, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity",
      mainbar.y.label =NULL)
grid.text("Alternatively Spliced",x = 0.15, y=0.85, gp=gpar(fontsize=20))
grid.text("Intersection Size",x = 0.35, y=0.60, gp=gpar(fontsize=20), rot = 90)

# number of overlap, number of one category, number of all genes - other category, other category
# Genes with exp and meth data = 15341

# Biased genes proms
sum(upset_data1$`Female Hypermethylated`) #293
sum(upset_data1$`Unbiased`) #72
a1 <- phyper(1, 72, 15048, 293, lower.tail = F) # p = 0.4081699
sum(upset_data1$`Male Biased`) #99
b1 <- phyper(1, 99, 15048, 293, lower.tail = F) # p = 0.5738667

sum(upset_data1$`Male Hypermethylated`) #30
c1 <- phyper(1, 99, 15317, 30, lower.tail = F) # p =0.01579959***




head(exon_data)
upset_data <- exon_data[,c(1,6,7)]

exon_female_hyper <- subset(upset_data, diff_meth_feature == "exon_female_hyper")
exon_male_hyper <- subset(upset_data, diff_meth_feature == "exon_male_meth")

# These are all alt spliced
Mbias <- subset(upset_data, bias == "Mbias")
Fbias <- subset(upset_data, bias == "Fbias")
Unbias <- subset(upset_data, bias == "Unbias")


upset_data1 <- as.data.frame(unique(prom_data$gene_id))
colnames(upset_data1) <- "gene_id"

upset_data1$exon_female_hyper <- 0
upset_data1$exon_female_hyper[upset_data1$gene_id %in% exon_female_hyper$gene_id] <- 1
upset_data1$exon_male_hyper <- 0
upset_data1$exon_male_hyper[upset_data1$gene_id %in% exon_male_hyper$gene_id] <- 1
upset_data1$Mbias <- 0
upset_data1$Mbias[upset_data1$gene_id %in% Mbias$gene_id] <- 1
upset_data1$Fbias <- 0
upset_data1$Fbias[upset_data1$gene_id %in% Fbias$gene_id] <- 1
upset_data1$Unbias <- 0
upset_data1$Unbias[upset_data1$gene_id %in% Unbias$gene_id] <- 1

row.names(upset_data1) <- upset_data1$gene_id
upset_data1 <- upset_data1[,-1]

colnames(upset_data1) <- c("Female Hypermethylated","Male Hypermethylated",
                           "Male Biased","Female Biased","Unbiased")

upset(upset_data1, nsets =5, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity",
      mainbar.y.label =NULL)
grid.text("Alternatively Spliced",x = 0.15, y=0.85, gp=gpar(fontsize=20))
grid.text("Intersection Size",x = 0.35, y=0.60, gp=gpar(fontsize=20), rot = 90)

# Biased genes exons
sum(upset_data1$`Male Hypermethylated`) #30
sum(upset_data1$`Unbiased`) #72
a2 <- phyper(1, 72, 15317, 30, lower.tail = F) # p = 0.008628158**

x <- c(a1,b1,c1,a2)
p.adjust(x, method="fdr", n=length(x)) 

#0.54422654 0.57386673 0.03159918 0.03159918




