# -----------------------------------------------
# Reltionship of differential methylation and Diff exp (using prom & exon as cetegory)
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

# -----------------------------------------------
# Read in all data
# -----------------------------------------------

# Methylation
methylation_all_data <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/weightedMeth_exons_promotors_only.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)

methylation_all_data_1 <- dcast(methylation_all_data, feature + gene_id ~ origin, value.var= "weightedMeth.mean")
methylation_all_data_1$meth_diff <- methylation_all_data_1$female - methylation_all_data_1$male

methylation_all_data <- melt(methylation_all_data_1, id.vars=c("feature", "gene_id","meth_diff"))
colnames(methylation_all_data) <- c("feature","gene_id","meth_diff","origin","weighted_meth")


diff_meth_common_promotor_exon_geneIDs <- read_csv("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/diff_meth_gene_lists/diff_meth_common_promotor_exon_geneIDs.txt", 
                                                   col_names = FALSE)
diff_meth_unique_exon_geneIDs <- read_csv("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/diff_meth_gene_lists/diff_meth_unique_exon_geneIDs.txt", 
                                          col_names = FALSE)
diff_meth_unique_promotor_geneIDs <- read_csv("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/diff_meth_gene_lists/diff_meth_unique_promotor_geneIDs.txt", 
                                              col_names = FALSE)

colnames(diff_meth_common_promotor_exon_geneIDs) <- "gene_id"
colnames(diff_meth_unique_exon_geneIDs) <- "gene_id"
colnames(diff_meth_unique_promotor_geneIDs) <- "gene_id"

diff_meth_common_promotor_exon_geneIDs$diff_meth_feature <- "prom_and_exon"
diff_meth_unique_exon_geneIDs$diff_meth_feature <- "exon"
diff_meth_unique_promotor_geneIDs$diff_meth_feature <- "promotor"

diff_meth_genes <- rbind(diff_meth_common_promotor_exon_geneIDs,diff_meth_unique_exon_geneIDs,diff_meth_unique_promotor_geneIDs)

diff_meth_genes_all_info <- merge(methylation_all_data, diff_meth_genes, all=T)
diff_meth_genes_all_info$diff_meth_feature[is.na(diff_meth_genes_all_info$diff_meth_feature)] <- "not_diff_meth"

# Expression
logFC_DEgenes <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/logFC_DEgenes.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
logFC_DEgenes <- logFC_DEgenes[,c(1,8,9)]
colnames(logFC_DEgenes) <- c("gene_id","logFC","bias_status")

FPKM_values_by_sex <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/FPKM_values_by_sex.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
diff_exp_data <- merge(logFC_DEgenes, FPKM_values_by_sex, by="gene_id")
diff_exp_data$logFPKM <- log(diff_exp_data$FPKM)


# Put them together
diffmeth_diffexp <- merge(diff_meth_genes_all_info, diff_exp_data, by=c("gene_id","origin"))
diffmeth_diffexp[is.na(diffmeth_diffexp)] <- 0

# NOTE: we have methylation data on many more genes compared to expression data
# This is why we don't see as many of the diff methylated genes on these expression
# comparison plots
length(unique(diff_meth_genes_all_info$gene_id)) #35782
length(unique(diff_exp_data$gene_id)) #15570
length(unique(diffmeth_diffexp$gene_id)) #15341


#write.table(diffmeth_diffexp, file="diff_meth_and_exp_together.txt", sep='\t',
#           col.names = T, row.names = F, quote = F)

# -----------------------------------------------
# Scatter plots of all differential data
# -----------------------------------------------

#Lets take a look, scatter plot all data
ggplot(diffmeth_diffexp, aes(x=meth_diff, y=logFC, colour=feature  ))+
  geom_point()+
  xlab("Weighted Methylation Difference")+
  ylab("log(Fold-Change)")+
  scale_colour_manual("",breaks=c("exon_first3", "promotors_2000bp"),
                      values = c("#9933FF","#00CC66"),
                      labels= c("Exons 1-3", "Promotors"))+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_blank())+
  geom_vline(xintercept=0.15)+
  geom_vline(xintercept=-0.15)+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  xlim(-0.5,0.6)


# Break it down to exons and promotors
exon_data <- subset(diffmeth_diffexp, feature == "exon_first3")
exon_data$diff_meth_feature[exon_data$diff_meth_feature=="promotor"]<-"not_diff_meth"

promotor_data <- subset(diffmeth_diffexp, feature == "promotors_2000bp")
promotor_data$diff_meth_feature[promotor_data$diff_meth_feature=="exon"]<-"not_diff_meth"

# As exon methylation has been averaged here for plotting some significant exon overall differences have been
# reduced to > 0.15, remove these to make plotting neater
exon_data$diff_meth_feature[exon_data$meth_diff < 0.15 & exon_data$diff_meth > -0.15] <-"not_diff_meth"
promotor_data$diff_meth_feature[promotor_data$meth_diff < 0.15 & promotor_data$diff_meth > -0.15] <-"not_diff_meth"


ggplot(exon_data, aes(x=meth_diff, y=logFC, colour=diff_meth_feature))+
  geom_point()+
  xlab("Weighted Methylation Difference")+
  ylab("log(Fold-Change)")+
  ggtitle("Exons 1-3")+
  scale_colour_manual("",
                      limits=c("exon","prom_and_exon" ,"not_diff_meth"),
                      values = c("#9933FF","#CCCC00","#999999"),
                      labels= c("Exon Only", "Exon and Promotor","Not Differentially\nMethylated"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        title = element_text(size=16))+
  geom_vline(xintercept=0.15)+
  geom_vline(xintercept=-0.15)+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  xlim(-0.5,0.6)


ggplot(promotor_data, aes(x=meth_diff, y=logFC, colour=diff_meth_feature))+
  geom_point()+
  xlab("Weighted Methylation Difference")+
  ylab("log(Fold-Change)")+
  ggtitle("Promotors")+
  scale_colour_manual("",
                      limits=c("promotor","prom_and_exon" ,"not_diff_meth"),
                      values = c("#00CC66","#CCCC00","#999999"),
                      labels= c("Promotor Only", "Exon and Promotor","Not Differentially\nMethylated"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        title = element_text(size=16))+
  geom_vline(xintercept=0.15)+
  geom_vline(xintercept=-0.15)+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  xlim(-0.5,0.6)



# -----------------------------------------------
# Binned diff/non-diff methylated gene, violin plots
# -----------------------------------------------

condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

head(diffmeth_diffexp)

# Urgh, if I want the mean as a red dot I need to make a new plotting column
diffmeth_diffexp$plot_column <- paste0(diffmeth_diffexp$origin,"_",diffmeth_diffexp$diff_meth_feature)

ggplot(diffmeth_diffexp, aes(x=plot_column, y=logFPKM)) + 
  geom_violin( aes(fill=origin))+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.02, dotsize = 0.3)+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  xlab("Gene Group")+
  ylab("log(FPKM)")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_not_diff_meth", "male_not_diff_meth","female_prom_and_exon",
                           "male_prom_and_exon","female_promotor" ,"male_promotor" ,"female_exon" ,"male_exon"),
                   labels=c("Non-DM", "Non-DM",
                            "Promotor\nand Exon ", "Promotor\nand Exon",
                            " Promotor\nOnly","Promotor\nOnly","Exon\nOnly","Exon\nOnly"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16))

# Stats
diffmeth_diffexp$diff_meth_feature <- as.factor(diffmeth_diffexp$diff_meth_feature)
model1<-lm(FPKM ~ origin * diff_meth_feature, data=diffmeth_diffexp)
model2<-lm(FPKM ~ origin + diff_meth_feature, data=diffmeth_diffexp)
anova(model1,model2) # No sig interaction
summary.lm(model2) 

summary(glht(model2, linfct=mcp(diff_meth_feature="Tukey")))
# sig higher exp in non-diff meth genes compared to those with diff meth exon
# and diff meth exons & proms but not compared to diff meth proms

# -----------------------------------------------
# Binned diff/non-diff exp gene, violin plots
# -----------------------------------------------

exon_data <- subset(diffmeth_diffexp, feature=="exon_first3")
head(exon_data)
exon_data$log_meth <- log(exon_data$weighted_meth)
exon_data$plot_column <- paste0(exon_data$origin,"_",exon_data$bias_status)

promotor_data <- subset(diffmeth_diffexp, feature=="promotors_2000bp")
head(promotor_data)
promotor_data$log_meth <- log(promotor_data$weighted_meth)
promotor_data$plot_column <- paste0(promotor_data$origin,"_",promotor_data$bias_status)


ggplot(exon_data, aes(x=plot_column, y=weighted_meth)) + 
  geom_violin( aes(fill=origin))+
  #geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.02, dotsize = 0.3)+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  ggtitle("Exons 1-3")+
  xlab("Gene Group")+
  ylab("log(Average Weighted Methylation)")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_Unbias", "male_Unbias","female_Fbias",
                            "male_Fbias","female_Mbias" ,"male_Mbias"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        title = element_text(size=16))


ggplot(promotor_data, aes(x=plot_column, y=weighted_meth)) + 
  geom_violin( aes(fill=origin))+
  #geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.02, dotsize = 0.3)+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  ggtitle("Promotors")+
  xlab("Gene Group")+
  ylab("log(Average Weighted Methylation)")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_Unbias", "male_Unbias","female_Fbias",
                            "male_Fbias","female_Mbias" ,"male_Mbias"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        title = element_text(size=16))


# -----------------------------------------------
# Same binned diff exp but as bar plots
# -----------------------------------------------

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

summary_exon<-summarySE(exon_data, measurevar = "weighted_meth", 
                        groupvars = c("plot_column","origin"))

ggplot(summary_exon, aes(x=plot_column, y=weighted_meth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-ci, ymax=weighted_meth+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Exons 1-3")+
  xlab("Gene Group")+
  ylab("log(Average Weighted Methylation)")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_Unbias", "male_Unbias","female_Fbias",
                            "male_Fbias","female_Mbias" ,"male_Mbias"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        title = element_text(size=16))



summary_proms<-summarySE(promotor_data, measurevar = "weighted_meth", 
                        groupvars = c("plot_column","origin"))

ggplot(summary_proms, aes(x=plot_column, y=weighted_meth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-ci, ymax=weighted_meth+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Promotors")+
  xlab("Gene Group")+
  ylab("log(Average Weighted Methylation)")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_Unbias", "male_Unbias","female_Fbias",
                            "male_Fbias","female_Mbias" ,"male_Mbias"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        title = element_text(size=16))

# the two above graphs are basically the same, so pool prom and exon data
diffmeth_diffexp$plot_column <- paste0(diffmeth_diffexp$origin,"_",diffmeth_diffexp$bias_status)

summary_all<-summarySE(diffmeth_diffexp, measurevar = "weighted_meth", 
                         groupvars = c("plot_column","origin"))

ggplot(summary_all, aes(x=plot_column, y=weighted_meth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-ci, ymax=weighted_meth+ci),
                width=.2,
                position = position_dodge(.9))+
  xlab("Gene Group")+
  ylab("log(Average Weighted Methylation)")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_Unbias", "male_Unbias","female_Fbias",
                            "male_Fbias","female_Mbias" ,"male_Mbias"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16))

# Stats
exon_data$bias_status <- as.factor(exon_data$bias_status)
model1<-lm(weighted_meth ~ origin * bias_status, data=exon_data)
model2<-lm(weighted_meth ~ origin + bias_status, data=exon_data)
anova(model1,model2) # Sig interaction effect
summary.lm(model1) 

exon_data$SHD<-interaction(exon_data$origin,exon_data$bias_status)
model1_new<-lm(weighted_meth~-1+SHD, data=exon_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) # all comparisons


promotor_data$bias_status <- as.factor(promotor_data$bias_status)
model1<-lm(weighted_meth ~ origin * bias_status, data=promotor_data)
model2<-lm(weighted_meth ~ origin + bias_status, data=promotor_data)
anova(model1,model2) # Sig interaction effect
summary.lm(model1) 

promotor_data$SHD<-interaction(promotor_data$origin,promotor_data$bias_status)
model1_new<-lm(weighted_meth~-1+SHD, data=promotor_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey")))

# The results look the same so can just pool exon and promotor as for graph
diffmeth_diffexp$bias_status <- as.factor(diffmeth_diffexp$bias_status)
model1<-lm(weighted_meth ~ origin * bias_status, data=diffmeth_diffexp)
model2<-lm(weighted_meth ~ origin + bias_status, data=diffmeth_diffexp)
anova(model1,model2) # Sig interaction effect
summary.lm(model1) 

diffmeth_diffexp$SHD<-interaction(diffmeth_diffexp$origin,diffmeth_diffexp$bias_status)
model1_new<-lm(weighted_meth~-1+SHD, data=diffmeth_diffexp)
summary(glht(model1_new,linfct=mcp(SHD="Tukey")))
# Everything sig diff to everything else, expect male.unbias is not sig compared
# to male.Fbias

# -----------------------------------------------
# Number of genes diff meth and diff exp
# -----------------------------------------------

# On second thought this is kind of pointless, really need to label up the proms and exons
# as being female methylated or male methylated and not worry if they occur in both
# exons and promotors for the same gene, this will give three categories as with the
# expression data, allowing much better comparisons to be made

head(diffmeth_diffexp)
upset_data <- diffmeth_diffexp[,c(1,6,8)]

# Urgh can't think of a short way of doing this
not_diff_meth <- subset(upset_data, diff_meth_feature == "not_diff_meth")
prom_and_exon <- subset(upset_data, diff_meth_feature == "prom_and_exon")
promotor <- subset(upset_data, diff_meth_feature == "promotor")
exon <- subset(upset_data, diff_meth_feature == "exon")

Unbias <- subset(upset_data, bias_status == "Unbias")
Mbias <- subset(upset_data, bias_status == "Mbias")
Fbias <- subset(upset_data, bias_status == "Fbias")

upset_data1 <- as.data.frame(unique(diffmeth_diffexp$gene_id))
colnames(upset_data1) <- "gene_id"

#upset_data1$not_diff_meth <- 0
#upset_data1$not_diff_meth[upset_data1$gene_id %in% not_diff_meth$gene_id] <- 1
upset_data1$prom_and_exon <- 0
upset_data1$prom_and_exon[upset_data1$gene_id %in% prom_and_exon$gene_id] <- 1
upset_data1$promotor <- 0
upset_data1$promotor[upset_data1$gene_id %in% promotor$gene_id] <- 1
upset_data1$exon <- 0
upset_data1$exon[upset_data1$gene_id %in% exon$gene_id] <- 1
upset_data1$Unbias <- 0
upset_data1$Unbias[upset_data1$gene_id %in% Unbias$gene_id] <- 1
upset_data1$Mbias <- 0
upset_data1$Mbias[upset_data1$gene_id %in% Mbias$gene_id] <- 1
upset_data1$Fbias <- 0
upset_data1$Fbias[upset_data1$gene_id %in% Fbias$gene_id] <- 1
row.names(upset_data1) <- upset_data1$gene_id
upset_data1 <- upset_data1[,-1]

upset(upset_data1, nsets = 6, order.by = "freq",
      text.scale = 2,
      point.size = 4,
      scale.sets = "identity")



