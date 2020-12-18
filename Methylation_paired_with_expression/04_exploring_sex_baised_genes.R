# Make graphs etc for genes with ONLY male or female expression

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/meth_paired_with_exp")

library(readr)
library(reshape2)
library(Hmisc)
library(dplyr)
library(tidyr)
library(FSA)
library(ggplot2)
library(multcomp)

# -----------------------------------------------
# Read in data
# -----------------------------------------------

weightedMeth_exons_promotors_only <- read_delim("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation/weightedMeth_exons_promotors_only.txt", 
                                                "\t", escape_double = FALSE, trim_ws = TRUE)

FPKMs_logFC_bias_catergory <- read_delim("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/transcription/FPKMs_logFC_bias_catergory.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
FPKM_values_by_sex <- melt(FPKMs_logFC_bias_catergory, id.vars=c("gene_id","log2FC","bias"))
colnames(FPKM_values_by_sex)[4]<-"origin"
colnames(FPKM_values_by_sex)[5]<-"FPKM"

meth_exp <- merge(FPKM_values_by_sex, weightedMeth_exons_promotors_only, by=c("gene_id","origin"), all=T)
meth_exp <- meth_exp[!is.na(meth_exp$weightedMeth.mean) & (!is.na(meth_exp$FPKM)),]
meth_exp$logFPKM <- log(meth_exp$FPKM)
meth_exp$logFPKM[meth_exp$logFPKM == -Inf] <- 0


head(meth_exp)

# Low 0.5-10, medium 11-1000, high > 1000 https://www.ebi.ac.uk/gxa/FAQ.html
#hist(meth_exp$FPKM[(meth_exp$FPKM >0 & meth_exp$FPKM <0.5)])
#meth_exp$exp_status <- "Low"
#meth_exp$exp_status[meth_exp$FPKM == 0] <- "Not Expressed"
#meth_exp$exp_status[meth_exp$FPKM >10 & meth_exp$FPKM <1000] <- "Medium"
#meth_exp$exp_status[meth_exp$FPKM >1000] <- "High"

# -----------------------------------------------
# Summary function
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

# -----------------------------------------------
# Bar plots of a gross overview
# -----------------------------------------------
meth_exp$diff_exp_cat <- "diff"
meth_exp$diff_exp_cat[meth_exp$bias=="Unbias"] <- "not_diff"

#range(meth_exp$weightedMeth.mean[meth_exp$diff_exp_cat=="diff"])

prom_data <- meth_exp[meth_exp$feature=="promotors_2000bp",]
exon_data <- meth_exp[meth_exp$feature=="exon_first3",]

summary_exon<-summarySE(exon_data, measurevar = "weightedMeth.mean", 
                        groupvars = c("diff_exp_cat","origin"))

ggplot(summary_exon, aes(x=diff_exp_cat, y=weightedMeth.mean, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth.mean-ci, ymax=weightedMeth.mean+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Exons 1-3")+
  xlab("Gene Group")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("not_diff","diff"),
                   labels=c("Not\nDifferentially\nExpressed", "Differentially\nExpressed"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

model1<-lm(weightedMeth.mean ~ origin * diff_exp_cat, data=exon_data)
model2<-lm(weightedMeth.mean ~ origin + diff_exp_cat, data=exon_data)
anova(model1,model2) # Sig interaction
summary.lm(model1) # All sig
exon_data$SHD<-interaction(exon_data$origin,exon_data$diff_exp_cat)
model1_new<-lm(weightedMeth.mean~-1+SHD, data=exon_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) # all dig

summary_proms<-summarySE(prom_data, measurevar = "weightedMeth.mean", 
                         groupvars = c("diff_exp_cat","origin"))

ggplot(summary_proms, aes(x=diff_exp_cat, y=weightedMeth.mean, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth.mean-ci, ymax=weightedMeth.mean+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Promotors")+
  xlab("Gene Group")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("not_diff","diff"),
                   labels=c("Not\nDifferentially\nExpressed", "Differentially\nExpressed"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

model1<-lm(weightedMeth.mean ~ origin * diff_exp_cat, data=prom_data)
model2<-lm(weightedMeth.mean ~ origin + diff_exp_cat, data=prom_data)
anova(model1,model2) # Sig interaction
summary.lm(model1) # Low and no-meth different to others, sex not sig
prom_data$SHD<-interaction(prom_data$origin,prom_data$diff_exp_cat)
model1_new<-lm(weightedMeth.mean~-1+SHD, data=prom_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey")))

# -----------------------------------------------
# Bar plots of it all
# -----------------------------------------------

prom_data <- meth_exp[meth_exp$feature=="promotors_2000bp",]
exon_data <- meth_exp[meth_exp$feature=="exon_first3",]

summary_exon<-summarySE(exon_data, measurevar = "weightedMeth.mean", 
                        groupvars = c("bias","origin"))

ggplot(summary_exon, aes(x=bias, y=weightedMeth.mean, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth.mean-ci, ymax=weightedMeth.mean+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Exons 1-3")+
  xlab("Gene Expression Level")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("Unbias", "Fbias",
                            "Mbias","extreme_female_biased","extreme_male_biased",
                            "female_limited_exp","male_limited_exp"),
                   labels=c("Unbiased","Female\nBiased","Male\nBiased",
                            "Extreme\nFemale\nBiased","Extreme\nMale\nBiased",
                            "Female\nLimited","Male\nLimited"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

model1<-lm(weightedMeth.mean ~ origin * bias, data=exon_data)
model2<-lm(weightedMeth.mean ~ origin + bias, data=exon_data)
anova(model1,model2) # Sig interaction
summary.lm(model1) # Low and no-meth different to others, sex not sig
exon_data$SHD<-interaction(exon_data$origin,exon_data$bias)
model1_new<-lm(weightedMeth.mean~-1+SHD, data=exon_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey")))

summary_proms<-summarySE(prom_data, measurevar = "weightedMeth.mean", 
                        groupvars = c("bias","origin"))

ggplot(summary_proms, aes(x=bias, y=weightedMeth.mean, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth.mean-ci, ymax=weightedMeth.mean+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Promotors")+
  xlab("Gene Expression Level")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("Unbias", "Fbias",
                            "Mbias","extreme_female_biased","extreme_male_biased",
                            "female_limited_exp","male_limited_exp"),
                   labels=c("Unbiased","Female\nBiased","Male\nBiased",
                            "Extreme\nFemale\nBiased","Extreme\nMale\nBiased",
                            "Female\nLimited","Male\nLimited"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

model1<-lm(weightedMeth.mean ~ origin * bias, data=prom_data)
model2<-lm(weightedMeth.mean ~ origin + bias, data=prom_data)
anova(model1,model2) # Sig interaction
summary.lm(model1)
prom_data$SHD<-interaction(prom_data$origin,prom_data$bias)
model1_new<-lm(weightedMeth.mean~-1+SHD, data=prom_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey")))












# -----------------------------------------------
# Plot of meth of just sex-biased genes (IGNORE THE BELOW)
# -----------------------------------------------

head(meth_exp)

female_notexp_genes <- meth_exp[meth_exp$origin=="female" & meth_exp$FPKM == 0,]
male_notexp_genes <- meth_exp[meth_exp$origin=="male" & meth_exp$FPKM == 0,]

female_not_exp_geneIDs <- unique(female_notexp_genes$gene_id)
male_not_exp_geneIDs <- unique(male_notexp_genes$gene_id)

female_not_exp_geneIDs_unique <- female_not_exp_geneIDs[!female_not_exp_geneIDs %in% male_not_exp_geneIDs]
male_not_exp_geneIDs_unique <- male_not_exp_geneIDs[!male_not_exp_geneIDs %in% female_not_exp_geneIDs]


all_data_female_not_exp <- meth_exp[meth_exp$gene_id %in% female_not_exp_geneIDs_unique,]
all_data_male_not_exp <- meth_exp[meth_exp$gene_id %in% male_not_exp_geneIDs_unique,]

all_data_female_not_exp$bias_status <- "female_not_exp"
all_data_male_not_exp$bias_status <- "male_not_exp"

all_data <- rbind(all_data_female_not_exp, all_data_male_not_exp)
all_prom_data <- all_data[all_data$feature=="promotors_2000bp",]
all_exon_data <- all_data[all_data$feature=="exon_first3",]

summary_prom_bias_genes<-summarySE(all_prom_data, measurevar = "weightedMeth.mean", 
                        groupvars = c("bias_status","origin"))

ggplot(summary_prom_bias_genes, aes(x=bias_status, y=weightedMeth.mean, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth.mean-ci, ymax=weightedMeth.mean+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Promotors")+
  xlab("Biased Gene Category")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_not_exp", "male_not_exp"),
                   labels=c("Silent in\nFemales","Silent in\nMales"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        title = element_text(size=16))

summary_exon_bias_genes<-summarySE(all_exon_data, measurevar = "weightedMeth.mean", 
                                   groupvars = c("bias_status","origin"))

ggplot(summary_exon_bias_genes, aes(x=bias_status, y=weightedMeth.mean, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth.mean-ci, ymax=weightedMeth.mean+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Exons 1-3")+
  xlab("Biased Gene Category")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_not_exp", "male_not_exp"),
                   labels=c("Silent in\nFemales","Silent in\nMales"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        title = element_text(size=16))
