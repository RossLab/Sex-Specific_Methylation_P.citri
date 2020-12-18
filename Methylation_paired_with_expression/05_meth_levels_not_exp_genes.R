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
meth_exp <- meth_exp[!is.na(meth_exp$weightedMeth.mean),]
meth_exp[is.na(meth_exp)] <- 0

head(meth_exp)
colnames(meth_exp)[7]<-"weightedMeth"
# -----------------------------------------------
not_exp_only <- meth_exp[meth_exp$bias=="0",]
not_exp_only_proms <- not_exp_only[not_exp_only$feature=="promotors_2000bp",]
not_exp_only_exons <- not_exp_only[not_exp_only$feature=="exon_first3",]

ggplot(not_exp_only_exons, aes(x=weightedMeth, fill = origin))+
  geom_histogram()+
  ggtitle("Exons 1-3")+
  xlab("Weighted Methylation Level")+
  ylab("Number of Genes")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

ggplot(not_exp_only_proms, aes(x=weightedMeth, fill = origin))+
  geom_histogram()+
  ggtitle("Promoters")+
  xlab("Weighted Methylation Level")+
  ylab("Number of Genes")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

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
meth_exp$diff_exp_cat[meth_exp$bias==0] <- "not_exp"

prom_data <- meth_exp[meth_exp$feature=="promotors_2000bp",]
exon_data <- meth_exp[meth_exp$feature=="exon_first3",]

summary_exon<-summarySE(exon_data, measurevar = "weightedMeth", 
                        groupvars = c("diff_exp_cat","origin"))

ggplot(summary_exon, aes(x=diff_exp_cat, y=weightedMeth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth-ci, ymax=weightedMeth+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Exons 1-3")+
  xlab("Gene Group")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("not_exp","not_diff","diff"),
                   labels=c("Not\n Expressed","Not\nDifferentially\nExpressed", "Differentially\nExpressed"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

model1<-lm(weightedMeth ~ origin * diff_exp_cat, data=exon_data)
model2<-lm(weightedMeth ~ origin + diff_exp_cat, data=exon_data)
anova(model1,model2) # Sig interaction
summary.lm(model1) # only non sig: interaction with male diff exp is not diff to non-diff exp
exon_data$SHD<-interaction(exon_data$origin,exon_data$diff_exp_cat)
model1_new<-lm(weightedMeth~-1+SHD, data=exon_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) # only non sig: male.not_diff - male.diff 

# -----------------------------------------------

summary_proms<-summarySE(prom_data, measurevar = "weightedMeth", 
                         groupvars = c("diff_exp_cat","origin"))

ggplot(summary_proms, aes(x=diff_exp_cat, y=weightedMeth, fill=origin))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weightedMeth-ci, ymax=weightedMeth+ci),
                width=.2,
                position = position_dodge(.9))+
  ggtitle("Promotors")+
  xlab("Gene Group")+
  ylab("Average Weighted Methylation")+
  scale_fill_manual("", breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("not_exp","not_diff","diff"),
                   labels=c("Not\n Expressed","Not\nDifferentially\nExpressed", "Differentially\nExpressed"))+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        plot.title = element_text(size=20))

model1<-lm(weightedMeth ~ origin * diff_exp_cat, data=prom_data)
model2<-lm(weightedMeth ~ origin + diff_exp_cat, data=prom_data)
anova(model1,model2) # Sig interaction
summary.lm(model1) # only non sig: interaction with male diff exp is not diff to non-diff exp
prom_data$SHD<-interaction(prom_data$origin,prom_data$diff_exp_cat)
model1_new<-lm(weightedMeth~-1+SHD, data=prom_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey")))# only non sig: male.not_diff - male.diff 


