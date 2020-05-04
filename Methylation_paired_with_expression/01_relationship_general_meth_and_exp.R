# Making line and scatter graphs to show relationship of gene exp and meth

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/meth_paired_with_exp")

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
weighted_meth <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/PCITRI_weighted_meth_annotation_by_sex.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
#weighted_meth_wide <- spread(weighted_meth, origin, weightedMeth.mean)

#weighted_meth_wide$meth_diff <- weighted_meth_wide$female - weighted_meth_wide$male # + value = higher in female
weighted_meth$gene_id <- gsub("ID=","",weighted_meth$gene_id)
weighted_meth$gene_id <- gsub(".t.*","",weighted_meth$gene_id)
weighted_meth <- weighted_meth[(weighted_meth$feature == "exon_first3" | 
                                  weighted_meth$feature == "promotors_2000bp"),]
weighted_meth <- weighted_meth[,c(2,3,7,8)]
weighted_meth <- aggregate(weightedMeth.mean ~ feature+gene_id+origin, data=weighted_meth,
                            FUN = mean)
#write.table(weighted_meth, file="weightedMeth_exons_promotors_only.txt", sep="\t", quote = F,
#                      col.names = T, row.names = F)


setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/stevie_RSEM_counts")
file.list <- list.files("./", pattern = "*genes.results")
read_file1 <- function(x){
  read_delim(x, "\t", escape_double = FALSE, trim_ws = TRUE)
}
exp_samples <- lapply(file.list, read_file1)
setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/meth_paired_with_exp")

for (i in seq_along(exp_samples)){
  exp_samples[[i]] <- exp_samples[[i]][,c(1,7)]
}

colnames(exp_samples[[1]]) <- c("gene_id","F1_FPKM")
colnames(exp_samples[[2]]) <- c("gene_id","F2_FPKM")
colnames(exp_samples[[3]]) <- c("gene_id","F3_FPKM")
colnames(exp_samples[[4]]) <- c("gene_id","M1_FPKM")
colnames(exp_samples[[5]]) <- c("gene_id","M2_FPKM")
colnames(exp_samples[[6]]) <- c("gene_id","M3_FPKM")

exp_data <-Reduce(function(...) merge(..., by="gene_id", all=TRUE), exp_samples)
exp_data$female <- (exp_data$F1_FPKM + exp_data$F2_FPKM + exp_data$F3_FPKM)/3
exp_data$male <- (exp_data$M1_FPKM + exp_data$M2_FPKM + exp_data$M3_FPKM)/3
#exp_data$FPKM_diff <- exp_data$female_mean_FPKM - exp_data$male_mean_FPKM # + value = higher in female
exp_data <- exp_data[,c(1,8,9)]

exp_data_long <- melt(exp_data)
colnames(exp_data_long) <- c("gene_id","origin","FPKM")
#write.table(exp_data_long, file="FPKM_values_by_sex.txt", sep="\t", quote = F,
 #           col.names = T, row.names = F)

# ----------------------------------------------
# General meth vs general exp full scatter
# -----------------------------------------------
meth_exp <- merge(exp_data_long, weighted_meth, by=c("gene_id","origin"), all=T)
meth_exp <- meth_exp[!is.na(meth_exp$weightedMeth.mean & !is.na(meth_exp$FPKM)),]
#meth_exp$FPKM[meth_exp$FPKM==0] <- 1
meth_exp$logFPKM <- log(meth_exp$FPKM)
meth_exp$logFPKM[meth_exp$logFPKM == -Inf] <- 0

prom_data <- meth_exp[meth_exp$feature=="promotors_2000bp",]

ggplot(prom_data, aes(x=weightedMeth.mean, y=logFPKM, colour = origin))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  scale_colour_manual("",breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  ggtitle("Promotors")+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text = element_text(size=18),
        plot.title = element_text(size=18))

prom_data_males <- prom_data[prom_data$origin=="male",]
ggplot(prom_data_males, aes(x=weightedMeth.mean, y=logFPKM))+
  geom_point(colour="steelblue1",size=2)+
  geom_smooth(method = "lm",size=2)+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Promotors: Male")+
  theme_bw()+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26),
        plot.title = element_text(size=26))+
  xlim(0,1.0)

prom_data_females <- prom_data[prom_data$origin=="female",]
ggplot(prom_data_females, aes(x=weightedMeth.mean, y=logFPKM))+
  geom_point(colour="pink1",size=2)+
  geom_smooth(method = "lm", colour="deeppink",size=2)+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Promotors: Female")+
  theme_bw()+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26),
        plot.title = element_text(size=26))+
  xlim(0,1.0)




exon_data <- meth_exp[meth_exp$feature=="exon_first3",]
ggplot(exon_data, aes(x=weightedMeth.mean, y=logFPKM, colour = origin))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  scale_colour_manual("",breaks=c("female", "male"),
                      values = c("pink1","steelblue1"),
                      labels= c("Female", "Male"))+
  ggtitle("Exons 1-3")+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text = element_text(size=18),
        plot.title = element_text(size=18))

exon_data_males <- exon_data[exon_data$origin=="male",]
ggplot(exon_data_males, aes(x=weightedMeth.mean, y=logFPKM))+
  geom_point(colour="steelblue1",size=2)+
  geom_smooth(method = "lm",size=2)+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Exons 1-3: Male")+
  theme_bw()+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26),
        plot.title = element_text(size=26))+
  xlim(0,1.0)

exon_data_females <- exon_data[exon_data$origin=="female",]
ggplot(exon_data_females, aes(x=weightedMeth.mean, y=logFPKM))+
  geom_point(colour="pink1",size=2)+
  geom_smooth(method = "lm", colour="deeppink",size=2)+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Exons 1-3: Female")+
  theme_bw()+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26),
        plot.title = element_text(size=26))+
  xlim(0,1.0)


# Stats
model1<-lm(FPKM~origin*weightedMeth.mean, data=prom_data)
model2<-lm(FPKM~origin+weightedMeth.mean, data=prom_data)
anova(model1,model2) # No interaction
summary.lm(model2) # Sex not significiant but meth does predict exp


model1<-lm(FPKM~origin*weightedMeth.mean, data=exon_data)
model2<-lm(FPKM~origin+weightedMeth.mean, data=exon_data)
anova(model1,model2) # Sex not significiant but meth does predict exp
summary.lm(model2)


# ----------------------------------------------
# General meth vs general exp binned scatter
# -----------------------------------------------

# For exon methylation
meth_only <- subset(meth_exp, feature == "exon_first3")
female_data <- subset(meth_only, origin=="female")
male_data <- subset(meth_only, origin=="male")

female_data$bins <- as.numeric(cut2(female_data$weightedMeth.mean , g=100))
male_data$bins<-as.numeric(cut2(male_data$weightedMeth.mean , g=100))

female<-as.data.frame(aggregate(female_data$FPKM, by=list(female_data$bins), mean))
female$logfpkm <- log10(female$x)
colnames(female)<-c("meth_bin","FPKM","logFPKM")
female$status <- "female"
#plot(female$meth_bin~female$FPKM)

male<-as.data.frame(aggregate(male_data$FPKM, by=list(male_data$bins), mean))
male$logfpkm <- log10(male$x)
colnames(male)<-c("meth_bin","FPKM","logFPKM")
male$status <- "male"
#plot(male$meth_bin~male$FPKM)

final_data <- merge(male, female, by="meth_bin")
final_data <- rbind(male, female)


ggplot(data=final_data, aes(x=meth_bin, y=logFPKM, colour=status))+
  geom_point(size=2)+
  geom_smooth(method="loess", size=2)+
  xlab("Methylation Rank (Low to High)")+
  ylab("log(FPKM)")+
  ggtitle("Exons 1-3")+
  scale_colour_manual("", breaks=c("female","male"),
                      values = c("pink1","steelblue1"),
                      labels=c("Female","Male"))+
  theme_bw() +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26),
        legend.text=element_text(size=24),
        plot.title = element_text(size=26))


# Stats for exons
data1<-melt(final_data,id=c("meth_bin","status"))

model1<-lm(value~status*meth_bin, data=data1)
model2<-lm(value~status+meth_bin, data=data1)
anova(model1,model2)# Slight sig interaction
summary.lm(model1) # sex non-sig, sig FPKM and onlu just sig interaction with sex


# For promotor methylation
meth_only <- subset(meth_exp, feature == "promotors_2000bp")
female_data <- subset(meth_only, origin=="female")
male_data <- subset(meth_only, origin=="male")

female_data$bins <- as.numeric(cut2(female_data$weightedMeth.mean, g=100))
male_data$bins<-as.numeric(cut2(male_data$weightedMeth.mean, g=100))

female<-as.data.frame(aggregate(female_data$FPKM, by=list(female_data$bins), mean))
female$logfpkm <- log10(female$x)
colnames(female)<-c("meth_bin","FPKM","logFPKM")
female$status <- "female"
#plot(female$meth_bin~female$FPKM)

male<-as.data.frame(aggregate(male_data$FPKM, by=list(male_data$bins), mean))
male$logfpkm <- log10(male$x)
colnames(male)<-c("meth_bin","FPKM","logFPKM")
male$status <- "male"
#plot(male$meth_bin~male$FPKM)

final_data <- merge(male, female, by="meth_bin")
final_data <- rbind(male, female)


ggplot(data=final_data, aes(x=meth_bin, y=logFPKM, colour=status))+
  geom_point( size=2)+
  geom_smooth(method="loess", size=2)+
  xlab("Methylation Rank (Low to High)")+
  ylab("log(FPKM)")+
  ggtitle("Promotors")+
  scale_colour_manual("", breaks=c("female","male"),
                      values = c("pink1","steelblue1"),
                      labels=c("Female","Male"))+
  theme_bw() +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26),
        legend.text=element_text(size=24),
        plot.title = element_text(size=26))

# Line stats (maybe not applicable as did stats on scatters???)
data1<-melt(final_data,id=c("meth_bin","status"))

model1<-lm(value~status*meth_bin, data=data1)
model2<-lm(value~status+meth_bin, data=data1)
anova(model1,model2)# significant interaction
summary.lm(model1)#sex non-sig, sig FPKM and onlu just sig interaction with sex



# ----------------------------------------------
# Binned violin plots (as in Trine's paper)
# -----------------------------------------------
head(meth_exp)

meth_exp$bins <- "no_meth"
meth_exp$bins[meth_exp$weightedMeth.mean > 0.7] <- "high"
meth_exp$bins[meth_exp$weightedMeth.mean < 0.3 &
                meth_exp$weightedMeth.mean > 0] <- "low"
meth_exp$bins[meth_exp$weightedMeth.mean > 0.3 & 
                meth_exp$weightedMeth.mean < 0.7] <- "medium"

meth_exp$both_info <- as.factor(paste(meth_exp$origin, meth_exp$bins, sep = "_"))

condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

prom_data <- meth_exp[meth_exp$feature=="promotors_2000bp",]
exon_data <- meth_exp[meth_exp$feature=="exon_first3",]

ggplot(prom_data, aes(x=both_info, y=logFPKM))+
  geom_violin(aes(fill=origin))+
 # geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.02, dotsize = 0.2)+
  xlab("Methylation Bin")+
  ylab("log(FPKM)")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual("",breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_no_meth", "male_no_meth",
                            "female_low", "male_low",
                            "female_medium", "male_medium",
                            "female_high", "male_high"),
                   labels=c("None", "None",
                            "Low", "Low",
                            "Medium", "Medium","High", "High"))+
  ggtitle("Promotors")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=26),
        legend.text = element_text(size=24),
        plot.title = element_text(size=26))


ggplot(exon_data, aes(x=both_info, y=logFPKM))+
  geom_violin(aes(fill=origin))+
  #geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.025, dotsize = 0.3)+
  xlab("Methylation Bin")+
  ylab("log(FPKM)")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual("",breaks=c("female", "male"),
                    values = c("pink1","steelblue1"),
                    labels= c("Female", "Male"))+
  scale_x_discrete(limits=c("female_no_meth", "male_no_meth",
                            "female_low", "male_low",
                            "female_medium", "male_medium",
                            "female_high", "male_high"),
                   labels=c("None", "None",
                            "Low", "Low",
                            "Medium", "Medium","High", "High"))+
  ggtitle("Exons 1-3")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=26),
        legend.text = element_text(size=24),
        plot.title = element_text(size=26))

# Stats
meth_exp$origin <- as.factor(meth_exp$origin)
meth_exp$bins <- as.factor(meth_exp$bins)
prom_data <- meth_exp[meth_exp$feature=="promotors_2000bp",]
exon_data <- meth_exp[meth_exp$feature=="exon_first3",]

model1<-lm(FPKM ~ origin * bins, data=exon_data)
model2<-lm(FPKM ~ origin + bins, data=exon_data)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Low different to others, sex not sig
exon_data$SHD<-interaction(exon_data$origin,exon_data$bins)
model1_new<-lm(FPKM~-1+SHD, data=exon_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) # Low is diff to high and medium but not no-meth

model1<-lm(FPKM ~ origin * bins, data=prom_data)
model2<-lm(FPKM ~ origin + bins, data=prom_data)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Low and no-meth different to others, sex not sig
prom_data$SHD<-interaction(prom_data$origin,prom_data$bins)
model1_new<-lm(FPKM~-1+SHD, data=prom_data)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) # Low diff to high and medium but not no-meth as above
# Also female no-meth (but not male) diff to both medium and both high 


