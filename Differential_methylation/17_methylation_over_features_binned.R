## -------------------------------------------------------------------------
# Re-make meth over feature graph for diff levels of methylation
## -------------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation")
library(readr)
library(doBy)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Hmisc)
library(scales)
## -------------------------------------------------------------------------
annotation <- read_delim("PCITRI_weighted_meth_annotation_by_sex.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
annotation$feature <- as.factor(annotation$feature)
annotation_males <- subset(annotation, origin =="male")
annotation_females <- subset(annotation, origin =="female")

both <- merge(annotation_males, annotation_females, by=c("gene_id","start","end"))
both <- both[,-c(7,9,10,11,12)]
colnames(both) <- c("gene_id","start","end","scaffold","feature","cpg_count","male_mean_weightedMeth","female_mean_weightedMeth")

both$methylated <- "no"
both$methylated[both$male_mean_weightedMeth > 0.1 | both$female_mean_weightedMeth > 0.1] <- "yes"
methylated_stuff <- both[,c(5,7,8)]
melted_meth_stuff <- melt(methylated_stuff)
colnames(melted_meth_stuff) <- c("Feature","Sex","Weighted_Meth")
melted_meth_stuff$Sex <- gsub("_mean_weightedMeth","",melted_meth_stuff$Sex)
melted_meth_stuff <- melted_meth_stuff[!melted_meth_stuff$Feature=="whole_gene",]

## -------------------------------------------------------------------------
#### Define summary function (ref:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/)
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
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 
## -------------------------------------------------------------------------
# Normal graph with all information

summary_all<-summarySE(melted_meth_stuff, measurevar = "Weighted_Meth", 
                       groupvars = c("Feature","Sex"))

ggplot(summary_all, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Meth-ci, ymax=Weighted_Meth+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Weighted Methylation Level")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

# Stats
library(multcomp)
model1<-lm(Weighted_Meth ~ Sex * Feature , data=melted_meth_stuff)
model2<-lm(Weighted_Meth ~ Sex + Feature , data=melted_meth_stuff)
anova(model1,model2) # Sig interaction
summary.lm(model1) # Everything sig

for_stats <- melted_meth_stuff
for_stats$SHD<-interaction(for_stats$Sex,
                           for_stats$Feature)
model1_new<-lm(Weighted_Meth~-1+SHD, data=for_stats)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) 
# Everything is different except: 
# male intergenic and female exon 1-3
# male TE and male exon 1-3
# male intron and female intergenic
# male and female promotors

## -------------------------------------------------------------------------
# Different categories
head(melted_meth_stuff)
meth_low <- melted_meth_stuff[melted_meth_stuff$Weighted_Meth <= 0.3,]
meth_medium <- melted_meth_stuff[melted_meth_stuff$Weighted_Meth > 0.3 &
                                   melted_meth_stuff$Weighted_Meth < 0.7,]
meth_high <- melted_meth_stuff[melted_meth_stuff$Weighted_Meth >= 0.7,]

summary_all_high<-summarySE(meth_high, measurevar = "Weighted_Meth", 
                       groupvars = c("Feature","Sex"))

a1 <- ggplot(summary_all_high, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Meth-ci, ymax=Weighted_Meth+ci),
                width=.2,
                position = position_dodge(.9))+
  #geom_text(aes(label = summary_all$N, x = Feature), 
   #         position = position_dodge(width = 0.9), vjust = -1)+ 
  theme_bw()+
 # xlab("")+
 # ylab("Mean Weighted Methylation Level per Feature")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

summary_all_med<-summarySE(meth_medium, measurevar = "Weighted_Meth", 
                            groupvars = c("Feature","Sex"))
a2 <- ggplot(summary_all_med, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Meth-ci, ymax=Weighted_Meth+ci),
                width=.2,
                position = position_dodge(.9))+
  #geom_text(aes(label = summary_all$N, x = Feature), 
  #         position = position_dodge(width = 0.9), vjust = -1)+ 
  theme_bw()+
  #xlab("Genomic Feature")+
  #ylab("Mean Weighted Methylation Level per Feature")+
  ggtitle("Medium Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

summary_all_low<-summarySE(meth_low, measurevar = "Weighted_Meth", 
                           groupvars = c("Feature","Sex"))
a3 <- ggplot(summary_all_low, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Meth-ci, ymax=Weighted_Meth+ci),
                width=.2,
                position = position_dodge(.9))+
  #geom_text(aes(label = summary_all$N, x = Feature), 
  #         position = position_dodge(width = 0.9), vjust = -1)+ 
  theme_bw()+
  #xlab("Genomic Feature")+
  #ylab("Mean Weighted Methylation Level per Feature")+
  ggtitle("Low Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

## -------------------------------------------------------------------------
# Frequency of features being high/low methylated
head(melted_meth_stuff)
meth_low <- melted_meth_stuff[melted_meth_stuff$Weighted_Meth < 0.3,]
meth_medium <- melted_meth_stuff[melted_meth_stuff$Weighted_Meth > 0.3 &
                                   melted_meth_stuff$Weighted_Meth < 0.7,]
meth_high <- melted_meth_stuff[melted_meth_stuff$Weighted_Meth > 0.7,]
meth_none <- melted_meth_stuff[melted_meth_stuff$Weighted_Meth ==0,]

melted_meth_stuff$bins<-"low"
melted_meth_stuff$bins[melted_meth_stuff$Weighted_Meth > 0.3 &
                    melted_meth_stuff$Weighted_Meth < 0.7] <-"medium"
melted_meth_stuff$bins[melted_meth_stuff$Weighted_Meth > 0.7] <-"high"
melted_meth_stuff$bins[melted_meth_stuff$Weighted_Meth ==0] <-"none"

melted_meth_stuff$combined <- paste0(melted_meth_stuff$Feature, "_", melted_meth_stuff$Sex)
melted_meth_stuff_2 <- melted_meth_stuff[,-c(3,5)]
melted_meth_stuff_2$counts <- with(melted_meth_stuff_2, 
                                  ave(bins, Feature, Sex, bins, FUN=length))
plot_data <- melted_meth_stuff_2[!duplicated(melted_meth_stuff_2),]

plot_data_prom <- subset(plot_data, bins =="low")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b1<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("Genomic Feature")+
#  ylab("Count")+
  ggtitle("Low Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

plot_data_prom <- subset(plot_data, bins =="medium")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b2<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  #xlab("Genomic Feature")+
  #ylab("Count")+
  ggtitle("Medium Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

plot_data_prom <- subset(plot_data, bins =="high")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b3<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("")+
  ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title.y=element_text(size=12),
        axis.title.x = element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Count")+
 # ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))


## -------------------------------------------------------------------------
# Make one nice figure (or two?)
library(ggpubr)
levels_across_feature <- ggarrange(a1, a2, a3,
                        ncol=3, nrow=1, common.legend = TRUE, legend="right")

annotate_figure(levels_across_feature, 
                left = text_grob("Mean Weighted Methylation Level per Feature", 
                                 color = "black", rot = 90, size=14),
                bottom = text_grob("Genomic Feature", 
                                   color = "black", size =12,
                                   hjust = 0.75 ))

## -------------------------------------------------------------------------
plot_data_prom <- subset(plot_data, bins =="none")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b4 <- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("Genomic Feature")+
  #  ylab("Count")+
  ggtitle("No Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))

plot_data_prom <- subset(plot_data, bins =="high")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b3_1<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("")+
 # ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs","Intergenic"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE","intergenic"))


## -------------------------------------------------------------------------
# Slightly diff figures
levels_across_feature <- ggarrange(b3_1,b2,b1,b4,
                                   ncol=2, nrow=2, common.legend = TRUE, legend="right")

annotate_figure(levels_across_feature, 
                 left = text_grob("Count", 
                                 color = "black", rot = 90, size=14),
                bottom = text_grob("Genomic Feature", 
                                   color = "black", size =12,
                                   hjust = 0.75 ))
