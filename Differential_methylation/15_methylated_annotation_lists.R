## -------------------------------------------------------------------------
# Take average weighted methylation level of feature across bio replicates
## -------------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation")
library(readr)
library(doBy)
library(dplyr)

# Determine threshold of saying something is methylated, then apply that threshold to get methyated lists 
# of different annotations such as exons

annotation <- read_delim("PCITRI_weighted_meth_annotation_by_sex.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
annotation$feature <- as.factor(annotation$feature)
annotation_males <- subset(annotation, origin =="male")
annotation_females <- subset(annotation, origin =="female")

hist(annotation_males$weightedMeth.mean[annotation_males$feature == "exon_first3"], breaks=30)
hist(annotation_females$weightedMeth.mean[annotation_females$feature == "exon_first3"], breaks=30)

test <- annotation_females[annotation_females$weightedMeth.mean > 0.1,] #10,000-6,000 between nearest bars
hist(test$weightedMeth.mean[test$feature == "exon_first3"],breaks=30)

test2 <- annotation_females[annotation_females$weightedMeth.mean > 0.00,] #30,000-10,000
hist(test2$weightedMeth.mean[test2$feature == "exon_first3"],breaks=30) # Looking back at original histogram this looks like a good cut off
# Given how weighted meth work, this means 5% of all reads in a CpG context for this whole regions
# show a C and are therefore representarive of methylation 

### AHHHH what if we use the genome-wide level to determine the cuf-off (think about the logic...)
# This is 8.28% in females and 9.7 in males.
# We still get a little bar below 0.1 when we use this
par(mfrow=c(1,2))
test3 <- annotation_females[annotation_females$weightedMeth.mean > 0.1,] 
hist(test3$weightedMeth.mean[test3$feature == "exon_first3"],breaks=30, main="female exons 1-3",
     xlab = "Mehtylation Level")
test4 <- annotation_males[annotation_males$weightedMeth.mean > 0.1,] 
hist(test4$weightedMeth.mean[test4$feature == "exon_first3"],breaks=30, main="male exons 1-3",
     xlab = "Mehtylation Level")

# I'm thinking this works!!! We get a peack at the low end but then a drop off before we get the big peak again
# at 0, which makes me think most sites have this level or greater if they're methylated, anything lower
# could just be error before reaching the peak of 0, otherwise is hould also be high, no?

# Have a read at other literature that uses weighted meth to see what they use as a cut-off. 
# Can't find anything, mostly arbituary values

# Check this cut-off is also appropriate for promotors (hmmm maybe should be 0.05)
test3 <- annotation_females[annotation_females$weightedMeth.mean > 0.1,] 
hist(test3$weightedMeth.mean[test3$feature == "promotors_2000bp"],breaks=30, main="female promotors",
     xlab = "Mehtylation Level")
test4 <- annotation_males[annotation_males$weightedMeth.mean > 0.1,] 
hist(test4$weightedMeth.mean[test4$feature == "promotors_2000bp"],breaks=30, main="male promotors",
     xlab = "Mehtylation Level")


test3 <- annotation_females[annotation_females$weightedMeth.mean > 0.1,] 
hist(test3$weightedMeth.mean[test3$feature == "TE"],breaks=30, main="female TEs",
     xlab = "Mehtylation Level")
test4 <- annotation_males[annotation_males$weightedMeth.mean > 0.1,] 
hist(test4$weightedMeth.mean[test4$feature == "TE"],breaks=30, main="male TEs",
     xlab = "Mehtylation Level")


#  Keep only annotations found both in male and female
both <- merge(annotation_males, annotation_females, by=c("gene_id","start","end"))
both <- both[,-c(7,9,10,11,12)]
colnames(both) <- c("gene_id","start","end","scaffold","feature","cpg_count","male_mean_weightedMeth","female_mean_weightedMeth")

both$methylated <- "no"
both$methylated[both$male_mean_weightedMeth > 0.1 | both$female_mean_weightedMeth > 0.1] <- "yes"

write.table(both, file = "PCITRI_weighted_meth_annotation_common_both_sexes.txt", sep="\t",
            col.names = T, row.names = F, quote = F)


## -------------------------------------------------------------------------
# Make some boxplots of methylation level of diff genomic features
## -------------------------------------------------------------------------

# This is for features that are methylated in either sex (>0.1 in females or males)
library(ggplot2)
library(reshape2)

methylated_stuff <- subset(both, methylated == "yes")
methylated_stuff <- methylated_stuff[,c(5,7,8)]
melted_meth_stuff <- melt(methylated_stuff)
colnames(melted_meth_stuff) <- c("Feature","Sex","Weighted_Meth")
melted_meth_stuff$Sex <- gsub("_mean_weightedMeth","",melted_meth_stuff$Sex)
melted_meth_stuff <- melted_meth_stuff[!melted_meth_stuff$Feature=="promotors_1000bp",]
melted_meth_stuff <- melted_meth_stuff[!melted_meth_stuff$Feature=="whole_gene",]


ggplot(melted_meth_stuff, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_boxplot()+
  theme_bw()+
  ylab("Mean Weighted Methylation Level per Feature")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size = 16),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE") )


# Throw in some stats
res.1 <- aov(Weighted_Meth ~ Feature * Sex, data = melted_meth_stuff)
summary(res.1) # Interaction is significant

TukeyHSD(res.1)
# Overall males < female, p < 0.000 

# Exons 1-3: male < female, p < 0.000
# Promtors: male < female, p < 0.000
# Exons 3+: male < female, p < 0.000
# Introns: not significant, p = 0.995
# TEs: male < female, p < 0.000



# Boxplot for meth of all features, methylated in one sex or not 
methylated_stuff <- both[,c(5,7,8)]
melted_meth_stuff <- melt(methylated_stuff)
colnames(melted_meth_stuff) <- c("Feature","Sex","Weighted_Meth")
melted_meth_stuff$Sex <- gsub("_mean_weightedMeth","",melted_meth_stuff$Sex)
melted_meth_stuff <- melted_meth_stuff[!melted_meth_stuff$Feature=="promotors_1000bp",]
melted_meth_stuff <- melted_meth_stuff[!melted_meth_stuff$Feature=="whole_gene",]

library(ggforce)

ggplot(melted_meth_stuff, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_boxplot()+
  theme_bw()+
  ylab("Mean Weighted Methylation Level per Feature")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE"))+
  facet_zoom(ylim = c(0, 0.1))



# Throw in some stats
res.1 <- aov(Weighted_Meth ~ Feature * Sex, data = melted_meth_stuff)
summary(res.1) # Interaction is significant

TukeyHSD(res.1)
# Overall males > female, p < 0.000

# Exons 1-3: male > female, p < 0.000
# Promtors: not significant p = 0.991
# Exons 3+: male > female, p < 0.000
# Introns: male > female, p < 0.000
# TEs: male > female, p < 0.000

## -------------------------------------------------------------------------
# Make the same plots but as bar plots
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

summary_all<-summarySE(melted_meth_stuff, measurevar = "Weighted_Meth", 
                       groupvars = c("Feature","Sex"))

ggplot(summary_all, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Meth-ci, ymax=Weighted_Meth+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Mean Weighted Methylation Level per Feature")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE"))


methylated_stuff <- subset(both, methylated == "yes")
methylated_stuff <- methylated_stuff[,c(5,7,8)]
melted_meth_stuff <- melt(methylated_stuff)
colnames(melted_meth_stuff) <- c("Feature","Sex","Weighted_Meth")
melted_meth_stuff$Sex <- gsub("_mean_weightedMeth","",melted_meth_stuff$Sex)
melted_meth_stuff <- melted_meth_stuff[!melted_meth_stuff$Feature=="promotors_1000bp",]
melted_meth_stuff <- melted_meth_stuff[!melted_meth_stuff$Feature=="whole_gene",]

summary_all<-summarySE(melted_meth_stuff, measurevar = "Weighted_Meth", 
                       groupvars = c("Feature","Sex"))

ggplot(summary_all, aes(x=Feature, y=Weighted_Meth, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Meth-ci, ymax=Weighted_Meth+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Mean Weighted Methylation Level per Feature")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE"),
                   labels = c("Promtors","Exons 1-3","Exons 4+","Introns","TEs"),
                   limits =c("promotors_2000bp","exon_first3","exon_notFirst3","intron","TE"))
