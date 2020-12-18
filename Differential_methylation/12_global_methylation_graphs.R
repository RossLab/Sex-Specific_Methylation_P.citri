## -------------------------------------------------------------------------
## Making Fancy Genome-Methyation Differences Graphs
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Read in the subsetted file from methylkit 
objectmethbase1 <- read_delim("./differential_meth_methylkit/F_vs_M_objectmethbase.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase1$chrBase <- paste(objectmethbase1$chr, ".", objectmethbase1$start, sep="")
objectmethbase1 <- objectmethbase1[,-3]
objectmethbase1$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase1$strand)


## -------------------------------------------------------------------------

# Setset out each sample (c4 and c6)
a <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs1","numTs1")]
b <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs2","numTs2")]
c <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs3","numTs3")]
d <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs4","numTs4")]
e <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs5","numTs5")]
f <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs6","numTs6")]
g <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs7","numTs7")]
h <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs8","numTs8")]
j <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs9","numTs9")]

## -------------------------------------------------------------------------

# Write out each sample
all_files <- list(a,b,c,d,e,f,g,h,j)

for(i in seq_along(all_files)){
  colnames(all_files[[i]])[c(3,5,6,7)] <- c("base","coverage","numCs","numTs")
  all_files[[i]]$freqC <- round((all_files[[i]]$numCs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]]$freqT <- round((all_files[[i]]$numTs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]] <- all_files[[i]][-c(6,7)]
  myfile <- file.path("./", paste0(i,"_","subsetted_final.txt"))
  write.table(all_files[[i]], file=myfile, quote=F, sep="\t", row.names=F)
}

# Urgh need to go and rename each file with sample name and condition,
# check the methylkit scritps for the order of samples

## -------------------------------------------------------------------------
# Use methylkit to get the data all togther and plottable 
setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation/differential_meth_methylkit")
file.listA <- list("F1_subsetted_final.txt","F2_subsetted_final.txt","F3_subsetted_final.txt",
                   "F4_subsetted_final.txt","F5_subsetted_final.txt","M1_subsetted_final.txt",
                   "M3_subsetted_final.txt","M4_subsetted_final.txt","M5_subsetted_final.txt")

sample_list <- list("F1", "F2", "F3", "F4", "F5", "M1", "M3", "M4", "M5")

# Make a list of all genotypes in the right order = genotype_list

raw_data <- methRead(file.listA,
                     sample.id = sample_list,
                     treatment = c(rep(0, 5), rep(2, 4)),
                     assembly="PCITRI_v0", 
                     context="CpG")

meth_all_data <- unite(raw_data)#2774656 CpGs mmmm


## -------------------------------------------------------------------------
# Make a nice dendogram

clusterSamples(meth_all_data, dist="correlation", method="ward", plot=TRUE)

hc <- clusterSamples(meth_all_data, dist="correlation", method="ward", plot=FALSE)

data1 <- dendro_data(hc)
labs <- label(data1)

labs$Sex<- c(rep("Female", 5), rep("Male", 4))

ggplot(segment(data1)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=1,
               show.legend = F)+
  geom_text(data=labs,
            aes(label=label, x=x, y=-0.009, colour=Sex,fontface="bold"),
            show.legend = T, angle = 90, size = 7, hjust = 1)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text = element_blank(),
        legend.position = "none")+
  scale_colour_manual(values=c("pink1", "steelblue1"))

## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(meth_all_data, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)

PCA_data1$Sex <- c(rep("Female", 5), rep("Male", 4))


percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC3, PC4, colour=Sex))+
  geom_point(size=14)+
  geom_text_repel(aes(label=sample), size=14,show.legend=FALSE, 
                  point.padding = 0.85, box.padding = 0.25)+
  theme_bw()+
  xlab(paste0("PC1: ",percentage[3],"variance")) +
  ylab(paste0("PC2: ",percentage[4],"variance")) +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.text=element_text(size=30),
        legend.title=element_blank())+
  scale_colour_manual(values=c("pink1", "steelblue1"))

# Include a coefficent of variation to measure the dispersion across individual CpGs in males and females
# Function from: https://rcompanion.org/rcompanion/c_02.html
summary.list = function(x)list(
  N.with.NA.removed= length(x[!is.na(x)]),
  Count.of.NA= length(x[is.na(x)]),
  Mean=mean(x, na.rm=TRUE),
  Median=median(x, na.rm=TRUE),
  Max.Min=range(x, na.rm=TRUE),
 # Range=max(Data$ Fish, na.rm=TRUE) - min(Data$ Fish, na.rm=TRUE),
  Variance=var(x, na.rm=TRUE),
  Std.Dev=sd(x, na.rm=TRUE),
  Coeff.Variation.Prcnt=sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)*100,
  Std.Error=sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)])),
  Quantile=quantile(x, na.rm=TRUE)
)

meth_levels_by_sex <- read_delim("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/methylation/regional_genome_browser/meth_levels_by_sex.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
head(meth_levels_by_sex)
summary.list(meth_levels_by_sex$Male) #138.21
summary.list(meth_levels_by_sex$Female)#217.4
# Female more vairable but I guess this makes sense due to the bimodal nature of female data

# Should probably do it by sample
F_vs_M_objectmethbase <- read_delim("F_vs_M_objectmethbase.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
head(F_vs_M_objectmethbase)
F_vs_M_objectmethbase$female_1 <- 1 - ((F_vs_M_objectmethbase$coverage1 - F_vs_M_objectmethbase$numCs1) / F_vs_M_objectmethbase$coverage1)
F_vs_M_objectmethbase$female_2 <- 1 - ((F_vs_M_objectmethbase$coverage2 - F_vs_M_objectmethbase$numCs2) / F_vs_M_objectmethbase$coverage2)
F_vs_M_objectmethbase$female_3 <- 1 - ((F_vs_M_objectmethbase$coverage3 - F_vs_M_objectmethbase$numCs3) / F_vs_M_objectmethbase$coverage3)
F_vs_M_objectmethbase$female_4 <- 1 - ((F_vs_M_objectmethbase$coverage4 - F_vs_M_objectmethbase$numCs4) / F_vs_M_objectmethbase$coverage4)
F_vs_M_objectmethbase$female_5 <- 1 - ((F_vs_M_objectmethbase$coverage5 - F_vs_M_objectmethbase$numCs5) / F_vs_M_objectmethbase$coverage5)

F_vs_M_objectmethbase$male_1 <- 1 - ((F_vs_M_objectmethbase$coverage6 - F_vs_M_objectmethbase$numCs6) / F_vs_M_objectmethbase$coverage6)
F_vs_M_objectmethbase$male_3 <- 1 - ((F_vs_M_objectmethbase$coverage7 - F_vs_M_objectmethbase$numCs7) / F_vs_M_objectmethbase$coverage7)
F_vs_M_objectmethbase$male_4 <- 1 - ((F_vs_M_objectmethbase$coverage8 - F_vs_M_objectmethbase$numCs8) / F_vs_M_objectmethbase$coverage8)
F_vs_M_objectmethbase$male_5 <- 1 - ((F_vs_M_objectmethbase$coverage9 - F_vs_M_objectmethbase$numCs9) / F_vs_M_objectmethbase$coverage9)

summary.list(F_vs_M_objectmethbase$female_1) #218.16
summary.list(F_vs_M_objectmethbase$female_2) #220.66
summary.list(F_vs_M_objectmethbase$female_3) #214.63
summary.list(F_vs_M_objectmethbase$female_4) #224.69
summary.list(F_vs_M_objectmethbase$female_5) #240.93

summary.list(F_vs_M_objectmethbase$male_1) #150.38
summary.list(F_vs_M_objectmethbase$male_3) #143.98
summary.list(F_vs_M_objectmethbase$male_4) #154.03
summary.list(F_vs_M_objectmethbase$male_5) #157.07
## -------------------------------------------------------------------------
# Make a plot to show the methylation level overall for males and females
# Levels taken from the alignment output bismark files
library(reshape2)

Female <- c(8.2, 7.9,8,7.6,7.3)
Male <- c(9.5,9.3,8.9,9.4,NA)

levels <-data.frame(Female, Male)
levels2 <- melt(levels)
levels2 <- levels2[!is.na(levels2$value),]
colnames(levels2) <- c("Sex","Methylation")

ggplot(levels2, aes(y=Methylation, x=Sex,fill=Sex))+
  geom_boxplot()+
  theme_bw()+
  xlab("Sex") +
  ylab("CpG Methylation Level") +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.text=element_text(size=30),
        legend.position = "none")+
  scale_fill_manual(values=c("pink1", "steelblue1"))
head(levels2)

t.test(Methylation~Sex,data=levels2)
# t = -7.1724, df = 6.9889, p-value = 0.0001831
