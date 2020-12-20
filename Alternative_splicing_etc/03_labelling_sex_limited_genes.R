#---------------------------------------------------------
# Sorting out lists of genes which are sex-limited etc.
#---------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/transcription")

library(readr)
library(reshape2)
library(ggplot2)
#---------------------------------------------------------

FPKM_values_by_sex <- read_delim("FPKM_values_by_sex.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
FPKMs <- dcast(FPKM_values_by_sex, gene_id ~ origin)
head(FPKMs)

logFC_DEgenes <- read_delim("logFC_DEgenes.txt", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
head(logFC_DEgenes)
logFC_DEgenes <- logFC_DEgenes[,c(1,8,9)]
colnames(logFC_DEgenes)<-c("gene_id","log2FC","bias")

FPKMS_in_data <- merge(FPKMs, logFC_DEgenes, by= "gene_id")

#---------------------------------------------------------

male_limited_exp <- FPKMS_in_data[FPKMS_in_data$female ==0,] #204
female_limited_exp <- FPKMS_in_data[FPKMS_in_data$male ==0,] #150

FPKMS_in_data$bias[FPKMS_in_data$gene_id %in% male_limited_exp$gene_id] <- "male_limited_exp"
FPKMS_in_data$bias[FPKMS_in_data$gene_id %in% female_limited_exp$gene_id] <- "female_limited_exp"
unique(FPKMS_in_data$bias)

FPKMS_in_data$bias[FPKMS_in_data$log2FC < -10 &
                     !(FPKMS_in_data$gene_id %in% male_limited_exp$gene_id)] <- "extreme_male_biased"
FPKMS_in_data$bias[FPKMS_in_data$log2FC > 10 &
                     !(FPKMS_in_data$gene_id %in% female_limited_exp$gene_id)] <- "extreme_female_biased"
#write.table(FPKMS_in_data, file="FPKMs_logFC_bias_catergory.txt", sep="\t", quote = F, col.names = T,
 #           row.names = F)
table(FPKMS_in_data$bias)

# Goodness of fit
observed = c(150, 204)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) 



#---------------------------------------------------------
FPKMS_in_data$main_category <- "Unbias"
FPKMS_in_data$main_category[FPKMS_in_data$bias == "Mbias" |
                              FPKMS_in_data$bias == "extreme_male_biased" |
                              FPKMS_in_data$bias == "male_limited_exp"] <- "Mbias"
FPKMS_in_data$main_category[FPKMS_in_data$bias == "Fbias" |
                              FPKMS_in_data$bias == "extreme_female_biased" |
                              FPKMS_in_data$bias == "female_limited_exp"] <- "Fbias"
FPKMS_in_data$plotbias <- "main_colour"
FPKMS_in_data$plotbias[FPKMS_in_data$bias == "extreme_female_biased" |
                     FPKMS_in_data$bias=="extreme_male_biased"] <- "extreme_colour"
FPKMS_in_data$plotbias[FPKMS_in_data$bias == "female_limited_exp" |
                         FPKMS_in_data$bias=="male_limited_exp"] <- "limited_colour"
FPKMS_in_data_plot <- FPKMS_in_data[!FPKMS_in_data$bias=="Unbias",]

ggplot(FPKMS_in_data_plot, aes(x=main_category, fill=plotbias))+
  geom_bar()+
  theme_bw()+
  xlab("Gene Expression Category")+
  ylab("Number of Genes")+
  scale_fill_manual("",breaks=c("main_colour", "extreme_colour","limited_colour"),
                      labels = c("Fold-change > 1.5","Fold-change > 10","Sex-Limited"),
                      values = c("grey75", "grey45","grey14"))+
  theme_bw()+
  scale_x_discrete(limits=c("Fbias", "Mbias"),
                   labels=c("Female Biased", "Male Biased"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))

#---------------------------------------------------------
# Make something like Andrew's graph
FPKM_values_by_sex <- read_delim("FPKM_values_by_sex.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
FPKMs <- dcast(FPKM_values_by_sex, gene_id ~ origin)
head(FPKMs)

FPKMs <- FPKMs[FPKMs$male>1 | FPKMs$female>1 ,] #19273
FPKMs$total <- FPKMs$female + FPKMs$male
FPKMs$prop_female <- 1- (FPKMs$total - FPKMs$female)/FPKMs$total

ggplot(FPKMs, aes ( x= prop_female))+
  geom_histogram(colour="black", bins=50)+
  xlab("Proportion of Female Expression")+
  ylab("Number of Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))

# Try Andrew's method
FPKMs$female_sq <- FPKMs$female*FPKMs$female
FPKMs$male_sq <- FPKMs$male*FPKMs$male
FPKMs$SPM <- FPKMs$female_sq/(FPKMs$female_sq +FPKMs$male_sq )

ggplot(FPKMs, aes ( x= SPM))+
  geom_histogram(colour="black", bins=50)+
  xlab("SPM in Females")+
  ylab("Number of Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))

#---------------------------------------------------------
# Hmmm
FPKMS_in_data$both <- FPKMS_in_data$male + FPKMS_in_data$female
FPKMS_in_data$prop <-  FPKMS_in_data$female / FPKMS_in_data$both
hist(FPKMS_in_data$prop)

FPKMs$both <- FPKMs$male + FPKMs$female
FPKMs <- FPKMs[FPKMs$male>1 |FPKMs$female>1 ,]
FPKMs$prop <-  FPKMs$female / FPKMs$both
hist(FPKMs$prop)

look <- FPKMs[!FPKMs$gene_id %in% logFC_DEgenes$gene_id,]
hist(look$prop) # What are these genes which aren't in the final diff exp output ... they range in values ...
