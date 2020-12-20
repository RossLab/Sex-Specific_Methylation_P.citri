#---------------------------------------------------------
# Alternative Splicing Analysis: Sex-Specific Mealybug
#---------------------------------------------------------

# NOTE: used new isoform counts from RSEM and new PASA annotation

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/transcription/new_alt_splice_analysis/input_files")
library(IsoformSwitchAnalyzeR)
library(readr)
library(doBy)
library(dplyr)
library(tidyr)

file.list <- list.files("./", pattern = "*isoforms.results")

quant_data <- importIsoformExpression(sampleVector = file.list)

myDesign <- data.frame(
  sampleID = colnames(quant_data$abundance)[-1],
  condition = c("female","female","female","male","male","male"))

# Deal with the gff file
# grep "exon" pcitri_v0_updated_with_pasa_no_hashes.gff3 > exons
# sed 's/ID/gene_id/g' exons > exons1
# sed 's/=/ \"/g' exons1 > exons2
# sed 's/Parent/ transcript_id /g' exons2 > exons3
# sed 's/\.t.*;/\";/1' exons3 > exons4
# sed 's/$/\"/' exons4 > pcitri_v0_updated_with_pasa_dexseq.gtf 

aSwitchList <- importRdata(
  isoformCountMatrix   = quant_data$counts,
  isoformRepExpression = quant_data$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "pcitri_v0_updated_with_pasa_dexseq.gtf",
  showProgress = TRUE)

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 3,
  removeSingleIsoformGenes = TRUE,
  dIFcutoff = 0.25) 
# Round 1 with initial annotation pcitri_v0
# 1235 isoforms left (as vast majortity of genes not annotaed with alternate transcripts for P.citri)
# 28586 transcripts (93.13%) are single isoform genes anyway
# 847 transcripts filtered on expression levels 

# Round 2 with new annotation pcitri_v0 updated with pasa
# The filtering removed 31514 ( 94.29% of ) transcripts which are single isoform
# 1210 isoforms left LOL
# 1920 left after just removing single isoforms
# 1210 left on expression filtering

# Analyse remaining isoforms with DEXSeq for differntial usage
SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE)

extractSwitchSummary(SwitchListAnalyzed)
#Comparison nrIsoforms nrGenes : round 1
#1 female vs male        423      209
# 423 isoforms from 209 genes differentially used
# out of 1235 isforms testes

#       Comparison nrIsoforms nrSwitches nrGenes: round 2
# 1 female vs male        408        227     203

switchingIso <- extractTopSwitches( 
  SwitchListAnalyzed, 
  filterForConsequences = F, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = T,    # when FALSE isoforms are returned
  sortByQvals = TRUE)

write.table(file="list_diff_iso_genes_round2.txt",switchingIso$gene_id,
            sep = "\t", quote = F, col.names = T, row.names = F)

# Compare to round 1:
list_diff_iso_genes_round2 <- read_csv("list_diff_iso_genes_round2.txt")
list_diff_iso_genes <- read_csv("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/transcription/stevie_RSEM_counts/list_diff_iso_genes.txt")

look <- merge (list_diff_iso_genes_round2, list_diff_iso_genes, by="x")
# 192 in common out of # 209 round 1 and 203 round 2

# Find out the top two isoforms by significant qvalue
top_switches <- extractTopSwitches(
  SwitchListAnalyzed, 
  filterForConsequences = F, 
  n = 2, 
  sortByQvals = TRUE)

# Make a plot of one of the isoforms
switchPlot(SwitchListAnalyzed, gene = 'g3376')

# Make a plot of the top 10 isoforms (automatically outputs as pdf)
switchPlotTopSwitches(
  switchAnalyzeRlist = SwitchListAnalyzed, 
  n = 203,
  filterForConsequences = FALSE, 
  splitFunctionalConsequences = F)


write.table(file="FPKM_of_all_sig_isoforms.txt",SwitchListAnalyzed$isoformRepExpression,
            sep = "\t", quote = F, col.names = T, row.names = F)

# How many of these are also just generally diff exp genes?
head(switchingIso)
sig_isoforms <- as.data.frame(switchingIso$gene_id)
colnames(sig_isoforms)<-"gene_id"
FPKMs_logFC_bias_catergory <- read_delim("../FPKMs_logFC_bias_catergory.txt", 
                                         "\t", escape_double = FALSE, col_names = TRUE, 
                                         trim_ws = TRUE)
diff_exp_data <- FPKMs_logFC_bias_catergory[,c(1,5)]
unique(diff_exp_data$bias)

both <- merge(sig_isoforms, diff_exp_data, by="gene_id") #missing 19 genes ...
table(both$bias)
#write.table(both, file="alt_splice_with_exp_bias.txt",quote = F, col.names = T, row.names = F, sep = "\t")

# Majority have male biased exp!!!!! VERY INTERESTING
# The missing 19 genes are just not sig diff exp, see Stevie's email, so add to below
df <- data.frame(
  group = c("Male Biased", "Unbiased","Female Biased"),
  value = c(104, 97, 8)
)
head(df)

library(scales)
library(ggplot2)

ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Set3")+
  scale_fill_manual("",breaks=c("Male Biased", "Unbiased","Female Biased"),
                    values = c("steelblue1","#BEBADA","pink1"))+
  geom_text(aes( label = value), color = "black", size=20,position = position_stack(vjust = 0.5))+
  theme_void()+
  ggtitle("Alternatively Spliced Genes")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=36),
        plot.title = element_text(size=30,vjust = -5))

# For male bias, pot = total genes (19282), total male bias, total alt splice(209) and overlap
# Total bias = 10548 
phyper(112,10548,19078,209, lower.tail = F) #3.510096e-08 sig number of alt splice are also diff exp

observed = c(104, 8)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) #X-squared = 82.286, df = 1, p-value < 2.2e-16

# Have a look to see if males use more isoforms than females
iso_exp_data <- SwitchListAnalyzed$isoformRepIF
tail(iso_exp_data)
iso_exp_data$gene <- iso_exp_data$isoform_id
iso_exp_data$gene <- substr(iso_exp_data$gene,1,nchar(iso_exp_data$gene)-3)
iso_exp_data$female <- (iso_exp_data$PC_F1 + iso_exp_data$PC_F2 + iso_exp_data$PC_F3)/3
iso_exp_data$male <- (iso_exp_data$PC_M1 + iso_exp_data$PC_M2 + iso_exp_data$PC_M3)/3
iso_exp_data <- iso_exp_data[,c(1,8,9,10)]

# asking of all isoforms, how many not exp in either male/female
nrow(iso_exp_data[iso_exp_data$female==0,]) #12
nrow(iso_exp_data[iso_exp_data$male==0,]) #7
observed = c(12, 7)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions
chisq.test(x = observed,
           p = expected) # X-squared = 1.3158, df = 1, p-value = 0.2513


# IGNORE BELOW
# Have a look to see if males use more isoforms than females
iso_exp_data <- SwitchListAnalyzed$isoformRepIF
tail(iso_exp_data)

iso_exp_data$gene <- iso_exp_data$isoform_id
iso_exp_data$gene <- substr(iso_exp_data$gene,1,nchar(iso_exp_data$gene)-3)

iso_exp_data$female <- (iso_exp_data$PC_F1 + iso_exp_data$PC_F2 + iso_exp_data$PC_F3)/3
iso_exp_data$male <- (iso_exp_data$PC_M1 + iso_exp_data$PC_M2 + iso_exp_data$PC_M3)/3

iso_exp_data <- iso_exp_data[,c(1,8,9,10)]

iso_exp_data$prop_male_exp <- iso_exp_data$male / (iso_exp_data$female + iso_exp_data$male) 
hist(iso_exp_data$prop_male_exp)

iso_exp_data$total_iso <- 1
total_iso <- summaryBy(total_iso ~ gene, data=iso_exp_data, FUN=sum)

iso_exp_data$male_over <- ifelse(iso_exp_data$prop_male_exp > 0.5, 1, 0)
iso_exp_data$female_over <- ifelse(iso_exp_data$prop_male_exp < 0.5, 1, 0)

number_male_over <- summaryBy(male_over ~ gene, data=iso_exp_data, FUN=sum)
head(number_male_over)
sum(number_male_over$male_over.sum)

number_female_over <- summaryBy(female_over ~ gene, data=iso_exp_data, FUN=sum)
head(number_female_over)
sum(number_female_over$female_over.sum)

#Answer is no, female use just as many isoforms, they're just expressed much less

fpkm <- SwitchListAnalyzed$isoformRepExpression
head(fpkm)
fpkm$gene <- substr(fpkm$isoform_id,1,nchar(fpkm$isoform_id)-3)
tail(both)

males <- both[both$bias =="Mbias",]

fpkms_male <- fpkm[fpkm$gene %in% both$gene_id,]
fpkms_male$female_exp <- (fpkms_male$PC_F1 + fpkms_male$PC_F2 + fpkms_male$PC_F3)/3
fpkms_male$male_exp <- (fpkms_male$PC_M1 + fpkms_male$PC_M2 + fpkms_male$PC_M3)/3
mean(fpkms_male$male_exp)
mean(fpkms_male$female_exp)





