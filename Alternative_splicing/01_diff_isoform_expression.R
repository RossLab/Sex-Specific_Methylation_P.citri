#---------------------------------------------------------
# Alternative Splicing Analysis: Sex-Specific Mealybug
#---------------------------------------------------------

# NOTE: used Stevie's RSEM isoform results

setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/stevie_RSEM_counts")

library(IsoformSwitchAnalyzeR)
library(readr)

file.list <- list.files("./", pattern = "*.results")

quant_data <- importIsoformExpression(sampleVector = file.list)

myDesign <- data.frame(
  sampleID = colnames(quant_data$abundance)[-1],
  condition = c("female","female","female","male","male","male"))

aSwitchList <- importRdata(
  isoformCountMatrix   = quant_data$counts,
  isoformRepExpression = quant_data$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "/Users/holliemarshall/Dropbox/Edinburgh/Sex-specific-mealybugs/genome_annotation/useful_files/annotation_files/PCITRI.assembly.v0.braker.planococcus_citri.gt_edit.gff3",
  showProgress = TRUE)

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 3,
  removeSingleIsoformGenes = TRUE) 
# 1235 isoforms left (as vast majortity of genes not annotaed with alternate transcripts for P.citri)
# 28586 transcripts (93.13%) are single isoform genes anyway
# 847 transcripts filtered on expression levels 
# NOTE: if left expression and isoform cut off at 1 and 0, does't make big difference

# Analyse remaining isoforms with DEXSeq for differntial usage
SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE)

extractSwitchSummary(SwitchListAnalyzed)
#Comparison nrIsoforms nrGenes
#1 female vs male        423      209
# 423 isoforms from 209 genes differentially used


switchingIso <- extractTopSwitches( 
  SwitchListAnalyzed, 
  filterForConsequences = F, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = T,    # when FALSE isoforms are returned
  sortByQvals = TRUE)

write.table(file="list_diff_iso_genes.txt",switchingIso$gene_id,
            sep = "\t", quote = F, col.names = T, row.names = F)


# Find out the top two isoforms by significant qvalue
top_switches <- extractTopSwitches(
  SwitchListAnalyzed, 
  filterForConsequences = F, 
  n = 2, 
  sortByQvals = TRUE)

# Make a plot of one of the isoforms
switchPlot(SwitchListAnalyzed, gene = 'g19348')

# Make a plot of the top 10 isoforms (automatically outputs as pdf)
switchPlotTopSwitches(
  switchAnalyzeRlist = SwitchListAnalyzed, 
  n = 10,
  filterForConsequences = FALSE, 
  splitFunctionalConsequences = F)





