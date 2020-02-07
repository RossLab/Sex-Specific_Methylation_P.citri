## -------------------------------------------------------------------------
# Weighted methylation per annotation for each sex
## -------------------------------------------------------------------------

library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*final_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

# Read in gene with start/end and total CpGs per gene
annotation_with_total_cpgs <- read_table2("merged_annotations.txt")

## -------------------------------------------------------------------------
registerDoParallel(cores = 10)

# Calculate weighted meth for each gene for each sample
foreach(i = seq_along(samples)) %dopar% {
  df <- samples[[i]]
  df <- subset(df, total_coverage > 5)
  output <- sqldf("SELECT sample.chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.scaffold,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.cpg_count,
                    annot.feature
                    FROM df AS sample
                    LEFT JOIN annotation_with_total_cpgs AS annot
                    ON sample.chr = annot.scaffold
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$gene_id),]
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ scaffold + feature + gene_id + start + end + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$cpg_count*check$count_c.sum)/(check$cpg_count*check$total_coverage.sum)
  myfile <- file.path("./", paste0(i,"_","with_weighted_meth.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}


## -------------------------------------------------------------------------
