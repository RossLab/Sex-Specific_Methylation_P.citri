library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(sqldf)

file.list <- list("trim_F1_deduplicated_sorted.bam", "trim_F2_deduplicated_sorted.bam",
"trim_F3_deduplicated_sorted.bam", "trim_F4_deduplicated_sorted.bam",
"trim_F5_deduplicated_sorted.bam", "trim_M1_deduplicated_sorted.bam",
"trim_M3_deduplicated_sorted.bam", "trim_M4_deduplicated_sorted.bam", 
"trim_M5_deduplicated_sorted.bam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("F1", "F2", "F3", "F4", "F5", "M1", "M3", "M4", "M5"),
				treatment = c(0,0,0,0,0,1,1,1,1),
				assembly="Pcitri_v0", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)
