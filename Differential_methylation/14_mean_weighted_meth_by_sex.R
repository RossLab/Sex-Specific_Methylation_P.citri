## -------------------------------------------------------------------------
# Take average weighted methylation level of feature across bio replicates
## -------------------------------------------------------------------------

library(sqldf)
library(readr)
library(doBy)
library(dplyr)


# Make one file covering all samples
file.list = list.files(("./"),pattern="*with_weighted_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)


# Make one dataframe for each population
females <- samples[1:5]
males <- samples[6:9]

for(i in seq_along(females)){
  females[[i]]$origin <- "female"
}
females_all <- as.data.frame(bind_rows(females))
females_merged <- summaryBy(weightedMeth ~ scaffold + feature + gene_id + start + end +
                              cpg_count + origin, data = females_all, FUN=mean)

for(i in seq_along(males)){
  males[[i]]$origin <- "male"
}
males_all <- as.data.frame(bind_rows(males))
males_merged <- summaryBy(weightedMeth ~ scaffold + feature + gene_id + start + end +
                            cpg_count + origin, data = males_all, FUN=mean)



all_data <- rbind(females_merged, males_merged)
write.table(all_data, file="PCITRI_weighted_meth_annotation_by_sex.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")



