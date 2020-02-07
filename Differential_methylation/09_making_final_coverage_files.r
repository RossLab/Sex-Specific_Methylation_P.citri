library(readr)

file.list = list.files(("./"),pattern="*_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", col_names=F)
}

samples <- lapply(file.list, read_file1)

for(i in seq_along(samples)){
    colnames(samples[[i]]) <- c("chr", "cpg", "count_c", "count_t")
    samples[[i]]$total_coverage <- samples[[i]]$count_c + samples[[i]]$count_t
    samples[[i]] <- samples[[i]][,-4]
    final_file <- samples[[i]]
    myfile <- file.path("./", paste0(i,"_","final_coverage.txt"))
    write.table(final_file, file=myfile, quote=F, sep="\t", row.names=F)
}

# NOTE: will need to rename files with corresponding sample name, should improve loop above for this.
# F1
# F2
# F3
# F4
# F5
# M1
# M3
# M4
# M5
