##-------------------------------------------------------------------------
# Re-scaling bismark m-bias plots for species with low methylation
##-------------------------------------------------------------------------
setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation/m-bias")
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)
library(readr)


# Before reading in samples needed to take the M-bias .txt output from bismark_methylation_extractor
# open file in excel and remove the headers throughout, add a column for read 1 or 2 and 
# remove the non-CpG context data
temp = list.files(pattern="*.txt")
read_file <- function(x){
  read_delim(x, "\t", escape_double = FALSE, trim_ws = TRUE)}
files = lapply(temp, read_file)
names(files) <- c("F1","F2","F3","F4","F5","M1","M3","M4","M5")

for (i in seq_along(files)){
  files[[i]]$read <- as.factor(files[[i]]$read)
}

# Plot just read 1
for (i in seq_along(files)){
  files[[i]] <- subset(files[[i]], read ==1)
}


# Making the plots
BGI_R2 <- lapply(names(files), function(d) ggplot(data = files[[d]], 
                                            aes(x=position, y=`% methylation`, group=read, 
                                                colour=read))+
                                            geom_line()+
                                            ylab("Percentage Methylation")+
                                            xlab("Read Position")+
                                            ylim(7,11)+
                                            ggtitle(d)+
               theme(axis.text=element_text(size=8),
                     axis.title=element_text(size=8),
                     plot.title = element_text(size=8),
                     legend.title = element_text(size=8),
                     legend.text=element_text(size=8)))

pdf("BGI_read1.pdf")
do.call("grid.arrange", BGI_R2)
dev.off()

