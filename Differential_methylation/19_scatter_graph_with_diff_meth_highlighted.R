## -------------------------------------------------------------------------
# Making scatter with diff meth genes highlighted
## -------------------------------------------------------------------------
setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/methylation")
library(readr)
library(ggplot2)
library(reshape2)

all_with_meth <- read_delim("weightedMeth_exons_promotors_only.txt", 
                                                "\t", escape_double = FALSE, trim_ws = TRUE)

diff_meth_exon_geneIDs <- read_csv("diff_meth_gene_lists/diff_meth_exon_geneIDs.txt", 
                                   col_names = FALSE)
colnames(diff_meth_exon_geneIDs) <- "gene_id"
diff_meth_promotor_geneIDs <- read_csv("diff_meth_gene_lists/diff_meth_promotor_geneIDs.txt", 
                                       col_names = FALSE)
colnames(diff_meth_promotor_geneIDs) <- "gene_id"
## -------------------------------------------------------------------------

all_with_meth_wide <- dcast(all_with_meth, feature + gene_id ~ origin, value.var="weightedMeth.mean")

exon_data <- all_with_meth_wide[all_with_meth_wide$feature=="exon_first3" ,]
prom_data <- all_with_meth_wide[all_with_meth_wide$feature=="promotors_2000bp" ,]

exon_data$diff <- "no"
exon_data$diff[exon_data$gene_id %in% diff_meth_exon_geneIDs$gene_id] <- "yes"
#write.table(exon_data, file="all_exon_information.txt", sep="\t",
 #            quote = F, col.names = T, row.names = F)

prom_data$diff <- "no"
prom_data$diff[prom_data$gene_id %in% diff_meth_promotor_geneIDs$gene_id] <- "yes"
#write.table(prom_data, file="all_prom_information.txt", sep="\t",
 #           quote = F, col.names = T, row.names = F)

## -------------------------------------------------------------------------

ggplot(prom_data, aes(x=female, y=male))+
  geom_point(aes(colour=diff), size=2.25)+
  ylim(0,1.0)+
  geom_abline(intercept = 0)+
  geom_smooth(method = "loess", level = 0.99)+
  theme_bw()+
  xlab("Female Weighted Methylation")+
  ylab("Male Weighted Methylation")+
  ggtitle("Promotors")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title=element_text(size = 22),
        legend.position = "none")+
  scale_colour_manual(breaks = c("no","yes"),
                    values=c("black","red"))+
  annotate("text", x = 0.2, y = 0.9, size=8,label = "Spearman's rho:\n0.64")

cor.test(prom_data$female, prom_data$male,  method = "spearman")

ggplot(exon_data, aes(x=female, y=male))+
  geom_point(aes(colour=diff), size=2.25)+
  ylim(0,1.0)+
  geom_abline(intercept = 0)+
  geom_smooth(method = "loess", level = 0.99)+
  theme_bw()+
  xlab("Female Weighted Methylation")+
  ylab("Male Weighted Methylation")+
  ggtitle("Exons 1-3")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title=element_text(size = 22),
        legend.position = "none")+
  scale_colour_manual(breaks = c("no","yes"),
                      values=c("black","red"))+
  annotate("text", x = 0.2, y = 0.9, size=8,label = "Spearman's rho:\n0.57")

cor.test(exon_data$female, exon_data$male,  method = "spearman")


# Mealybugs why you be so weird ... not sure what to make of this ...
# I'm happy the stats of calling a feature diff methylated are good





