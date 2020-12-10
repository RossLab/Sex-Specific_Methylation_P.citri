#---------------------------------------------------------
# Expression levels of DNMT1 gene
#---------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Sex-specific-mealybugs/transcription")

library(readr)
library(reshape2)
library(ggplot2)
#---------------------------------------------------------

P_citri_FPKMs_trimmed <- read_csv("P_citri_FPKMs_trimmed.csv")

# g14868 is the DNMT1 (from Stevie's thesis)

head(P_citri_FPKMs_trimmed)

dnmt1 <- P_citri_FPKMs_trimmed[P_citri_FPKMs_trimmed$Gene=="g14868",]
dnmt1 <- melt(dnmt1)
dnmt1$Sex <- c("Female","Female","Female","Male","Male","Male")

ggplot(dnmt1, aes(x=Sex, y=value, fill=Sex))+
  geom_boxplot()+
  geom_jitter(size=5)+
  theme_bw()+
  xlab("Sex")+
  ylab("DNMT1 Expression Level (FPKM)")+
  scale_fill_manual("",breaks=c("Female", "Male"),
                    values = c("grey75", "grey45","grey14"))+
  theme_bw()+
  scale_fill_manual(limits=c("Female", "Male"),
                   values = c("pink1","steelblue1"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none")




