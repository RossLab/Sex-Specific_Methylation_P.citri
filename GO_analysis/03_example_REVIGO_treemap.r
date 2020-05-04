
setwd("~/Dropbox/Edinburgh/Sex-specific-mealybugs/GO_analyses")
# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006468","protein phosphorylation",4.137,2.1536,0.681,0.000,"protein phosphorylation"),
c("GO:0044267","cellular protein metabolic process",14.293,1.4165,0.747,0.536,"protein phosphorylation"),
c("GO:0070647","protein modification by small protein conjugation or removal",0.821,1.8839,0.750,0.592,"protein phosphorylation"),
c("GO:0006796","phosphate-containing compound metabolic process",13.110,2.1296,0.815,0.664,"protein phosphorylation"),
c("GO:0009987","cellular process",63.780,1.4735,0.967,0.000,"cellular process"),
c("GO:0042592","homeostatic process",1.661,1.7526,0.608,0.000,"homeostatic process"),
c("GO:0072507","divalent inorganic cation homeostasis",0.111,1.8249,0.539,0.600,"homeostatic process"),
c("GO:0070838","divalent metal ion transport",0.358,2.8927,0.673,0.000,"divalent metal ion transport"),
c("GO:0006887","exocytosis",0.210,1.8249,0.693,0.242,"divalent metal ion transport"),
c("GO:0051641","cellular localization",2.041,1.3449,0.764,0.289,"divalent metal ion transport"),
c("GO:0007010","cytoskeleton organization",0.786,1.8839,0.882,0.041,"cytoskeleton organization"),
c("GO:0017144","drug metabolic process",0.058,1.3249,0.877,0.056,"drug metabolism"),
c("GO:1901564","organonitrogen compound metabolic process",17.886,1.9195,0.889,0.091,"organonitrogen compound metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="./treemaps_and_annotations/treemaps/alt_spliced.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "Significantly Enriched GO Terms, p-val<0.05",
  inflate.labels = T,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 1,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCC00",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "bottom",
  title.legend = "",
  fontsize.labels = c(0,5),
  #force.print.labels=TRUE,
  
  fontsize.legend = 8
  ,palette = "Set3"
)

dev.off()
