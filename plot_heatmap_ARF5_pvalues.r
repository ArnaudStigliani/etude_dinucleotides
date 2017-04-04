rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Cairo)


pvalues <- read.csv("p_values_ARF5.csv",header=FALSE,sep="\t")
names_pvalues <- pvalues[,1]
pvalues <- data.matrix(pvalues[,-1])
rownames(pvalues) <- names_pvalues
colnames(pvalues) <- c(1:99)
pvalues <- data.matrix(pvalues)

pvalues <-  log(pvalues+10^(-200),10)


my_palette <- colorRampPalette(c("red","blue"))(n = 400)
my_palette <- c(my_palette,colorRampPalette(c("white"))(n = 5))




 Cairo(width = 1300, height = 900 , file="ARF5_pvalues_part1.png", type="png", pointsize=18, 
       bg = "white", canvas = "white", units = "px", dpi = "auto")

heatmap.2(pvalues[1:55,],
  ## cellnote = bound,  # same data set for cell labels
  ## main = "Correlation", # heat map title
  ## notecol="black",      # change font color of cell labels to black
  density.info="none",
  trace="none",                                      # turns off density plot inside color legend
  ## trace="none",
  Rowv=FALSE,                      # turns off trace lines inside the heat map
  symkey=FALSE,
  symbreaks=FALSE,
  margins =c(5,15),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
 #  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  Colv="NA",
  keysize=1
)            # turn off column clustering


dev.off()               # close the PNG device

 Cairo(width = 1300, height = 900, file="ARF5_pvalues_part2.png", type="png", pointsize=18, 
       bg = "white", canvas = "white", units = "px", dpi = "auto")
heatmap.2(pvalues[56:110,],
  ## cellnote = bound,  # same data set for cell labels
  ## main = "Correlation", # heat map title
  ## notecol="black",      # change font color of cell labels to black
  density.info="none",
  trace="none",                                      # turns off density plot inside color legend
  ## trace="none",
  Rowv=FALSE,                      # turns off trace lines inside the heat map
  symkey=FALSE,
  symbreaks=FALSE,
  margins =c(5,15),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  ## #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  Colv="NA",
  keysize=1
  )            # turn off column clustering
dev.off()               # close the PNG device
