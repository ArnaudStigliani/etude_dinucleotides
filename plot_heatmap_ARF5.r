rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Cairo)

bound <- read.csv("property_ARF5_bound_brut.csv",header=FALSE,sep="\t")
names_bound <- bound[,1]
bound <- data.matrix(bound[,-1])
rownames(bound) <- names_bound
colnames(bound) <- c(1:99)
bound <- data.matrix(bound)

bound <- read.csv("property_ARF5_unbound_brut.csv",header=FALSE,sep="\t")
names_bound <- bound[,1]
bound <- data.matrix(bound[,-1])
rownames(bound) <- names_bound
colnames(bound) <- c(1:99)
bound <- data.matrix(bound)


bound <- t(scale(t(bound)))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

 Cairo(width = 1300, height = 900 , file="ARF5_part1.png", type="png", pointsize=18, 
       bg = "white", canvas = "white", units = "px", dpi = "auto")
heatmap.2(boundl,
  ## cellnote = bound,  # same data set for cell labels
  ## main = "Correlation", # heat map title
  ## notecol="black",      # change font color of cell labels to black
  density.info="none",
  trace="none",                                      # turns off density plot inside color legend
  ## trace="none",
  Rowv=FALSE,                      # turns off trace lines inside the heat map
  margins =c(5,15),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  ## #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  Colv="NA",
  keysize=1
  )            # turn off column clustering
dev.off()               # close the PNG device

