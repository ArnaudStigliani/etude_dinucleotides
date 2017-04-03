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

bound <- t(scale(t(bound)))
## bound <- t(scale(t(cbind(bound[,1:42],bound[,54:99]))))
## bound <- cbind(bound[,1:42],matrix(0,ncol=11,nrow=110),bound[,(55-11):(99-11)])
## bound <- t(scale(t(cbind(bound[,1:38],bound[,50:]))))
## bound <- cbind(bound[,1:38],matrix(0,ncol=12,nrow=110),bound[,(55-12):(99-12)])
#bound <- t(scale(t(cbind(bound[,1:39],bound[,62:99]))))
#bound <- t(scale(t(bound[,47:53])))


boundl <- bound*4
sign <- ifelse(bound<0,-1,1)
boundl <- log(abs(boundl)+1)
boundl <- boundl * sign

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)


col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

 Cairo(width = 1300, height = 900 , file="ARF5_part1.png", type="png", pointsize=18, 
       bg = "white", canvas = "white", units = "px", dpi = "auto")

heatmap.2(boundl[1:55,],
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

 Cairo(width = 1300, height = 900, file="ARF5_part2.png", type="png", pointsize=18, 
       bg = "white", canvas = "white", units = "px", dpi = "auto")

heatmap.2(boundl[56:110,],
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
