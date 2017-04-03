rm(list=ls())
library(Biostrings)

nRegion <- 600
library(stringr)

pwm <- matrix(as.integer(c(1,5,9,13,0,1,2,3)),byrow=FALSE,nrow=4)
rownames(pwm) <- c("A","C","G","T")

#----------------------- read table ------------------------------------

dinuc <- read.table("table.txt",header=TRUE,sep="\t")
names <- dinuc[,2]

#----------------------- read fasta ------------------------------------

ARF_ER7 <- readDNAStringSet("ER7_unbound.fasta")
ARF_ER8 <- readDNAStringSet("ER8_unbound.fasta")
#ARF <- narrow(ARF,start=35,width=30)
reg_size_ER7<- mean(width(ARF_ER7))
reg_size_ER8<- mean(width(ARF_ER8))

#------------------------- compute  Dinucleotide  -------------------------

seq_ER7 <- as.character(ARF_ER7)
seq_ER8 <- as.character(ARF_ER8)


scores_ER7<- as.data.frame(lapply(seq_ER7,FUN=PWMscoreStartingAt,pwm=pwm,starting.at=(1:(reg_size_ER7-1))))
dim_score_ER7<- dim(scores_ER7)
scores_ER7 <- matrix(as.vector(unlist(scores_ER7)),nrow=dim_score_ER7[1])


scores_ER8<- as.data.frame(lapply(seq_ER8,FUN=PWMscoreStartingAt,pwm=pwm,starting.at=(1:(reg_size_ER8-1))))
dim_score_ER8<- dim(scores_ER8)
scores_ER8 <- matrix(as.vector(unlist(scores_ER8)),nrow=dim_score_ER8[1])



tableau_ER7<- matrix(0,125,reg_size_ER7-1)
rownames(tableau_ER7) <- dinuc[,2]
for (i in 1:(dim_score_ER7[2]))
{
    for (j in 1:(dim_score_ER7[1]))
    {
        tableau_ER7[,j] <- tableau_ER7[,j] + dinuc[,scores_ER7[j,i]+2]
    }
}
tableau_ER7<- tableau_ER7/dim_score_ER7[2]



tableau_ER8<- matrix(0,125,reg_size_ER8-1)
rownames(tableau_ER8) <- dinuc[,2]
for (i in 1:(dim_score_ER8[2]))
{
    for (j in 1:(dim_score_ER8[1]))
    {
        tableau_ER8[,j] <- tableau_ER8[,j] + dinuc[,scores_ER8[j,i]+2]
    }
}
tableau_ER8<- tableau_ER8/dim_score_ER8[2]


tableau_ER7<- tableau_ER7[!sapply(names,FUN=str_detect,pattern="RNA"),]
tableau_ER8<- tableau_ER8[!sapply(names,FUN=str_detect,pattern="RNA"),]

#------------------------print tables ER7/ER8 separately-----------------------


write.table(tableau_ER7,"property_ARF2_ER7_unbound_brut.csv", sep="\t",col.names=FALSE,quote=FALSE)
write.table(tableau_ER8,"property_ARF2_ER8_unbound_brut.csv", sep="\t",col.names=FALSE,quote=FALSE)

#------------------------print combined tables --------------------------------


tableau_ER7_temp <- cbind(tableau_ER7[,1:49],rep(0,dim(tableau_ER7)[1]))
tableau_ER7 <- cbind (tableau_ER7_temp,tableau_ER7[,50:dim(tableau_ER7)[2]])
tableau_ER8[,50] <- tableau_ER8[,50]*2
tableau <- (tableau_ER8*length(ARF_ER8)+tableau_ER7*length(ARF_ER7))/(length(ARF_ER7)+length(ARF_ER8))




## png("heatmap_ARF2_ER7.png",1800,1200)
## tab_heatmap <- heatmap(tableau_ER7, Rowv=NA, Colv=NA, col = rainbow(256,start=3/6,end=4/6), margins=c(5,10))
## dev.off()

## png("heatmap_ARF2_ER8.png",1800,1200)
## tab_heatmap <- heatmap(tableau_ER8, Rowv=NA, Colv=NA, col = rainbow(256,start=3/6,end=4/6), margins=c(5,10))
## dev.off()

png("heatmap_ARF2_unbound.png",1800,1200)
tab_heatmap <- heatmap(tableau, Rowv=NA, Colv=NA, col = rainbow(256,start=3/6,end=4/6), margins=c(5,10))
dev.off()


write.table(tableau,"property_ARF2_unbound_brut.csv", sep="\t",col.names=FALSE,quote=FALSE)

