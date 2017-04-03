rm(list=ls())
library(Biostrings)
library(plotly)
nRegion <- 600
library(stringr)

pwm <- matrix(c(1,5,9,13,0,1,2,3),byrow=FALSE,nrow=4)
rownames(pwm) <- c("A","C","G","T")

#----------------------- read table ------------------------------------

dinuc <- read.table("table.txt",header=TRUE,sep="\t")

#----------------------- read fasta ------------------------------------

ARF <- readDNAStringSet("ARF5_DR26.fas")
#ARF <- narrow(ARF,start=35,width=30)
reg_size <- mean(width(ARF))

#------------------------- compute  Dinucleotide  -------------------------

seq <- as.character(ARF)


scores <- as.data.frame(lapply(seq,FUN=PWMscoreStartingAt,pwm=pwm,starting.at=(1:(reg_size-1))))
dim_score<- dim(scores)
scores <- matrix(as.vector(unlist(scores)),nrow=dim_score[1])

tableau <- matrix(0,125,reg_size-1)
rownames(tableau) <- dinuc[,2]
for (i in 1:(dim_score[2]))
{
    for (j in 1:(dim_score[1]))
    {
        tableau[,j] <- tableau[,j] + dinuc[,scores[j,i]+2]
    }
}
tableau <- tableau/dim_score[2]



## min_l <- apply(tableau,FUN=min,1)
## tableau <- tableau-min_l
## max_l<- apply(tableau,FUN=max,1)
## somme<- apply(tableau,FUN=sum,1)
## tableau <- abs(tableau)/abs(somme)
names <- rownames(tableau)
tableau <- tableau[!sapply(names,FUN=str_detect,pattern="RNA"),]


## png("heatmap_ARF5_max.png",1800,1200)
## tab_heatmap <- heatmap(tableau, Rowv=NA, Colv=NA, col = rainbow(256,start=3/6,end=4/6), margins=c(5,10))
## dev.off()

write.table(tableau,"property_ARF5_DR26_bound_brut.csv", sep="\t",col.names=FALSE,quote=FALSE)




##############################"Normalisation#############################################




bound <- read.csv("property_ARF5_bound_brut.csv",header=FALSE,sep="\t")
names_bound <- bound[,1]
bound <- data.matrix(bound[,-1])
rownames(bound) <- names_bound
colnames(bound) <- c(1:99)

unbound <- read.table("property_ARF5_unbound_brut.csv",header=FALSE,sep="\t")
names_unbound <- unbound[,1] 
unbound <- data.matrix(unbound[,-1])
rownames(unbound) <- names_unbound
colnames(unbound) <- c(1:99)



tab_2 <- cbind(bound)#,unbound)
tab_shift <- tab_2
mini <- apply(tab_2,FUN=min,1)
tab_2 <- tab_2-mini
maxi <- apply(tab_2,FUN=max,1)

bound <- (bound-mini)/maxi
unbound <- (unbound-mini)/maxi



write.table(bound,"property_ARF5_DR26_bound_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)
write.table(unbound,"property_ARF5_unbound_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)


## tab_shift_m <- matrix(c(tab_shift),nrow=2,byrow=TRUE)


## png("heatmap_Adrien_max.png",1800,1200)
## tab_heatmap <- heatmap(tab_shift_m, Rowv=NA, Colv=NA, col = rainbow(256,start=3/6,end=4/6), margins=c(5,10))
## dev.off()
## write.table(tab_shift_m_2,"property_ARF5_Adrien.csv", sep="\t",col.names=FALSE,quote=FALSE)

## ref <- bound-unbound
## colnames(ref) <- c(1:99)

## p <- plot_ly(z = tab_shift, type = "heatmap")

## dev.off()

