rm(list=ls())
library(Biostrings)
library(plotly)
library(stringr)
library(gplots)
library(RColorBrewer)
library(Cairo)


pwm <- matrix(c(1,5,9,13,0,1,2,3),byrow=FALSE,nrow=4)
rownames(pwm) <- c("A","C","G","T")

#----------------------- read table ------------------------------------

dinuc <- read.table("table.txt",header=TRUE,sep="\t")
dinuc <- dinuc[!sapply(dinuc[,2],FUN=str_detect,pattern="RNA"),]

#----------------------- read fasta ------------------------------------

ARF <- readDNAStringSet("ARF5_bound.fasta")
ARF_neg <- readDNAStringSet("ARF5_unbound.fasta")
reg_size <- mean(width(ARF))
reg_size_neg<- mean(width(ARF_neg))

#------------------------- compute  Dinucleotide  -------------------------

#------------------------------------pos----------------------------

seq <- as.character(ARF)

scores <- as.data.frame(lapply(seq,FUN=PWMscoreStartingAt,pwm=pwm,starting.at=(1:(reg_size-1))))
dim_score<- dim(scores)
scores <- matrix(as.vector(unlist(scores)),nrow=dim_score[1])


M <- list()
tableau <- matrix(0,dim(dinuc)[1],reg_size-1)
rownames(tableau) <- dinuc[,2]
for (i in 1:(dim_score[2]))
{
    for (j in 1:(dim_score[1]))
    {
        tableau[,j] <- dinuc[,scores[j,i]+2]
    }
    M[[i]] <- tableau
    tableau <- matrix(0,dim(dinuc)[1],reg_size-1)
    rownames(tableau) <- dinuc[,2]
}

#-------------------------------------neg----------------------------

seq_neg <- as.character(ARF_neg)

scores_neg <- as.data.frame(lapply(seq_neg,FUN=PWMscoreStartingAt,pwm=pwm,starting.at=(1:(reg_size_neg-1))))
dim_score_neg<- dim(scores_neg)
scores_neg<- matrix(as.vector(unlist(scores_neg)),nrow=dim_score_neg[1])


M_neg<- list()
tableau_neg <- matrix(0,dim(dinuc)[1],reg_size_neg-1)
rownames(tableau_neg) <- dinuc[,2]
for (i in 1:(dim_score_neg[2]))
{
    for (j in 1:(dim_score_neg[1]))
    {
        tableau_neg[,j] <- dinuc[,scores_neg[j,i]+2]
    }
    M_neg[[i]] <- tableau_neg
    tableau_neg <- matrix(0,dim(dinuc)[1],reg_size_neg-1)
    rownames(tableau_neg) <- dinuc[,2]
}

#-------------------------------------compute heatmap----------------------------

## unbound <- read.csv("property_ARF5_unbound_brut.csv",header=FALSE,sep="\t")
## unbound <- data.matrix(unbound[,-1])
## moy <- apply(FUN=mean,unbound,1)


p_value <- matrix(0,dim(dinuc)[1],reg_size-1)
rownames(p_value) <- dinuc[,2]
colnames(p_value) <- c(1:(reg_size-1))
for (i in 1:(dim(p_value)[1]))
{
    for (j in 1:(dim(p_value)[2]))
    {
        values <- NULL
        values_neg <- NULL
        for (k in 1:length(M))
        {
            values[k] <- M[[k]][i,j]
            values_neg[k] <- M_neg[[k]][i,j]
        }
        p_value[i,j] <- t.test(x=values,y=values_neg)$p.value
    }
    print(i)
}

write.table(p_value,"p_values_ARF5.csv",sep="\t",quote=FALSE)
