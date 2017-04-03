rm(list=ls())

################################## Normalisation ############################

bound <- read.csv("property_ARF2_bound_brut.csv",header=FALSE,sep="\t")
names_bound <- bound[,1]
bound <- data.matrix(bound[,-1])
rownames(bound) <- names_bound
colnames(bound) <- c(1:99)

unbound <- read.table("property_ARF2_unbound_brut.csv",header=FALSE,sep="\t")
names_unbound <- unbound[,1] 
unbound <- data.matrix(unbound[,-1])
rownames(unbound) <- names_unbound
colnames(unbound) <- c(1:99)



tab_2 <- cbind(bound,unbound)
tab_shift <- tab_2
mini <- apply(tab_2,FUN=min,1)
tab_2 <- tab_2-mini
maxi <- apply(tab_2,FUN=max,1)

bound <- (bound-mini)/maxi
unbound <- (unbound-mini)/maxi

options(scipen=999)

write.table(format(round(bound,5),nsmall=5),"property_ARF2_bound_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)
write.table(format(round(unbound,5),nsmall=5),"property_ARF2_unbound_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)

#----------------------------------ER7 and ER8 separated---------------------

rm(list=ls())
################################## Normalisation ############################
bound_ER7<- read.csv("property_ARF2_ER7_bound_brut.csv",header=FALSE,sep="\t")
names <-  bound_ER7[,1]
bound_ER7<- data.matrix(bound_ER7[,-1])
rownames(bound_ER7) <- names
colnames(bound_ER7) <- c(1:98)


bound_ER8<- read.csv("property_ARF2_ER8_bound_brut.csv",header=FALSE,sep="\t")
bound_ER8 <- data.matrix(bound_ER8[,-1])
rownames(bound_ER8) <- names
colnames(bound_ER8) <- c(1:99)

unbound_ER7 <- read.table("property_ARF2_ER7_unbound_brut.csv",header=FALSE,sep="\t")
unbound_ER7 <- data.matrix(unbound_ER7[,-1])
rownames(unbound_ER7) <- names
colnames(unbound_ER7) <- c(1:98)

unbound_ER8 <- read.table("property_ARF2_ER8_unbound_brut.csv",header=FALSE,sep="\t")
unbound_ER8 <- data.matrix(unbound_ER8[,-1])
rownames(unbound_ER8) <- names
colnames(unbound_ER8) <- c(1:99)

tab_tot <- cbind(bound_ER7,unbound_ER7,bound_ER8,unbound_ER8)

mini <- apply(tab_tot,FUN=min,1)
tab_tot <- tab_tot - mini
maxi <- apply(tab_tot,FUN=max,1)

bound_ER7 <- (bound_ER7 - mini ) /maxi
bound_ER8 <- (bound_ER8 - mini ) /maxi
unbound_ER7 <- (unbound_ER7 - mini ) /maxi
unbound_ER8 <- (unbound_ER8 - mini ) /maxi

bound_ER7 <- cbind(bound_ER7[,1:49],rep(0.5,dim(bound_ER7)[1]),bound_ER7[,50:98])
unbound_ER7 <- cbind(unbound_ER7[,1:49],rep(0.5,dim(unbound_ER7)[1]),unbound_ER7[,50:98])


bound <- NULL
unbound <- NULL
names_7_8 <- NULL
for (i in 1:dim(bound_ER7)[1])
{
    bound <- rbind(bound,bound_ER7[i,],bound_ER8[i,])
    names_7_8 <- rbind(names_7_8,as.character(names[i]),as.character(names[i]))
    unbound <- rbind(unbound,unbound_ER7[i,],unbound_ER8[i,])
}
rownames(bound) <- names_7_8 
rownames(unbound) <- names_7_8 

options(scipen=999)


write.table(format(round(bound,5),nsmall=5),"property_ARF2_bound_normalized_-min_sur_max_ER7_8_separated.csv", sep="\t",col.names=FALSE,quote=FALSE)
write.table(format(round(unbound,5),nsmall=5),"property_ARF2_unbound_normalized_-min_sur_max_ER7_8_separated.csv", sep="\t",col.names=FALSE,quote=FALSE)
