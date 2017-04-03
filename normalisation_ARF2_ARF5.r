rm(list=ls())

ARF5 <- read.table("property_ARF5_bound_brut.csv",header=FALSE,sep="\t")
names_ARF5 <- ARF5[,1]
ARF5 <- data.matrix(ARF5[,-1])
rownames(ARF5) <- names_ARF5
colnames(ARF5) <- c(1:99)
ARF5_5_prime <- ARF5[,9:43]
ARF5_3_prime <- ARF5[,57:91]




ARF2 <- read.table("property_ARF2_bound_brut.csv",header=FALSE,sep="\t")
names_ARF2 <- ARF2[,1] 
ARF2 <- data.matrix(ARF2[,-1])
rownames(ARF2) <- names_ARF2
colnames(ARF2) <- c(1:99)
rev_ARF2 <- t(apply(FUN=rev,ARF2,1))
ARF2_mean <- (ARF2+rev_ARF2)/2
ARF2_5_prime <- ARF2_mean[,1:34]
ARF2_3_prime <- ARF2_mean[,65:99]


tab_5_prime<- cbind(ARF5_5_prime,ARF2_5_prime)
mini_5_prime <- apply(tab_5_prime,FUN=min,1)
tab_5_prime <- tab_5_prime - mini_5_prime
maxi_5_prime <- apply(tab_5_prime,FUN=max,1)

ARF5_5_prime <- (ARF5_5_prime -mini_5_prime)/maxi_5_prime
ARF2_5_prime <- (ARF2_5_prime -mini_5_prime)/maxi_5_prime

tab_3_prime<- cbind(ARF5_3_prime,ARF2_3_prime)
mini_3_prime <- apply(tab_3_prime,FUN=min,1)
tab_3_prime <- tab_3_prime - mini_3_prime
maxi_3_prime<- apply(tab_3_prime,FUN=max,1)

ARF5_3_prime <- (ARF5_3_prime - mini_3_prime)/ maxi_3_prime
ARF2_3_prime <- (ARF2_3_prime - mini_3_prime)/ maxi_3_prime


options(scipen=999)

write.table(format(round(ARF5_3_prime,5),nsmall=2),"property_ARF5_3_prime_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)
write.table(format(round(ARF5_5_prime,5),nsmall=2),"property_ARF5_5_prime_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)
write.table(format(round(ARF2_3_prime,5),nsmall=2),"property_ARF2_3_prime_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)
write.table(format(round(ARF2_5_prime,5),nsmall=2),"property_ARF2_5_prime_normalized_-min_sur_max.csv", sep="\t",col.names=FALSE,quote=FALSE)
