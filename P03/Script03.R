#
# Exercise on LD
#

#install.packages("genetics")
#install.packages("HardyWeinberg")
#install.packages("LDheatmap")

library(genetics)
library(HardyWeinberg)
library(LDheatmap)

rm(list=ls())

load("./P03/CHBChr2-2000.rda")
ls()
Alleles

class(X)
X[1:10,1:10]

X <- t(X)
dim(X)
dim(Alleles)

SNP12 <- X[,12]
SNP13 <- X[,13]
SNP1000 <- X[,1000]
SNP12

table(SNP12)
table(SNP13)
table(SNP13,useNA="always")


SNP12g <- genotype(SNP12,sep="")
SNP13g <- genotype(SNP13,sep="")
SNP1000g <- genotype(SNP1000,sep="")

summary(SNP12g)
summary(SNP13g)
summary(SNP1000g)

res <- LD(SNP12g,SNP13g)

table(SNP12,SNP1000)

res <- LD(SNP12g,SNP1000g)
res

attributes(res)

chisq <- res$"X^2"
chisq

pchisq(chisq,df=1,lower.tail=FALSE)

n <- nrow(X)
n

nmis <- function(x) {
  y <- sum(is.na(x))
  return(y)
}

table(X[,1],useNA="always")

nmis(X[,1])

dim(X)

nmissingpersnp <- apply(X,2,nmis)
nmissingpersnp

X <- X[,nmissingpersnp==0]
dim(X)

sum(is.na(X))

X100 <- X[,1:100]

RES <- data.frame(genotype(X100[,1],sep=""))
RES

for(i in 2:ncol(X100)) {
   snp <- genotype(X100[,i],sep="")
   RES <- cbind(RES,snp)
}

dim(RES)
top(RES)

output <- LD(RES)
attributes(output)

Dm <- output$D
Dp <- output$"D'"
R2 <- output$"R^2"
X2 <- output$"X^2"

Dm[1:10,1:10]

Dm <- Dm[upper.tri(Dm)]
Dp <- Dp[upper.tri(Dp)]
R2 <- R2[upper.tri(R2)]
X2 <- X2[upper.tri(X2)]

Dm
pairs(cbind(Dm,Dp,R2,X2))

?LDheatmap
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")

LDheatmap(RES,LDmeasure="D'",color=rgb.palette(18))

LDheatmap(RES,LDmeasure="r",color=rgb.palette(18))








