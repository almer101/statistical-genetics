#
# Exercise descriptive analysis of genetic markers.
#

#
# First part: SNP Data
#


filename <- url("http://www-eio.upc.es/~jan/Data/bsg/Chromosome1_CHBPopSubset.rda")
load(filename)

install.packages("genetics")
library(genetics)

dim(Ysub)
class(Ysub)
Ysub[1:5,1:5]
n <- nrow(Ysub)
p <- ncol(Ysub)
n
p
Ysub[Ysub=="NN"] <- NA
perc.mis <- 100*sum(is.na(Ysub))/(n*p)
perc.mis

SNP1 <- Ysub[,1]
SNP1
SNP1.g <- genotype(SNP1,sep="")
SNP1.g
summary(SNP1.g)

SNP2 <- Ysub[,2]
SNP2.g <- genotype(SNP2,sep="")
summary(SNP2.g)

SNP3 <- Ysub[,3]
SNP3.g <- genotype(SNP3,sep="")
summary(SNP3.g)

nmis <- function(x) {
  y <- sum(is.na(x))
  return(y)
}

nmis.per.ind <- apply(Ysub,1,nmis)
pmis.per.ind <- 100*nmis.per.ind/p

plot(1:n,pmis.per.ind,xlab="Individual",ylab="Perc. Missing")

nmis.per.snp <- apply(Ysub,2,nmis)
pmis.per.snp <- 100*nmis.per.snp/n

plot(1:p,pmis.per.snp,xlab="SNP",ylab="Perc. Missing")

sum(nmis.per.snp==n)

x <- table(SNP3)
x
sum(x)
sum(x,na.rm=TRUE)
n

pC <- (2*x[1]+x[2])/(2*sum(x,na.rm=TRUE))
pT <- (2*x[3]+x[2])/(2*sum(x,na.rm=TRUE))
pC
pT
pC+pT
summary(SNP3.g)

Y2 <- Ysub[,nmis.per.snp < n]

affirst <- function(x){
  x <- genotype(x,sep="")
  out <- summary(x)
  af1 <- out$allele.freq[1,2] 
  return(af1)
}

affirst(Ysub[,1])
affirst(Ysub[,3])

af.first.allele <- apply(Y2,2,affirst)
hist(af.first.allele)

maf <- function(x){
  x <- genotype(x,sep="")
  out <- summary(x)
  af1 <- min(out$allele.freq[,2],na.rm=TRUE)
  af1[af1==1] <- 0 
  return(af1)
}

maf.per.snp <- apply(Y2,2,maf)
hist(maf.per.snp)

#
# Second part: STR data
#

rm(list=ls())

filename <- url("http://www-eio.upc.es/~jan/Data/bsg/JapanaseSTRs.rda")
load(filename)

ls()

Japanese[1:5,1:10]
X <- Japanese[,6:ncol(Japanese)]
n <- nrow(X)/2
p <- ncol(X)
n
p

sum(X==-9)
X[X==-9] <- NA
sum(is.na(X))

class(X)
STR1 <- X[,1]
table(STR1,useNA="always")

length(unique(STR1))

n.alleles <- function(x) {
  y <- length(unique(x[!is.na(x)]))
  return(y)
}

n.alleles(STR1) # number of alleles
table(STR1) # allele counts

STR1 <- STR1[!is.na(STR1)] 
na <- length(STR1)
na

index.1 <- seq(1,na,2)
index.2 <- seq(2,na,2)

allele.1 <- STR1[index.1]
allele.2 <- STR1[index.2]

allele.1 <- pmin(allele.1,allele.2)
allele.2 <- pmax(allele.1,allele.2)

allele.1
allele.2

individuals <- paste(allele.1,allele.2,sep="/")
individuals
g.counts <- table(individuals) # genotype counts
g.counts
unique(names(g.counts))
names(g.counts)

sum(g.counts)
length(g.counts) # number of genotypes

K <- n.alleles(STR1) # number of alleles
K

n.genotypes <- 0.5*K*(K+1)
n.genotypes

n.alleles.per.STR <- apply(X,2,n.alleles)
barplot(table(n.alleles.per.STR))

