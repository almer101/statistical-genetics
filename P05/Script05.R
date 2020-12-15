#
# MDS exercise
#

rm(list=ls())

load(url("http://www-eio.upc.es/~jan/data/bsg/CHBChr2-2000.rda"))

ls()

X[1:5,1:5]
dim(X)

X <- t(X)
n <- nrow(X)
n

Alleles[1:10]
length(Alleles)

nmis <- function(x) {
  y <- sum(is.na(x))
  return(y)
}

missingpervariant <- apply(X,2,nmis)
missingpervariant
sum(missingpervariant==45)

X <- X[,!(missingpervariant==45)]
Alleles <- Alleles[!(missingpervariant==45)]
length(Alleles)

dim(X)

X[1:5,1:5]

install.packages("HardyWeinberg")
library(HardyWeinberg)

?recode

Y <- recode(X,Alleles)

Y[1:5,1:5]
max(Y,na.rm=TRUE)
min(Y,na.rm=TRUE)

?dist

D <- dist(Y)
class(D)

D <- as.matrix(dist(Y))
class(D)

head(D)

mds.out <- cmdscale(D,k=n-1,eig=TRUE)
attributes(mds.out)
mds.out$eig
mds.out$GOF
dim(mds.out$points)
?cmdscale

mds.out <- cmdscale(D,k=2,eig=TRUE)
attributes(mds.out)

mds.out$GOF
dim(mds.out$points)


?cmdscale
mds.out$eig

X <- mds.out$points[,1:2]

plot(X[,1],X[,2],asp=1,xlab="First principal axis",
ylab="Second principal axis")

Dest <- as.matrix(dist(X))

Dobs.vec <- D[lower.tri(D)]
Dest.vec <- Dest[lower.tri(Dest)]

plot(Dobs.vec,Dest.vec,xlab="Observed",ylab="Fitted")
abline(0,1,col="blue")
cor(Dobs.vec,Dest.vec)

library(MASS)
isoMDS

nmds.out <- isoMDS(D,k=2)

set.seed(105)

k <- 2
yinit <- matrix(runif(k*n),ncol=k)
nmds.out <- isoMDS(D,k=k,y=yinit)
s2 <- nmds.out$stress

k <- 1
yinit <- matrix(runif(k*n),ncol=k)
nmds.out <- isoMDS(D,k=k,y=yinit)
s1 <- nmds.out$stress

k <- 3
yinit <- matrix(runif(k*n),ncol=k)
nmds.out <- isoMDS(D,k=k,y=yinit)
s3 <- nmds.out$stress

k <- 4
yinit <- matrix(runif(k*n),ncol=k)
nmds.out <- isoMDS(D,k=k,y=yinit)
s4 <- nmds.out$stress

k <- 5
yinit <- matrix(runif(k*n),ncol=k)
nmds.out <- isoMDS(D,k=k,y=yinit)
s5 <- nmds.out$stress

stress <- c(s1,s2,s3,s4,s5)
plot(1:5,stress,type="b")
















Z <- nmds.out$points[,1:2]
plot(Z[,1],Z[,2],asp=1)

R <- cbind(X,Z)
colnames(R) <- c("MDS-1","MDS-2","NMDS-1","NMDS-2")

pairs(R)

round(cor(R),digits=2)


round(cor(R),digits=3)


