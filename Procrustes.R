#Daniel Shriner
#3/23/20
#implementation of Procrustes analysis using principal coordinate analysis
#comparison of Fst and geographic distance matrices of equal dimensions

dat1 <- read.table("fst.txt",header=FALSE)
dat2 <- read.table("geographic_distances.txt",header=FALSE)
reps <- 1e9
Xloc <- cmdscale(as.dist(dat1), k = nrow(dat1)-1,eig=TRUE,add=TRUE)
Yloc <- cmdscale(as.dist(dat2), k = nrow(dat2)-1,eig=TRUE,add=TRUE)
Xcentroid <- Xloc$points[,c(1:2)]
Ycentroid <- Yloc$points[,c(1:2)]
C <- t(Ycentroid)%*%Xcentroid
Csvd <- svd(C)
D <- Csvd$d
statistic <- sqrt((sum(D)^2)/(sum(diag(t(Xcentroid)%*%Xcentroid))*sum(diag(t(Ycentroid)%*%Ycentroid))))
nulldist <- vector("numeric",reps)
for (i in 1:reps) {
	o <- sample(nrow(dat1))
	Xperm <- Xcentroid[o,]
	C <- t(Ycentroid)%*%Xperm
	Csvd <- svd(C)
	D <- Csvd$d
	nulldist[i] <- sqrt((sum(D)^2)/(sum(diag(t(Xperm)%*%Xperm))*sum(diag(t(Ycentroid)%*%Ycentroid))))
}
length(which(nulldist>statistic))/length(nulldist)
