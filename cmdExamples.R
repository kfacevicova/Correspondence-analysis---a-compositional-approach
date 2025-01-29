source("progs.R")
library(rospca)

# Plots for paper:

# Spanish health data: Fig. 1 and 2:
d <- read.csv("health.csv")
x <- d[,-1]
dimnames(x) <- list(d[,1],c("Very good","Good","Regular","Bad","Very bad"))
res <- coda_ca(x)
w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
resw <- coda_ca(x, w = w.ind, weighted=TRUE)
pdf("Spanish1.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(4.3,4.3,3,2))
plotca(res,weighted=FALSE)
plotca(resw,weighted=TRUE)
dev.off()

set.seed(124)
r <- bootst(x,B=10000)
pdf("Spanish2.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(4.3,4.3,1,1))
boxplot(r$ang,r$angw,ylab="Angle",names = c("Unweighted","Weighted"),cex.lab=1.2)
boxplot((r$ang-r$angw)~r$minnb,notch=TRUE,ylab="Angle unweighted minus weighted",
        xlab="Smallest value in table",cex.lab=1.2)
abline(h=0,lty=3)
dev.off()

##########################################################################################
# Fig. 3 and 4:

# Stores data:
d <- read.csv("stores.csv")
x <- d[,-1]
dimnames(x) <- list(d[,1],c("16-24","25-34","35-49","50+"))
res <- coda_ca(x)
sdU.sto <- res$simp # stores
w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
resw <- coda_ca(x, w = w.ind, weighted = TRUE)
sdW.sto <- resw$simpw # stores

set.seed(123)
r.sto <- bootst(x,B=1000)

# Spanish health data aggregated
d <- read.csv("health.csv")
x <- d[,-1]
x[,4] <- x[,4]+x[,5]
x <- x[,1:4]
dimnames(x) <- list(d[,1],c("Very good","Good","Regular","Bad"))
res <- coda_ca(x)
sdU.h2 <- res$simp # health2
w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
resw <- coda_ca(x, w = w.ind, weighted = TRUE)
sdW.h2 <- resw$simp # health2

set.seed(123)
r.h2 <- bootst(x,B=1000)

# News data
d <- read.csv("news.csv") # News interest in Europe (chapter 19) of Joint Correspondence Analysis
x <- d[-c(1:18),-1]
dimnames(x)[[1]] <- d[-c(1:18),1]
res <- coda_ca(x)
sdU.news <- res$simp # news
w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
resw <- coda_ca(x, w = w.ind, weighted = TRUE)
sdW.news <- resw$simp # news

set.seed(123)
r.news <- bootst(x,B=1000)

# Roman glass cups:
library(easyCODA)
data(cups)
x <- cups*100
res <- coda_ca(x)
sdU.cups <- res$simp # cups
w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
resw <- coda_ca(x, w = w.ind, weighted = TRUE)
sdW.cups <- resw$simp # cups

set.seed(123)
r.cups <- bootst(x,B=1000)

# fish data
library(easyCODA)
data(fish)
x <- fish[,-c(1:3)]*100
res <- coda_ca(x)
sdU.fish <- res$simp # Fish
w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
resw <- coda_ca(x, w = w.ind, weighted = TRUE)
sdW.fish <- resw$simp # Fish

set.seed(123)
r.mor <- bootst(x,B=1000)

# Galton data:
galton=read.table("galton_data_categories.txt")
rownames(galton)=c("Above","72.5","71.5","70.5","69.5","68.5","67.5","66.5","65.5","64.5","Below")
colnames(galton)=c("Below","62.2","63.2","64.2","65.2","66.2","67.2","68.2","69.2","70.2","71.2","72.2","73.2","Above")
#reduced data
x.reduced <- rbind(galton[1,] + galton[2,], galton[3:9,], galton[10,] + galton[11,])
x.reduced <- cbind(x.reduced[,1] + x.reduced[,2], x.reduced[,3:12], x.reduced[,13] + x.reduced[,14])
colnames(x.reduced)[c(1, 12)] <- c("Below", "Above")
rownames(x.reduced)[c(1, 9)] <- c("Above", "Below")
x <- x.reduced
#imputation
x[x==0]=2/3
res <- coda_ca(x)
sdU.gal <- res$simp # Galton
w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
resw <- coda_ca(x, w = w.ind, weighted = TRUE)
sdW.gal <- resw$simp # Galton

set.seed(123)
r.gal <- bootst(x,B=1000)


pdf("boxplot_stores.pdf",width=2,height=4.5)
par(mar=c(0.5,4.2,2,0.5))
boxplot(r.sto$ang-r.sto$angw,ylab="Angle unweighted minus weighted",cex.lab=1.2,main="Stores")
abline(h=0,lty=3)
dev.off()

pdf("boxplot_health2.pdf",width=1.5,height=4.5)
par(mar=c(0.5,2,2,0.5))
boxplot(r.h2$ang-r.h2$angw,main="Health2",ylab="")
abline(h=0,lty=3)
dev.off()

pdf("boxplot_cups.pdf",width=1.5,height=4.5)
par(mar=c(0.5,2,2,0.5))
boxplot(r.cups$ang-r.cups$angw,ylab="",main="Cups")
abline(h=0,lty=3)
dev.off()

pdf("boxplot_news.pdf",width=1.5,height=4.5)
par(mar=c(0.5,2,2,0.5))
boxplot(r.news$ang-r.news$angw,ylab="",main="News")
abline(h=0,lty=3)
dev.off()

pdf("boxplot_galton.pdf",width=1.5,height=4.5)
par(mar=c(0.5,2,2,0.5))
boxplot(r.gal$ang-r.gal$angw,ylab="",main="Galton")
abline(h=0,lty=3)
dev.off()

pdf("boxplot_fish.pdf",width=1.5,height=4.5)
par(mar=c(0.5,2,2,0.5))
boxplot(r.mor$ang-r.mor$angw,ylab="",main="Fish")
abline(h=0,lty=3)
dev.off()

############################################
pdf("boxplotsd_stores.pdf",width=2,height=4.5)
par(mar=c(2,4.2,2,0.5))
boxplot(r.sto$simp,r.sto$simpw,ylab="Simplicial deviance",ylim=c(0,1),
        cex.lab=1.2,main="Stores",names=c("U","W"))
segments(0.5,sdU.sto,1.5,sdU.sto,col=6,lty=1,lwd=2)
segments(1.5,sdW.sto,2.5,sdW.sto,col=6,lty=1,lwd=2)
dev.off()

pdf("boxplotsd_health2.pdf",width=1.5,height=4.5)
par(mar=c(2,2,2,0.5))
boxplot(r.h2$simp,r.h2$simpw,ylab="",ylim=c(0,1),
        cex.lab=1.2,main="Health2",names=c("U","W"))
segments(0.5,sdU.h2,1.5,sdU.h2,col=6,lty=1,lwd=2)
segments(1.5,sdW.h2,2.5,sdW.h2,col=6,lty=1,lwd=2)
dev.off()

pdf("boxplotsd_cups.pdf",width=1.5,height=4.5)
par(mar=c(2,2,2,0.5))
boxplot(r.cups$simp,r.cups$simpw,ylab="",ylim=c(0,1),
        cex.lab=1.2,main="Cups",names=c("U","W"))
segments(0.5,sdU.cups,1.5,sdU.cups,col=6,lty=1,lwd=2)
segments(1.5,sdW.cups,2.5,sdW.cups,col=6,lty=1,lwd=2)
dev.off()

pdf("boxplotsd_news.pdf",width=1.5,height=4.5)
par(mar=c(2,2,2,0.5))
boxplot(r.news$simp,r.news$simpw,ylab="",ylim=c(0,1),
        cex.lab=1.2,main="News",names=c("U","W"))
segments(0.5,sdU.news,1.5,sdU.news,col=6,lty=1,lwd=2)
segments(1.5,sdW.news,2.5,sdW.news,col=6,lty=1,lwd=2)
dev.off()

pdf("boxplotsd_galton.pdf",width=1.5,height=4.5)
par(mar=c(2,2,2,0.5))
boxplot(r.gal$simp,r.gal$simpw,ylab="",ylim=c(0,1),
        cex.lab=1.2,main="Galton",names=c("U","W"))
segments(0.5,sdU.gal,1.5,sdU.gal,col=6,lty=1,lwd=2)
segments(1.5,sdW.gal,2.5,sdW.gal,col=6,lty=1,lwd=2)
dev.off()

pdf("boxplotsd_fish.pdf",width=1.5,height=4.5)
par(mar=c(2,2,2,0.5))
boxplot(r.mor$simp,r.mor$simpw,ylab="",ylim=c(0,1),
        cex.lab=1.2,main="Fish",names=c("U","W"))
segments(0.5,sdU.fish,1.5,sdU.fish,col=6,lty=1,lwd=2)
segments(1.5,sdW.fish,2.5,sdW.fish,col=6,lty=1,lwd=2)
dev.off()

