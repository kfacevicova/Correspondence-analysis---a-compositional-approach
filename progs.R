########################
### Correspondence analysis of a compositional table
########################

coda_ca <- function(x, w = NULL, weighted = FALSE)
{
  I = nrow(x)
  J = ncol(x)
  clr_x <- log(x/exp(1/(I*J)*sum(log(x))))
  
  #marginals
  yrow = apply(clr_x, 1, sum)/J
  ycol = apply(clr_x, 2, sum)/I
  #independent table - clr
  yind = yrow*matrix(1, nrow = I, ncol = J) + 
       t(ycol*t(matrix(1, nrow = I, ncol = J)))
  #interaction table - clr
  yint = clr_x - yind
  
  if(is.null(w)){   # unweighted CA
    #singular value decomposition
    yint_dec = svd(yint)
    u = yint_dec$u
    rownames(u) = rownames(x)
    v = yint_dec$v
    rownames(v) = colnames(x)
    
    # explained variances
    exp.variability.int.unweighted <- yint_dec$d^2/sum(yint_dec$d^2)*100 
    exp.variability.orig.unweighted <- yint_dec$d^2/norm(as.matrix(clr_x),type="F")^2*100 
    #simplicial deviance
    simp.deviance <- norm(as.matrix(yint),type="F")^2/norm(as.matrix(clr_x),type="F")^2
    return(list(u = u, v = v, 
                intu = exp.variability.int.unweighted, 
                oriu = exp.variability.orig.unweighted,
                intw = NULL, 
                oriw = NULL,
                simp = simp.deviance))
  } 
  else { # weighted CA
      # weighting is applied on the complete table
      x.w <- x/w
      # clr representation wrt. weighted measure
      clr_x.w = log(x.w/exp(1/(sum(w))*sum(w*log(x.w))))
      # marginals
      yrow.w = apply(w*clr_x.w, 1, sum)/apply(w, 1, sum)
      ycol.w = apply(w*clr_x.w, 2, sum)/apply(w, 2, sum)
      #independent table
      yind.w = yrow.w*matrix(1, nrow = I, ncol = J) +
               t(ycol.w*t(matrix(1, nrow = I, ncol = J)))
      #interaction table
      yint.w = clr_x.w - yind.w
      # moving back to uniform reference measure
      clr_x.w.unif <- sqrt(w)*clr_x.w
      yind.w.unif <- sqrt(w)*yind.w
      yint.w.unif <- sqrt(w)*yint.w
      #singular value decomposition
      yint_dec = svd(yint.w.unif)
      u = yint_dec$u
      rownames(u) = rownames(x)
      v = yint_dec$v
      rownames(v) = colnames(x)
      # explained variances
      exp.variability.int.unweighted <- yint_dec$d^2/norm(as.matrix(yint),type="F")^2*100 # with respect to variability in the interaction table
      exp.variability.orig.unweighted <- yint_dec$d^2/norm(as.matrix(clr_x),type="F")^2*100 # with respect to variability of the complete table
      exp.variability.int.weighted <- yint_dec$d^2/sum(yint_dec$d^2)*100 # with respect to variability in the (weighted) interaction table
      exp.variability.orig.weighted <- yint_dec$d^2/norm(as.matrix(clr_x.w.unif),type="F")^2*100 # with respect to variability of the (weighted) complete table
      #simplicial deviance
      simp.deviance.weighted <- norm(as.matrix(yint.w.unif),type="F")^2/norm(as.matrix(clr_x.w.unif),type="F")^2
      return(list(u = u, v = v, 
                  intu = exp.variability.int.unweighted, 
                  oriu = exp.variability.orig.unweighted,
                  intw = exp.variability.int.weighted, 
                  oriw = exp.variability.orig.weighted,
                  simpw = simp.deviance.weighted))
  } 
}

# plot results
plotca <- function(res,weighted=FALSE){
  rgx <- range(c(res$u[,1],res$v[,1]))
  rgy <- range(c(res$u[,2],res$v[,2]))
  if (weighted){
    tit <- "Weighted CACT"
    plot(res$u[,1:2], xlim=rgx*1.15, ylim=rgy*1.15,
         xlab=paste0("Component 1 (", round(res$intw[1], 1), "/", round(res$oriw[1], 1), "%)"),
         ylab=paste0("Component 2 (", round(res$intw[2], 1), "/", round(res$oriw[2], 1), "%)"),
         pch="",cex.lab=1.2)
    title(tit)
  }
  else {
    tit <- "Unweighted CACT"
    plot(res$u[,1:2], xlim=rgx*1.15, ylim=rgy*1.15,
         xlab=paste0("Component 1 (", round(res$intu[1], 1), "/", round(res$oriu[1], 1), "%)"),
         ylab=paste0("Component 2 (", round(res$intu[2], 1), "/", round(res$oriu[2], 1), "%)"),
         pch="",cex.lab=1.2)
    title(tit)
  }
  text(res$u[,1:2], rownames(res$u),col="blue")
  points(res$v[,1:2], pch="")
  text(res$v[,1:2], rownames(res$v),col="red")
  abline(h=0,lty=3)
  abline(v=0,lty=3)
}

# more efficient implementation of CoDa CA:
mycaw <- function(x,weighted=TRUE){
  x <- as.matrix(x)
  L <- as.matrix(log(x))
  r <- apply(x,1,sum)/sum(x)
  c <- apply(x,2,sum)/sum(x)
  Z <- (diag(length(r))-rep(1,length(r))%*%t(r))%*%L%*%(diag(length(c))-c%*%t(rep(1,length(c))))
  r <- apply(x,1,sum)/sum(x)
  c <- apply(x,2,sum)/sum(x)
  if (weighted){
    Z <- diag(sqrt(r))%*%Z%*%diag(sqrt(c))
  }
  S.svd <- svd(Z)
  # Std coord of rows:
  phi <- diag(1/sqrt(r))%*%S.svd$u
  #phi <- S.svd$u
  # Std coord of cols:
  gam <- diag(1/sqrt(c))%*%S.svd$v
  # Princ. coord of rows
  Fr <- phi%*%diag(S.svd$d)
  rownames(Fr) <- rownames(x)
  # Princ. coord of cols:
  Gc <- gam%*%diag(S.svd$d)
  rownames(Gc) <- colnames(x)
  varcomp <- S.svd$d^2/sum(S.svd$d^2)
  list(Fr=Fr,Gc=Gc,var=varcomp)
}

# Bootstrap tables:
bootst <- function(x,B=100){
  xori <- x
  x[x==0] <- 2/3
  res <- coda_ca(x)
  w.ind <- (apply(x, 1, sum)/ncol(x)) %*% t(apply(x,2,sum)/nrow(x))
  resw <- coda_ca(x, w = w.ind, weighted = TRUE)
  rescaw <- mycaw(x)
  ang <- angw <- angcaw <- minnb <- simp <- simpw <- rep(NA,B)
  x <- xori
  xmod <- as.data.frame(as.table(data.matrix(x)))
  xlong <- xmod[rep(1:nrow(xmod), xmod$Freq), 1:2]
  for (i in 1:B){
    xboot <- xlong[sample(nrow(xlong),replace=TRUE),]
    xb <- table(xboot)  
    minnb[i] <- min(xb)
    #imputation
    xb[xb==0] <- 2/3
    res1 <- coda_ca(xb)
    w1.ind <- (apply(xb, 1, sum)/ncol(xb)) %*% t(apply(xb,2,sum)/nrow(xb))
    res1w <- coda_ca(xb, w = w1.ind, weighted = TRUE)
    res1caw <- mycaw(xb)
    ang[i] <- max(angle(res$u,res1$u[,1:2]),angle(res$v,res1$v[,1:2]))
    angw[i] <- max(angle(resw$u,res1w$u[,1:2]),angle(resw$v,res1w$v[,1:2]))
    angcaw[i] <- max(angle(resw$u,res1w$u[,1:2]),angle(resw$v,res1w$v[,1:2]))
    simp[i] <- res1$simp
    simpw[i] <- res1w$simpw
  }
  list(ang=ang,angw=angw,angcaw=angcaw,minnb=minnb,simp=simp,simpw=simpw)
}

