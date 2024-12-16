########################
### Correspondence analysis of a compositional table
########################

# we have to define
# matrix x - original IxJ CoDa table
# matrix w - IxJ matrix of weights. If not defined, unweighted CA is computed
# to.weight - "orig.x" or "int.x" or "orig.x.dir" - defines if the original table or directly its
#             interactive part is weighted, eventually weighted is the original table without further 
#             decomposition

#for weighted tables, the explained variability is computed wrt. WEIGHTED objects


coda_ca <- function(x, w = NULL, to.weight = "orig.x")
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
  
  if(is.null(w))   # unweighted CA
  {
    #singular value decomposition
    yint_dec = svd(yint)
    
    u = yint_dec$u
    rownames(u) = rownames(x)
    v = yint_dec$v
    rownames(v) = colnames(x)
    
    # explained variances
    exp.variability.int.unweighted <- yint_dec$d^2/sum(yint_dec$d^2)*100 # with respect to variability in the interaction table
    exp.variability.orig.unweighted <- yint_dec$d^2/norm(as.matrix(clr_x),type="F")^2*100 # with respect to variability of the complete table
 
    #simplicial deviance
    simp.deviance <- norm(as.matrix(yint),type="F")^2/norm(as.matrix(clr_x),type="F")^2
       
    return(list(u = u, v = v, 
                intu = exp.variability.int.unweighted, 
                oriu = exp.variability.orig.unweighted,
                intw = NULL, 
                oriw = NULL,
                simp = simp.deviance))
  } else
    if(to.weight == "orig.x")
    {
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
  else
    if(to.weight == "orig.x.dir")
  {
    # weighting is applied on the complete table (without further decomposition)
    x.w <- x/w
    
    # clr representation wrt. weighted measure
    clr_x.w = log(x.w/exp(1/(sum(w))*sum(w*log(x.w))))
    
    # moving back to uniform reference measure
    clr_x.w.unif <- sqrt(w)*clr_x.w
    
    #singular value decomposition
    y_dec = svd(clr_x.w.unif)
    
    u = y_dec$u
    rownames(u) = rownames(x)
    v = y_dec$v
    rownames(v) = colnames(x)
    
    # explained variances
    exp.variability.int.unweighted <- NULL # with respect to variability in the interaction table
    exp.variability.orig.unweighted <- y_dec$d^2/norm(as.matrix(clr_x),type="F")^2*100 # with respect to variability of the complete table
    exp.variability.int.weighted <- NULL # with respect to variability in the (weighted) interaction table
    exp.variability.orig.weighted <- y_dec$d^2/norm(as.matrix(clr_x.w.unif),type="F")^2*100 # with respect to variability of the (weighted) complete table
    
    return(list(u = u, v = v, 
                intu = exp.variability.int.unweighted, 
                oriu = exp.variability.orig.unweighted,
                intw = exp.variability.int.weighted, 
                oriw = exp.variability.orig.weighted))
    
  }
  else
    {
      # weighting is applied on the interaction table
      xint.w <- exp(yint) / w
      
      # weighted interaction table in clr
      yint.w = log(xint.w/exp(1/(sum(w))*sum(w*log(xint.w))))
      
      # moving back to uniform reference measure
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
      exp.variability.orig.weighted <- NULL # with respect to variability of the (weighted) complete table
      
      return(list(u = u, v = v, 
                  intu = exp.variability.int.unweighted, 
                  oriu = exp.variability.orig.unweighted,
                  intw = exp.variability.int.weighted, 
                  oriw = exp.variability.orig.weighted))
      
    }
}


# plot result
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

