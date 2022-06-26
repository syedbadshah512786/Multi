#' Canonical Correlation of Covariance or Correlation Matrix
#' @param sigma must ba a covariance or correlation matrix
#' @description The Discrimination & Classification are multivariate techniques concered with separating distint sets of objects (or observations) and with allocating new objects (observations) to perviously define groups.
#' @references Inspired by Gold Medalist Dr Shahla Faisal who study me in GCUF
#' @author Syed Hammad Raza (Studies BS Statistics in GCUF)
#' @export
ccs <- function(sigma){
  cat("The sigma11 is: \n")
  sigma11<-sigma[1:2,1:2]
  print(sigma11)

  cat("The sigma12 is: \n")
  sigma12<-sigma[1:2,3:4]
  print(sigma12)

  cat("The sigma21 is: \n")
  sigma21<-sigma[3:4,1:2]
  print(sigma21)

  cat("The sigma22 is: \n")
  sigma22<-sigma[3:4,3:4]
  print(sigma22)

  cat("The inverse of sigma11 is: \n")
  sigma11.inv <- solve(sigma11)
  print(sigma11.inv)

  cat("The inverse of sigma22 is: \n")
  sigma22.inv <- solve(sigma22)
  print(sigma22.inv)

  E11 <- eigen(sigma11)
  sigma11.sqrt.inv <- (1/sqrt(E11$value[1])) * E11$vectors[,1] %*% t(E11$vectors[,1])+
    (1/sqrt(E11$value[2])) * E11$vectors[,2] %*% t(E11$vectors[,2])
  cat("The square root inverse of sigma11 is: \n")
  print(sigma11.sqrt.inv)

  E22 <- eigen(sigma22)
  sigma22.sqrt.inv <- (1/sqrt(E22$value[1])) * E22$vectors[,1] %*% t(E22$vectors[,1])+
    (1/sqrt(E22$value[2])) * E22$vectors[,2] %*% t(E22$vectors[,2])
  cat("The square root inverse of sigma22 is: \n")
  print(sigma22.sqrt.inv)


  sigma.group.1 <- sigma11.sqrt.inv %*% sigma12 %*% sigma22.inv %*% sigma21 %*% sigma11.sqrt.inv
  cat("The sigma of the Group 1 is: \n")
  print(sigma.group.1)

  sigma.group.2 <- sigma22.sqrt.inv %*% sigma21 %*% sigma11.inv %*% sigma12 %*% sigma22.sqrt.inv
  cat("The sigma of the Group 2 is: \n")
  print(sigma.group.2)

  e.group.1 <- eigen(sigma.group.1)
  e.group.2 <- eigen(sigma.group.2)

  cat("The eigen values & Vectors of Group 1 is: \n")
  print(e.group.1)
  cat("The eigen values & Vectors of Group 2 is: \n")
  print(e.group.2)

  #-------Canonical Correlation are

  p1 <- (t(e.group.1$vectors[,1]) %*% sigma12 %*% e.group.2$vectors[,1]) /
    (sqrt((t(e.group.1$vectors[,1]) %*% sigma11 %*% e.group.1$vectors[,1])) %*%
       sqrt((t(e.group.2$vectors[,1]) %*% sigma22 %*% e.group.2$vectors[,1])))
  cat("The correlation between U1 & V1 is: \n")
  print(p1)

  p2 <- (t(e.group.1$vectors[,2]) %*% sigma12 %*% e.group.2$vectors[,2]) /
    (sqrt((t(e.group.1$vectors[,2]) %*% sigma11 %*% e.group.1$vectors[,2])) %*%
       sqrt((t(e.group.2$vectors[,2]) %*% sigma22 %*% e.group.2$vectors[,2])))
  cat("The correlation between U2 & V2 is: \n")
  print(p2)
}



#' Canonical Correlation Analysis
#' @param X is a 1st Group
#' @param Y is a 2nd Group
#' @description The Discrimination & Classification are multivariate techniques concered with separating distint sets of objects (or observations) and with allocating new objects (observations) to perviously define groups.
#' @references Inspired by Gold Medalist Dr Shahla Faisal who study me in GCUF
#' @author Syed Hammad Raza (Studies BS Statistics in GCUF)
#' @export

cc <- function (X, Y)
{
  Xnames = dimnames(X)[[2]]
  Ynames = dimnames(Y)[[2]]
  ind.names = dimnames(X)[[1]]
  res = rcc(X, Y, 0, 0)
  return(res)
}



comput<-function (X, Y, res)
{
  X.aux = scale(X, center = TRUE, scale = FALSE)
  Y.aux = scale(Y, center = TRUE, scale = FALSE)
  X.aux[is.na(X.aux)] = 0
  Y.aux[is.na(Y.aux)] = 0
  xscores = X.aux %*% res$xcoef
  yscores = Y.aux %*% res$ycoef
  corr.X.xscores = cor(X, xscores, use = "pairwise")
  corr.Y.xscores = cor(Y, xscores, use = "pairwise")
  corr.X.yscores = cor(X, yscores, use = "pairwise")
  corr.Y.yscores = cor(Y, yscores, use = "pairwise")
  return(list(xscores = xscores, yscores = yscores, corr.X.xscores = corr.X.xscores,
              corr.Y.xscores = corr.Y.xscores, corr.X.yscores = corr.X.yscores,
              corr.Y.yscores = corr.Y.yscores))
}

rcc <- function (X, Y, lambda1, lambda2)
{
  Xnames <- dimnames(X)[[2]]
  Ynames <- dimnames(Y)[[2]]
  ind.names <- dimnames(X)[[1]]
  Cxx <- var(X, na.rm = TRUE, use = "pairwise") + diag(lambda1,
                                                       ncol(X))
  Cyy <- var(Y, na.rm = TRUE, use = "pairwise") + diag(lambda2,
                                                       ncol(Y))
  Cxy <- cov(X, Y, use = "pairwise")
  res <- geigen(Cxy, Cxx, Cyy)
  names(res) <- c("cor", "xcoef", "ycoef")
  scores <- comput(X, Y, res)
  return(list(cor = res$cor, names = list(Xnames = Xnames,
                                          Ynames = Ynames, ind.names = ind.names), xcoef = res$xcoef,
              ycoef = res$ycoef, scores = scores))
}

geigen<-function (Amat, Bmat, Cmat)
{
  Bdim <- dim(Bmat)
  Cdim <- dim(Cmat)
  if (Bdim[1] != Bdim[2])
    stop("BMAT is not square")
  if (Cdim[1] != Cdim[2])
    stop("CMAT is not square")
  p <- Bdim[1]
  q <- Cdim[1]
  s <- min(c(p, q))
  if (max(abs(Bmat - t(Bmat)))/max(abs(Bmat)) > 1e-10)
    stop("BMAT not symmetric.")
  if (max(abs(Cmat - t(Cmat)))/max(abs(Cmat)) > 1e-10)
    stop("CMAT not symmetric.")
  Bmat <- (Bmat + t(Bmat))/2
  Cmat <- (Cmat + t(Cmat))/2
  Bfac <- chol(Bmat)
  Cfac <- chol(Cmat)
  Bfacinv <- solve(Bfac)
  Cfacinv <- solve(Cfac)
  Dmat <- t(Bfacinv) %*% Amat %*% Cfacinv
  if (p >= q) {
    result <- svd(Dmat)
    values <- result$d
    Lmat <- Bfacinv %*% result$u
    Mmat <- Cfacinv %*% result$v
  }
  else {
    result <- svd(t(Dmat))
    values <- result$d
    Lmat <- Bfacinv %*% result$v
    Mmat <- Cfacinv %*% result$u
  }
  geigenlist <- list(values, Lmat, Mmat)
  names(geigenlist) <- c("values", "Lmat", "Mmat")
  return(geigenlist)
}





#' Leave-one-out (Lachenbruch "Holdout") Procedure for equal and unequal variance
#'
#' @param x is a 1st Group
#' @param y is a 2nd group
#' @param i is such a variable which you went to test whether it from 1st or 2nd Group
#' @param var is use to select a method whether it has for equal or unequal variance (Var="Equal" for equal and Var="NotEqual" for unequal variance)
#' @param g is use to define a group i.e the test value is belong from 1 group or from 2nd group (g=1 for 1st Group and g=2 for 2nd group)
#' @references Inspired by Gold Medalist Dr Shahla Faisal who study me in GCUF
#' @author Syed Hammad Raza (Studies BS Statistics in GCUF)
#' @export

s.holdout <-   function(x,y,i,var,g)
{
  if(var=="Equal")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Equal \n")
      xo <- x[i,]
      x1 <- x[-i,]
      x2 <- y

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Equal \n")
      xo <- y[i,]
      x1 <- x
      x2 <- y[-i,]

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 2nd Group is Misclassified")} else {cat(" The",i,"Observation of 2nd Group is Correctly classified") }
    }
  }
  if(var=="NotEqual")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Not Equal \n")
      x1 <- x[-i,]
      x2 <- y
      xo  <- x[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS >= RHS){cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Not Equal \n")
      x1  <- x
      x2  <- y[-i,]
      xo  <- y[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS <= RHS){cat( "The",i,"Observation of 2nd Group is Correctly classified")} else {cat(" The",i,"Observation of 2nd Group is Misclassified") }
    }
  }
}

#' Classification of Multivariate Data Procedure for equal and unequal variance
#'
#' @param x is a 1st Group
#' @param y is a 2nd group
#' @param i is such a variable which you went to test whether it from 1st or 2nd Group
#' @param var is use to select a method whether it has for equal or unequal variance (Var="Equal" for equal and Var="NotEqual" for unequal variance)
#' @param g is use to define a group i.e the test value is belong from 1 group or from 2nd group (g=1 for 1st Group and g=2 for 2nd group)
#' @references Inspired by Gold Medalist Dr Shahla Faisal who study me in GCUF
#' @author Syed Hammad Raza(Studies BS Statistics in GCUF)
#' @export

CF.method <-   function(x,y,i,var,g)
{
  if(var=="Equal")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Equal \n")
      xo <- x[i,]
      x1 <- x
      x2 <- y

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Equal \n")
      xo <- y[i,]
      x1 <- x
      x2 <- y

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 2nd Group is Misclassified")} else {cat(" The",i,"Observation of 2nd Group is Correctly classified") }
    }
  }
  if(var=="NotEqual")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Not Equal \n")
      x1 <- x
      x2 <- y
      xo  <- x[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS >= RHS){cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Not Equal \n")
      x1  <- x
      x2  <- y
      xo  <- y[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS <= RHS){cat( "The",i,"Observation of 2nd Group is Correctly classified")} else {cat(" The",i,"Observation of 2nd Group is Misclassified") }
    }
  }
}


