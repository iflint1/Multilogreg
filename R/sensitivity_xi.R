#' @title Sensitivity of xi
#' @description Sensitivity of xi
#' @param alphas Current value of alpha
#' @param xis Current value of xi (log transformed)
#' @param sigmas Current value of sigma^2 (log transformed)
#' @param phis Current value of phi (log transformed)
#' @param fi exp(beta^T z(u))
#' @param fj exp(beta^T z(v))
#' @param R distance
#' @param mi mark of u
#' @param mj mark of v
#' @param mval all marks
#' @return Sensitivity of xi
#' @export

sensitivity.xi=function(alphas,xis,sigmas,phis,R=0,fi=matrix(1,ncol = p),fj=matrix(1,ncol = p),mi=NULL,mj=NULL,mval=NULL){
  p=length(mval)
  if(!is.matrix(alphas)){alphas=matrix(alphas,nrow=p)}
  q = dim(alphas)[2]
  k0=0
  Ptmp <- matrix(NA,length(R),p^2)
  for(k1 in 1:p){
    for(k2 in 1:p){
      k0 = k0 + 1
      tmp1 <- 0
      for(l in 1:q){tmp1 <- tmp1+alphas[k1,l]*alphas[k2,l]*exp(-R/exp(xis[l]))}
      if(k1==k2){Ptmp[,k0] <-  fi[,k1]*fj[,k2]*exp(tmp1+exp(sigmas[k2])*exp(-R/exp(phis[k2])))}
      if(k1!=k2){Ptmp[,k0] <-  fi[,k1]*fj[,k2]*exp(tmp1)}
    }
  }
  aggPtmp <- rowSums(Ptmp)

  # Xi Xi
  h1 = matrix(0,q,q)
  for(l1 in 1:q){
    for(l2 in 1:q){
      grad.g1 <- numeric(1)

      k0 = 0
      for(i in 1:p){
        for(j in 1:p){
          k0 = k0 + 1
          tmp <- Ptmp[,k0]
          tmp2 <- alphas[i,l1]*alphas[j,l1]*exp(-R/exp(xis[l1]))*(R/(exp(xis[l1])^2))
          tmp3 <- alphas[i,l2]*alphas[j,l2]*exp(-R/exp(xis[l2]))*(R/(exp(xis[l2])^2))
          grad.g1 <- grad.g1 + sum(tmp2*tmp3*tmp/aggPtmp)
        }
      }
      h1[l1,l2] = grad.g1
    }
  }

  grad.xi = matrix(0,q,length(R))
  for(l in 1:q){
    k0 = 0
    for(k1 in 1:p){
      for(k2 in 1:p){
        k0 = k0 + 1
        grad.xi[l,] = grad.xi[l,] + alphas[k1,l]*alphas[k2,l]*exp(-R/exp(xis[l]))*(R/(exp(xis[l])^2))*Ptmp[,k0]/aggPtmp
      }
    }
  }
  H.xi = h1 - grad.xi%*%t(grad.xi)

  H.xi
}
