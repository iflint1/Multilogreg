#' @title Sensitivity of phi
#' @description Sensitivity of phi
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
#' @return Sensitivity of phi
#' @export

sensitivity.phis=function(alphas=alphas,xis,sigmas,phis,R=0,fi,fj,mi=NULL,mj=NULL,mval=NULL){
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

  # phi phi

  h1 = matrix(0,p,p)
  for(m in 1:p){
    k0 = 0
    grad.g1 = 0
    for(i in 1:p){
      for(j in 1:p){
        k0 = k0 +1
        if(i ==m & j == m){
          tmp <- Ptmp[,k0]
          grad.g1 = grad.g1 + sum((exp(sigmas[m])*exp(-R/exp(phis[m]))*(R/exp(phis[m])^2))^2*Ptmp[,k0]/aggPtmp)
        }
      }
    }
    h1[m,m] = grad.g1
  }

  grad.phi = matrix(0,p,length(R))
  for(m in 1:p){
    k0 = 0
    for(i in 1:p){
      for(j in 1:p){
        k0 = k0 +1
        if(m == i & m == j){
          tmp <- Ptmp[,k0]
          grad.phi[m,] = grad.phi[m,] + exp(sigmas[m])*exp(-R/exp(phis[m]))*(R/exp(phis[m])^2)*Ptmp[,k0]/aggPtmp
        }
      }
    }
  }

  H.phi = h1 - grad.phi%*%t(grad.phi)
  H.phi

}

