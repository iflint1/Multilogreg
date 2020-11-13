#' @title Sensitivity of alpha
#' @description Sensitivity of alpha
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
#' @return Sensitivity of alpha
#' @export

sensitivity.alpha1=function(alphas=NULL,xis=NULL,sigmas=NULL,phis=NULL,R=0,fi=NULL,fj=NULL,mi=NULL,mj=NULL,mval=NULL){
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

  h1 = matrix(0,p*q,p*q)
  m0 = 0
  for(m1 in 1:p){
    for(m2 in 1:p){
      m0 = m0 +1
      for(l1 in 1:q){
        for(l2 in 1:q){
          grad.g1 <- numeric(1)
          if(m1 == m2){
            k0 = 0
            for(i in 1:p){
              for(j in 1:p){
                k0 = k0 + 1
                tmp <- Ptmp[,k0]
                if(i==m1 & j==m1){
                  tmp2 <- 2*alphas[i,l1]*exp(-R/exp(xis[l1]))
                  tmp3 <- 2*alphas[j,l2]*exp(-R/exp(xis[l2]))
                  grad.g1 <- grad.g1 + sum(tmp2*tmp3*tmp/aggPtmp)
                }
                if(i==m1 & j!=m1){
                  tmp2 <- alphas[j,l1]*exp(-R/exp(xis[l1]))
                  tmp3 <- alphas[j,l2]*exp(-R/exp(xis[l2]))
                  grad.g1 <- grad.g1 + sum(tmp2*tmp3*tmp/aggPtmp)
                }
                if(i!=m1 & j==m1){
                  tmp2 <- alphas[i,l1]*exp(-R/exp(xis[l1]))
                  tmp3 <- alphas[i,l2]*exp(-R/exp(xis[l2]))
                  grad.g1 <- grad.g1 + sum(tmp2*tmp3*tmp/aggPtmp)
                }
              }
            }
          }
          if(m1 != m2){
            t1=m2+p*(m1-1)
            t2=m1+p*(m2-1)
            ptmp1 <- Ptmp[,t1]
            ptmp2 <- Ptmp[,t2]
            tmp2 <- alphas[m1,l2]*exp(-R/exp(xis[l2]))
            tmp3 <- alphas[m2,l1]*exp(-R/exp(xis[l1]))
            grad.g1 <- grad.g1 + sum(tmp2*tmp3*ptmp1/aggPtmp)
            grad.g1 <- grad.g1 + sum(tmp2*tmp3*ptmp2/aggPtmp)
          }
          t1=m1+p*(l1-1)
          t2=m2+p*(l2-1)
          h1[t1,t2] = grad.g1
        }
      }
    }
  }

  grad.g2 = matrix(0,p*q,length(R))
  for(m in 1:p){
    for(l in 1:q){
      t=m+p*(l-1)
      k0 = 0
      for(k1 in 1:p){
        for(k2 in 1:p){
          k0 = k0 + 1
          tmp <- Ptmp[,k0]*exp(-R/exp(xis[l]))
          if(k1==m & k2==m){grad.g2[t,] <- grad.g2[t,] + tmp*2*alphas[k1,l]/aggPtmp}
          if(k1==m & k2!=m){grad.g2[t,] <- grad.g2[t,] + tmp*alphas[k2,l]/aggPtmp}
          if(k2==m & k1!=m){grad.g2[t,] <- grad.g2[t,] + tmp*alphas[k1,l]/aggPtmp}
        }
      }
    }
  }
  H = h1 - grad.g2%*%t(grad.g2)
  H
}
