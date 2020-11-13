#' @title Second order conditional composite likelihood function
#' @description Second order conditional composite likelihood function.
#' @param psi Current value of alpha
#' @param xis Current value of xi (log transformed)
#' @param sigmas Current value of sigma^2 (log transformed)
#' @param phis Current value of phi (log transformed)
#' @param fi exp(beta^T z(u))
#' @param fj exp(beta^T z(v))
#' @param R distance
#' @param mi mark of u
#' @param mj mark of v
#' @param mval all marks
#' @param stz indicator
#' @return likelihood
#' @export

second_order_likelihood <- function(psi,xis,sigmas,phis,R,fi,fj,mi,mj,mval,stz=F){
  
   p = length(mval)
   if(stz==T){
    B=rbind(diag(p-1),rep(-1,p-1)) 
    alphas=B%*%psi
   }
   else{
     alphas=psi
   }
   q=dim(alphas)[2]

  if(!is.matrix(alphas)){alphas=matrix(alphas,nrow=p)}

  nR=length(R)
  
  Ptmp <- matrix(NA,nR,p^2) 
  ptmp <- numeric(nR)
  k0 = 0
  for(k1 in 1:p){
    for(k2 in 1:p){
      k0 = k0 +1
      ind = (mi==mval[k1])&(mj==mval[k2])
      tmp1 <- 0
      for(l in 1:q){tmp1 <- tmp1+alphas[k1,l]*alphas[k2,l]*exp(-R/exp(xis[l]))}
      if(k1==k2){Ptmp[,k0] <-  fi[,k1]*fj[,k2]*exp(tmp1+exp(sigmas[k1])*exp(-R/exp(phis[k1])))}
      else{Ptmp[,k0] <-  fi[,k1]*fj[,k2]*exp(tmp1)}
      ptmp[ind] = Ptmp[ind,k0]
    }
  }
  
  fval = ptmp/rowSums(Ptmp)
  fval = sum(log(fval))
  -fval
}
