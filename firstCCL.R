#' @title Semi-parametric multinomial logistic regression for multivariate point pattern data
#'
#' @description Semi-parametric multinomial logistic regression for multivariate point pattern data.
#'
#' @param X A multivariate point pattern.
#' @param Beta0 Initial value for Beta.
#' @param covariate A matrix or image of covariates.
#' @importFrom stats optim
#' @importFrom stats runif
#' @return Return estimated first order parameters of X
#' @author Kristian Bjørn Hessellund, Ganggang Xu, Yongtao Guan and Rasmus Waagepetersen.
#' @references Hessellund, K. B., Xu, G., Guan, Y. and Waagepetersen, R. (2020)
#' Semi-parametric multinomial logistic regression for multivariate point pattern data.
#' @examples
#' ### First order analysis of Washington D.C. street crimes ###
#'
#' nspecies <- length(unique(markedprocess$marks))
#' ncovs <- dim(Z_process)[2]
#'
#' # Initial parameters for street crime data
#' Beta0=matrix(runif(nspecies*ncovs,-0.5,0.5),nrow =nspecies)
#'
#' # Choosing the last type of street crime as the baseline
#' Beta0 <- as.matrix(t(scale(Beta0,center=Beta0[nspecies,],scale=FALSE)))
#'
#' # Parameter estimation
#' Betahat=FirstOrderCCL(X=markedprocess,Beta0 = Beta0,covariate = Z_process)
#'
#' # Parameter estimates
#' colnames(Betahat$betahat)=c("Burglary", "Assault w. weapon", "Motor v. theft",
#' "Theft f. auto", "Robbery", "Other theft (baseline)")
#' rownames(Betahat$betahat)=c("Intercept", "%African", "%Hispanic", "%Males",
#'  "Median income", "%Household","%Bachelor", "Distance to police station")
#' Betahat
#'
#' ###################################################################
#'
#' ### First order analysis of simulated LGCP ###
#'
#' # Size of the observation window
#' n.x <- n.y <- 1
#' xx=seq(0,n.x,length=100)
#' yy=seq(0,n.x,length=100)
#'
#' # Simulating a covariate
#' cov <- as.matrix(RFsimulate(RMexp(var=1,scale=0.05), x=xx, y=yy, grid=TRUE))
#'
#' # Simulating the background intensity
#' gamma <- 0.5
#' background <- as.matrix(RFsimulate(RMgauss(var=1,scale=0.2), x=xx, y=yy, grid=TRUE))*gamma-gamma^2/2
#'
#' #Set up parameters
#' beta1 <- c(0.1,0.2,0.3,0.4,0.5)
#' beta2 <- c(-0.1,-0.2,0,0.1,0.2)
#' beta2 <- as.matrix(beta2)
#'
#' # Parameters in the LGCP
#'
#' alpha <- matrix(c(0.5,-1,0.5,0,-1,0,0,0.5,0,0.5),nrow=5,byrow=TRUE)
#' xi    <- c(0.02,0.03)
#' sigma <- matrix(c(sqrt(0.5),sqrt(0.5),sqrt(0.5),sqrt(0.5),sqrt(0.5)),ncol=1)
#' phi   <- matrix(c(0.02,0.02,0.03,0.03,0.04),ncol=1)
#'
#' n.window <- n.x
#' n.points <- c(400,400,400,400,400)
#'
#' # Simulation of a multivariate LGCP
#' X <- sim_lgcp_multi(basecov=background,covariate=cov,betas=beta2,alphas=alpha,xi=xi,
#' sigma=sigma,phis=phi, n.window=n.window,n.points=n.points,beta0s=beta1)
#'
#' nspecies  <- dim(beta2)[1]
#' Beta.true <- cbind(beta1,beta2)
#' Beta.true <- as.matrix(t(scale(Beta.true,center=Beta.true[nspecies,],scale=FALSE)))
#' Betahat   <- FirstOrderCCL(X=X$markedprocess,Beta0=Beta.true,covariate = as.im(t(cov)))
#' Betahat
#'
#' @export
#'

FirstOrderCCL <- function(X, Beta0,covariate){

  n       = X$n
  mark.pp = sort(unique(X$marks))
  p       = length(mark.pp)
  process = list()

  for(i in 1:p){
    process[[i]]=X[mark.pp[i]==X$marks]
  }

  Z=list()
  if(is.null(covariate)==T){
    for(i in 1:p){
      x <- process[[i]]
      Z[[i]] <-   cbind(rep(1,x$n))
    }
    q=2
  }
  if(is.null(covariate)==F){
    if(is.im(covariate)==T){
      for(i in 1:p){
        Z[[i]] <-   cbind(1,covariate[process[[i]]])
      }
      q=2

    }
    if(is.matrix(covariate)==T){
      for(i in 1:p){
        tmp=mark.pp[i]
        ind=(rep(tmp,n)==X$marks)
        Z[[i]]=covariate[ind,]
      }
      q=dim(covariate)[2]
    }
  }
  #return(Z)
  ##Set up the initial values for betas
  par0 <- c(Beta0[,-p])

  ###Define the first order conditional likelihood
  loglikfun <- function(par){
    Beta <- cbind(matrix(par,nrow=q),rep(0,q))

    loglik <- 0
    for(i in 1:p){
      Ztmp <- Z[[i]]
      ptmp <- exp(Ztmp%*%Beta)
      ptmp <- ptmp[,i]/rowSums(ptmp)
      loglik <- loglik + sum(log(ptmp))
    }
    loglik
  }

  obj <- optim(par0,loglikfun,method = "BFGS",control=list(fnscale=-1))
  betahat <- matrix(obj$par,nrow = q)
  return(list(betahat=cbind(betahat,rep(0,q)),converg=obj$convergence))
}


