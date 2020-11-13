#' @title Simulation of semi-parametric multivariate log Gaussian Cox processes
#'
#' @description Simulation of semi-parametric multivariate log Gaussian Cox processes.
#'
#' @param basecov Background intensity rho_0.
#' @param covariate Optional. A simulated covariate. The covariate must be a matrix.
#' @param betas A matrix with covariates.
#' @param alphas Alpha parameters. Must be a matrix, where the number of rows correspond
#' to the number of components in the LGCP. The number of columns correspond to the number
#' of common latent field.
#' @param xis Correlation scale parameters for each common random field. The correlation functions
#' for the common latent fields are exponential.
#' @param sigmas Sigma parameters. The number of sigma parameters must correspond
#' to the number of components in the LGCP.
#' @param phis Correlation scale parameters for each type-specific random field. The correlation functions
#' for the type-specific random fields are exponential.
#' @param n.window window size.
#' @param n.points Expected number of point for each component in the LGCP.
#' The length of n.points must correspond to the number of components.
#' @param beta0s Intercepts. The length of beta0s must correspond to the number of components
#' in the LGCP.
#' @importFrom RandomFields RFsimulate
#'
#' @return Multivariate LGCP
#' @author Kristian Bjørn Hessellund, Ganggang Xu, Yongtao Guan and Rasmus Waagepetersen.
#' @references Hessellund, K. B., Xu, G., Guan, Y. and Waagepetersen, R. (2020)
#' Second order semi-parametric inference for multivariate log Gaussian Cox processes.
#'
#' @examples
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
#' plot(X$markedprocess)
#'

#' @export
#'
sim_lgcp_multi <- function(basecov,covariate,betas,alphas,xis,sigmas,phis, n.window,n.points,beta0s=NULL){

  nspecies <- dim(betas)[1]
  dimyx <- ncol(covariate)
  n.latent <- dim(alphas)[2]
  ###find the intercepts such that there are n.points[i] for each species
  logintensity.case <- matrix(0,nspecies,dimyx^2)

  ###If all intercepts are given, then use given beta0s; otherwise calculate beta0s so that
  ### average number of points are n.points[i] for process X_i

  if(is.null(beta0s)){
    beta0s <- numeric()
    for(i in 1:nspecies){
      beta0s[i] <- log(n.points[i])-mean(c(as.matrix(covariate*betas[i,1]+basecov)))
      logintensity.case[i,] <- c(t(as.matrix(beta0s[i]+covariate*betas[i,1]+basecov)))
    }
  }else{
    for(i in 1:nspecies){
      logintensity.case[i,] <- c(t(as.matrix(beta0s[i]+covariate*betas[i,1]+basecov+log(n.points[i]))))
    }
  }

  ###Create sqaure window for point processes
  gridmask <- as.mask(square(n.window), dimyx = dimyx)
  gridcord <- rasterxy.mask(gridmask)
  ###The uniq x and y seq
  xx <- unique(gridcord$x)
  yy <- unique(gridcord$y)

  ###Generate a Gaussian Random Field for the common latent intensity
  for (j in 1:n.latent){
    #Y1 <- RFsimulate(RMexp(var=1,scale=scales_c[j]), x=xx, y=yy, grid=TRUE)
    Y1 <- RFsimulate(RMexp(var=1,scale=xis[j]), x=xx, y=yy, grid=TRUE)
    for (l in 1:nspecies)
      logintensity.case[l,]=logintensity.case[l,]+alphas[l,j]*c(t(as.matrix(Y1)))
  }

  ###Generate a Gaussian Random Field for the individual latent intensity
  for (l in 1:nspecies){
    #U = RFsimulate(RMexp(var=1,scale=scales_res[l]), x=xx, y=yy, grid=TRUE)
    U = RFsimulate(RMexp(var=1,scale=phis[l]), x=xx, y=yy, grid=TRUE)
    logintensity.case[l,]=logintensity.case[l,]+sigmas[l]*c(t(as.matrix(U)))
  }

  mu <- process <- list()
  markedprocess <- NULL
  for (l in 1:nspecies){
    ###The intensity of LGCP
    intensity_lgcp <- exp(logintensity.case[l,]-sum(alphas[l,]^2)/2-sigmas[l]^2/2)
    intensity.df <- data.frame("x"= gridcord$x, "y" = gridcord$y, "value" = c(intensity_lgcp))
    mu[[l]] <- as.im(intensity.df)
    ###Simulate inhomogeneous poisson process from lgcp intensity
    process[[l]] <- rpoispp(lambda = mu[[l]])
  }


  ##Exact the pixel image for covariates
  pix <- as.im(data.frame("x"= gridcord$x, "y" = gridcord$y, "value" = c(t(covariate))))


  ###Aggregate all points into one marked point process
  win <- process[[1]]$window
  tmp <- NULL
  for (l in 1:nspecies){
    tmp <-rbind(tmp,cbind(process[[l]]$x,process[[l]]$y,l))
  }

  markedprocess <- ppp(x=tmp[,1],y=tmp[,2],marks=as.factor(tmp[,3]),window=win)


  list(process = process, mu = mu, pix = pix, beta0s=beta0s, markedprocess=markedprocess)
}

