#' Small Area Estimation using Hierarchical Bayesian under Zero Inflated Negative Binomial Distribution
#' @description This function is implemented to variable of interest \eqn{(y)} that assumed to be a Zero Inflated Negative Binomial Distribution. The range of data is \eqn{(y >= 0)}. This model can be used to handle overdispersion and excess zero in data.
#' @param formula Formula that describe the fitted model
#' @param iter.update Number of updates with default \code{3}
#' @param iter.mcmc Number of total iterations per chain with default \code{1100}
#' @param coef.nonzero Optional vector containing initial values \code{mu.b} for the mean of the prior distribution of the log model coefficients \eqn{(\beta)} with default \code{rep(0,nvar)}
#' @param coef.zero Optional vector containing initial values \code{mu.g} for the mean of the prior distribution of the logit model coefficients \eqn{(\gamma)} with default \code{rep(0,nvar)}
#' @param var.coef.nonzero Optional vector containing initial values \code{tau.b} for the variance of the prior distribution on the log model coefficients \eqn{(\beta)} with default \code{rep(1,nvar)}
#' @param var.coef.zero Optional vector containing initial values \code{tau.g} for the variance of the prior distribution of the logit model coefficients \eqn{(\gamma)} with default \code{rep(1,nvar)}
#' @param thin Thinning rate, must be a positive integer with default \code{1}
#' @param burn.in  Number of iterations to discard at the beginning with default \code{600}
#' @param tau.u Variance of random effect area for non-zero count of variable interest with default \code{1}
#' @param tau.v Variance of random effect area for zero count of variable interest with default \code{1}
#' @param data The data frame
#'
#' @return This function returns a list of the following objects:
#'    \item{Est}{A dataframe that contains the values, standar deviation, and quantile of Small Area mean Estimates using Hierarchical bayesian method}
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A data frame with the estimated model coefficient consist of \code{beta} (coefficient in the log model) and \code{gamma} (coefficient in the logit model)}
#'    \item{plot.beta}{Trace, Density, Autocorrelation Function Plot of MCMC samples \code{beta}}
#'    \item{plot.gamma}{Trace, Density, Autocorrelation Function Plot of MCMC samples \code{gamma}}
#' @export ZinbHB
#' @import coda
#' @import rjags
#' @import stats
#' @import stringr
#' @import graphics
#' @import grDevices
#' @examples
#' ## Compute Fitted Model
#' ## y ~ x1 +x2, nvar = 3
#'
#' ## For data without any nonsampled area
#' ## Load Dataset
#' data(dataZINB)
#' result <- ZinbHB(formula = y ~ x1 + x2, data = dataZINB)
#'
#' ## Result
#' result$Est                             # Small Area mean Estimates
#' result$refVar                          # refVar
#' result$coefficient                     # coefficient
#'
#'
#' # Load library 'coda' to execute the plot
#' # autocorr.plot(result$plot.beta[[3]])  # Generate ACF Plot beta
#' # plot(result$plot.beta[[3]])           # Generate Dencity and Trace plot beta
#' # autocorr.plot(result$plot.gamma[[3]]) # Generate ACF Plot gamma
#' # plot(result$plot.gamma[[3]])          # Generate Dencity and trace plot gamma
#'
#'
#' ## For data with nonsampled area use dataZINBNS
#'
#'
ZinbHB <- function(formula, iter.update = 3, iter.mcmc = 1100, coef.nonzero, coef.zero,
                   var.coef.nonzero, var.coef.zero, thin = 1, burn.in = 600,
                   tau.u = 1, tau.v = 1, data){


  result <- list(Est = NA, refVar = NA, coefficient = NA,
                 plot.beta = NA, plot.gamma = NA)


  formuladata <- model.frame(formula,data,na.action=NULL)
  if (any(is.na(formuladata[,-1])))
    stop("Auxiliary Variables contains NA values.")
  auxVar <- as.matrix(formuladata[,-1])
  nvar <- ncol(auxVar) + 1

  if (!missing(var.coef.nonzero)){

    if( length(var.coef.nonzero) != nvar ){
      stop("length of vector var.coef.nonzero does not match the number of regression coefficients, the length must be ",nvar)
    }

    tau.b.value = 1/var.coef.nonzero
  } else {
    tau.b.value = 1/rep(1,nvar)
  }

  if (!missing(var.coef.zero)){

    if( length(var.coef.zero) != nvar ){
      stop("length of vector var.coef.zero does not match the number of regression coefficients, the length must be ",nvar)
    }

    tau.g.value = 1/var.coef.zero
  } else {
    tau.g.value = 1/rep(1,nvar)
  }

  if (!missing(coef.nonzero)){
    if( length(coef.nonzero) != nvar ){
      stop("length of vector coef.nonzero does not match the number of regression coefficients, the length must be ",nvar)
    }
    mu.b.value = coef.nonzero
  } else {
    mu.b.value = rep(0,nvar)
  }

  if (!missing(coef.zero)){
    if( length(coef.zero) != nvar ){
      stop("length of vector coef.zero does not match the number of regression coefficients, the length must be ",nvar)
    }
    mu.g.value = coef.zero
  } else {
    mu.g.value = rep(0,nvar)
  }

  if (iter.update < 3){
    stop("the number of iteration updates at least 3 times")
  }

  #SAMPLED AREA
  if (!any(is.na(formuladata[,1]))){

    formuladata <- as.matrix(na.omit(formuladata))
    x <- model.matrix(formula,data = as.data.frame(formuladata))
    n <- nrow(formuladata)

    mu.b = mu.b.value
    mu.g = mu.g.value

    tau.b = tau.b.value
    tau.g = tau.g.value

    tau.ua = tau.ub = 1
    tau.va = tau.vb = 1

    a.var.u = a.var.v = 1

    tau.aa = tau.ab = 0.001
    tau.ba = tau.bb = 0.001

    for(iter in 1:iter.update){

      dat <- list("n"=n, "nvar"=nvar, "y"=formuladata[,1], "x"=as.matrix(x[,-1]),
                  "mu.b"=mu.b, "mu.g"=mu.g, "tau.b"=tau.b, "tau.g"=tau.g,
                  "tau.ua"=tau.ua, "tau.ub"=tau.ub,
                  "tau.va"=tau.va, "tau.vb"=tau.vb,
                  "tau.aa"=tau.aa, "tau.ab"=tau.ab,
                  "tau.ba"=tau.ba, "tau.bb"=tau.bb)

      inits <- list(u=rep(0,n), v=rep(0,n), b=mu.b, g=mu.g,
                    tau.u=tau.u, tau.v=tau.v, tau.pa=0.001, tau.pb=0.001)

      cat("model{
					for (i in 1:n) {
							#Likelihood
							y[i] ~ dnegbin(p.zinb[i],r[i])
							p.zinb[i] <- mu.eff[i]*phi[i]/(1+mu.eff[i]*phi[i])
							r[i] <- pow(mu.eff[i],2)*phi[i]
							mu.eff[i] <- (1-pi[i])*mu[i]

							#Regression Model
							log(mu[i]) <- b[1] + sum(b[2:nvar]*x[i,]) + u[i]
							zero[i] ~ dbern(pi[i])
							logit(pi[i]) <- g[1] + sum(g[2:nvar]*x[i,]) + v[i]

							#Random effect area
							u[i] ~ dnorm(0,tau.u)
							v[i] ~ dnorm(0,tau.v)

							phi[i] ~ dgamma(tau.pa,tau.pb)
					}

					#Priors
					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					    g[k] ~ dnorm(mu.g[k],tau.g[k])
					}

					tau.u ~ dgamma(tau.ua,tau.ub)
					tau.v ~ dgamma(tau.va,tau.vb)

					a.var.u <- 1 / tau.u
					a.var.v <- 1 / tau.v

					tau.pa ~ dgamma(tau.aa,tau.ab)
					tau.pb ~ dgamma(tau.ba,tau.bb)

			}", file="saeHBzinb.txt")

      jags.m <- jags.model(file = "saeHBzinb.txt", data=dat, inits=inits, n.chains=1, n.adapt=500)
      file.remove("saeHBzinb.txt")
      params <- c("mu.eff","a.var.u","a.var.v","b","g","tau.pa","tau.pb","tau.u","tau.v")
      samps  <- coda.samples(jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start = burn.in+1, end =iter.mcmc)
      result_samps = summary(samps1)

      a.var.u = result_samps$statistics[1]
      a.var.v = result_samps$statistics[2]

      beta = result_samps$statistics[3:(nvar+2),1:2]
      gama = result_samps$statistics[(nvar+3):(2*nvar+2),1:2]

      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
        mu.g[i]  = gama[i,1]
        tau.g[i] = 1/(gama[i,2]^2)
      }

      tau.aa = result_samps$statistics[2*nvar+n+3,1]^2/result_samps$statistics[2*nvar+n+3,2]^2
      tau.ab = result_samps$statistics[2*nvar+n+3,1]/result_samps$statistics[2*nvar+n+3,2]^2

      tau.ba = result_samps$statistics[2*nvar+n+4,1]^2/result_samps$statistics[2*nvar+n+4,2]^2
      tau.bb = result_samps$statistics[2*nvar+n+4,1]/result_samps$statistics[2*nvar+n+4,2]^2

      tau.ua = result_samps$statistics[2*nvar+n+5,1]^2/result_samps$statistics[2*nvar+n+5,2]^2
      tau.ub = result_samps$statistics[2*nvar+n+5,1]/result_samps$statistics[2*nvar+n+5,2]^2

      tau.va = result_samps$statistics[2*nvar+n+6,1]^2/result_samps$statistics[2*nvar+n+6,2]^2
      tau.vb = result_samps$statistics[2*nvar+n+6,1]/result_samps$statistics[2*nvar+n+6,2]^2

    }
    result_samps = summary(samps1)

    b.varnames <- list()
    g.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i]  <- str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
      idx.g.varnames <- as.character(i-1)
      g.varnames[i]  <- str_replace_all(paste("g[",idx.g.varnames,"]"),pattern=" ", replacement="")
    }

    result_mcmc_beta <- samps1[,c(3:(nvar+2))]
    result_mcmc_gama <- samps1[,c((nvar+3):(2*nvar+2))]
    colnames(result_mcmc_beta[[1]]) <- b.varnames
    colnames(result_mcmc_gama[[1]]) <- g.varnames

    a.var.u = result_samps$statistics[1,1:2]
    a.var.v = result_samps$statistics[2,1:2]
    a.var   = as.data.frame(rbind(a.var.u,a.var.v))

    beta = result_samps$statistics[3:(nvar+2),1:2]
    gama = result_samps$statistics[(nvar+3):(2*nvar+2),1:2]
    rownames(beta) <- b.varnames
    rownames(gama) <- g.varnames
    koef1 = as.data.frame(rbind(beta,gama))

    mu.eff      = result_samps$statistics[(2*nvar+3):(2*nvar+n+2),1:2]
    Estimation1 = as.data.frame(mu.eff)

    Quantiles <- as.data.frame(result_samps$quantiles)
    q_beta = (Quantiles[3:(nvar+2),])
    q_gama = (Quantiles[(nvar+3):(2*nvar+2),])
    rownames(q_beta) <- b.varnames
    rownames(q_gama) <- g.varnames
    koef2 = as.data.frame(rbind(q_beta,q_gama))

    q_mu.eff    = (Quantiles[(2*nvar+3):(2*nvar+n+2),])
    Estimation2 = as.data.frame(q_mu.eff)

    Estimation = cbind(Estimation1,q_mu.eff)
    koef       = cbind(koef1,koef2)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")
    colnames(koef)       <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")

  }

  else {
    #NONSAMPLED AREA
    formuladata <- as.data.frame(formuladata)

    #x <- as.matrix(formuladata[,2:nvar])
    n <- nrow(formuladata)

    mu.b = mu.b.value
    mu.g = mu.g.value

    tau.b = tau.b.value
    tau.g = tau.g.value

    tau.ua = tau.ub = 1
    tau.va = tau.vb = 1

    a.var.u = a.var.v = 1

    tau.aa = tau.ab = 0.001
    tau.ba = tau.bb = 0.001

    formuladata$idx <- rep(1:n)
    data_sampled    <- na.omit(formuladata)
    data_nonsampled <- formuladata[-data_sampled$idx,]

    rc = data_nonsampled$idx
    n1 = nrow(data_sampled)
    n2 = nrow(data_nonsampled)

    for(iter in 1:iter.update){
      dat <- list("n1"=n1, "n2"=n2, "nvar"=nvar, "y_sampled"=data_sampled[,1],
                  "x_sampled"=as.matrix(data_sampled[,2:nvar]), "x_nonsampled"=as.matrix(data_nonsampled[,2:nvar]),
                  "mu.b"=mu.b, "mu.g"=mu.g, "tau.b"=tau.b, "tau.g"=tau.g,
                  "tau.ua"=tau.ua, "tau.ub"=tau.ub,
                  "tau.va"=tau.va, "tau.vb"=tau.vb,
                  "tau.aa"=tau.aa, "tau.ab"=tau.ab,
                  "tau.ba"=tau.ba, "tau.bb"=tau.bb)

      inits <- list(u=rep(0,n1), v=rep(0,n1),
                    uT=rep(0,n2), vT=rep(0,n2),
                    b=mu.b, g=mu.g, tau.u=tau.u,
                    tau.v=tau.v, tau.pa=0.001, tau.pb=0.001)

      cat("model {
					for (i in 1:n1) {
							#Likelihood
							y_sampled[i] ~ dnegbin(p.zinb[i],r[i])
							p.zinb[i] <- mu.eff[i]*phi[i]/(1+mu.eff[i]*phi[i])
							r[i] <- pow(mu.eff[i],2)*phi[i]
							mu.eff[i] <- (1-pi[i])*mu[i]

							#Regression Model
							log(mu[i]) <- b[1] + sum(b[2:nvar]*x_sampled[i,]) + u[i]
							zero[i] ~ dbern(pi[i])
							logit(pi[i]) <- g[1] + sum(g[2:nvar]*x_sampled[i,]) + v[i]

							#Random effect area
							u[i] ~ dnorm(0,tau.u)
							v[i] ~ dnorm(0,tau.v)

							phi[i] ~ dgamma(tau.pa,tau.pb)
					}

					for (j in 1:n2) {
							#Model Regresi
							zeroT[j] ~ dbern(piT[j])
							log(muT[j]) <- mu.b[1] + sum(mu.b[2:nvar]*x_nonsampled[j,]) + u[j]
							logit(piT[j]) <- mu.g[1] + sum(mu.g[2:nvar]*x_nonsampled[j,]) + v[j]

							mu.eff.nonsampled[j] <- (1-piT[j])*muT[j]

							#Random effect area
							uT[j] ~ dnorm(0,tau.u)
							vT[j] ~ dnorm(0,tau.v)
          }

					#Priors
					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					    g[k] ~ dnorm(mu.g[k],tau.g[k])
					}

					tau.u ~ dgamma(tau.ua,tau.ub)
					tau.v ~ dgamma(tau.va,tau.vb)

					a.var.u <- 1 / tau.u
					a.var.v <- 1 / tau.v

					tau.pa ~ dgamma(tau.aa,tau.ab)
					tau.pb ~ dgamma(tau.ba,tau.bb)

			}", file="saeHBZINB.txt")

      jags.m <- jags.model( file = "saeHBZINB.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("saeHBZINB.txt")
      params <- c("mu.eff","mu.eff.nonsampled","a.var.u","a.var.v","b","g","tau.pa","tau.pb","tau.u","tau.v")
      samps  <- coda.samples(jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start = burn.in+1, end =iter.mcmc )
      result_samps = summary(samps1)

      a.var.u = result_samps$statistics[1]
      a.var.v = result_samps$statistics[2]

      beta = result_samps$statistics[3:(nvar+2),1:2]
      gama = result_samps$statistics[(nvar+3):(2*nvar+2),1:2]

      for (i in 1:nvar){
        mu.g[i]  = gama[i,1]
        mu.b[i]  = beta[i,1]
        tau.g[i] = 1/(gama[i,2]^2)
        tau.b[i] = 1/(beta[i,2]^2)
      }

      tau.aa = result_samps$statistics[2*nvar+n+3,1]^2/result_samps$statistics[2*nvar+n+3,2]^2
      tau.ab = result_samps$statistics[2*nvar+n+3,1]/result_samps$statistics[2*nvar+n+3,2]^2

      tau.ba = result_samps$statistics[2*nvar+n+4,1]^2/result_samps$statistics[2*nvar+n+4,2]^2
      tau.bb = result_samps$statistics[2*nvar+n+4,1]/result_samps$statistics[2*nvar+n+4,2]^2

      tau.ua = result_samps$statistics[2*nvar+n+5,1]^2/result_samps$statistics[2*nvar+n+5,2]^2
      tau.ub = result_samps$statistics[2*nvar+n+5,1]/result_samps$statistics[2*nvar+n+5,2]^2

      tau.va = result_samps$statistics[2*nvar+n+6,1]^2/result_samps$statistics[2*nvar+n+6,2]^2
      tau.vb = result_samps$statistics[2*nvar+n+6,1]/result_samps$statistics[2*nvar+n+6,2]^2
    }
    result_samps = summary(samps1)

    b.varnames <- list()
    g.varnames <- list()

    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i]  <- str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
      idx.g.varnames <- as.character(i-1)
      g.varnames[i]  <- str_replace_all(paste("g[",idx.g.varnames,"]"),pattern=" ", replacement="")
    }

    result_mcmc_beta <- samps1[,c(3:(nvar+2))]
    result_mcmc_gama <- samps1[,c((nvar+3):(2*nvar+2))]
    colnames(result_mcmc_beta[[1]]) <- b.varnames
    colnames(result_mcmc_gama[[1]]) <- g.varnames

    a.var.u = result_samps$statistics[1,1:2]
    a.var.v = result_samps$statistics[2,1:2]
    a.var   = as.data.frame(rbind(a.var.u,a.var.v))

    beta = result_samps$statistics[3:(nvar+2),1:2]
    gama = result_samps$statistics[(nvar+3):(2*nvar+2),1:2]
    rownames(beta) <- b.varnames
    rownames(gama) <- g.varnames
    koef1 = as.data.frame(rbind(beta,gama))

    mu.eff            = result_samps$statistics[(2*nvar+3):(2*nvar+n1+2),1:2]
    mu.eff.nonsampled = result_samps$statistics[(2*nvar+n1+3):(2*nvar+n+2),1:2]

    Estimation1       = matrix(rep(0,n),n,2)
    Estimation1[rc,]  = mu.eff.nonsampled
    Estimation1[-rc,] = mu.eff
    Estimation1       = as.data.frame(Estimation1)


    Quantiles       <- as.data.frame(result_samps$quantiles)
    q_beta          <- (Quantiles[3:(nvar+2),])
    q_gama          <- (Quantiles[(nvar+3):(2*nvar+2),])
    q_mu            <- (Quantiles[(2*nvar+3):(2*nvar+n1+2),])
    q_mu.nonsampled <- (Quantiles[(2*nvar+n1+3):(2*nvar+n+2),])
    q_Estimation    <- matrix(0,n,5)

    for (i in 1:5){
      q_Estimation[rc,i]  <- q_mu.nonsampled[,i]
      q_Estimation[-rc,i] <- q_mu[,i]
    }

    rownames(q_beta) <- b.varnames
    rownames(q_gama) <- g.varnames
    koef2=as.data.frame(rbind(q_beta,q_gama))

    Estimation = data.frame(Estimation1,q_Estimation)
    koef       = cbind(koef1,koef2)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")
    colnames(koef)       <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")
  }

  result$Est         = Estimation
  result$refVar      = a.var
  result$coefficient = koef
  result$plot.beta   = list(graphics.off() ,par(mar=c(2,2,2,2)),autocorr.plot(result_mcmc_beta,col="brown2",lwd=2),plot(result_mcmc_beta,col="brown2",lwd=2))
  result$plot.gamma  = list(graphics.off() ,par(mar=c(2,2,2,2)),autocorr.plot(result_mcmc_gama,col="brown2",lwd=2),plot(result_mcmc_gama,col="brown2",lwd=2))
  return(result)
}
