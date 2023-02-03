############################## Functions ###########################################

#-------------------------------------------------------------------------------------------
## Generalized Gompertz distribution
## El-Gohary A, Alshamrani A, Al-Otaibi AN. The generalized Gompertz distribution. Applied Mathematical Modelling. 2013;37:13-24.
## where alpha in R, lambda > 0 and theta > 0
#-------------------------------------------------------------------------------------------

### Function to return the value of mu_q  of the susceptibles as a function of quantile 0<q<1, alpha in R, lambda>0 and theta>0
quantil.mu <- function(q = 0.5,
                       alpha = -1,
                       lambda = 2,
                       theta = 2) {
  out <-
    (1 / alpha) * (log(lambda - alpha * log(1 - q ^ (1 / theta) * (
      1 - exp(lambda / alpha)
    ))) - log(lambda))
  return(out)
}

### Function to return 0<q<1 as a function of mu_q (mu) >0, alpha in R, lambda>0 and theta>0
quantil.q <- function(mu = 0.5,
                      alpha = -1,
                      lambda = 2,
                      theta = 2) {
  out <-
    exp(theta * log((1 - exp((lambda - exp(alpha * mu + log(lambda))) / alpha
    )) / (1 - exp(lambda / alpha))))
  return(out)
}

### Function to return the theta parameter, as a function of mu_q (mu) >0, alpha in R, lambda>0 and 0<q<1
theta.par <- function(mu = 0.5,
                      alpha = -1,
                      lambda = 2,
                      q = 0.5) {
  out <-
    log(q) / (log(1 - exp((
      lambda - exp(alpha * mu + log(lambda))
    ) / alpha)) - log(1 - exp(lambda / alpha)))
  return(out)
}

### Function that calculates the cure fraction################
p0.dGG <- function(alpha, lambda, theta) {
  out <- 1 - (1 - exp(lambda / alpha)) ^ theta
  return(out)
}

# dGGteste <- function(y, alpha=-1, lambda=2, theta=2){
#     fy = theta*lambda*exp(alpha*y-(lambda/alpha)*(exp(alpha*y)-1))*(1-exp(-(lambda/alpha)*(exp(alpha*y)-1)))^(theta-1)
#     return(fy)
# }
#
#

dGG <- function(y, alpha, lambda, mu, q, log = FALSE) {
  if (any(lambda < 0))
    stop(paste("lambda must be positive", "\n", ""))
  
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  
  if (any(y <= 0))
    stop(paste("y must be greater than 0 ", "\n", ""))
  
  if (any(q < 0) | any(q > 1))
    stop(paste("q must be in the interval (0,1) ", "\n", ""))
  
  theta <-
    -log(q) / (log(1 - exp(lambda / alpha)) - log(1 - exp(-(lambda / alpha) *
                                                            (exp(
                                                              alpha * mu
                                                            ) - 1))))
  logfy  = log(theta) + log(lambda) + (theta - 1) * log(1 - exp(-lambda *
                                                                  (exp(alpha * y) - 1) / alpha)) - (-alpha ^ 2 * y + lambda * exp(alpha *
                                                                                                                                    y) - lambda) / alpha
  
  if (log == FALSE)
    fy = exp(logfy)
  else
    fy = logfy
  fy
}

#----------------------------------------------------------------------------------------
pGG <-
  function(y,
           alpha,
           lambda,
           mu,
           q,
           lower.tail = TRUE,
           log.p = FALSE) {
    if (any(lambda < 0))
      stop(paste("lambda must be positive", "\n", ""))
    
    if (any(mu < 0))
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(y <= 0))
      stop(paste("y must be greater than 0 ", "\n", ""))
    
    if (any(q < 0) | any(q > 1))
      stop(paste("q must be in the interval (0,1) ", "\n", ""))
    
    theta <-
      -log(q) / (log(1 - exp(lambda / alpha)) - log(1 - exp(-(lambda / alpha) *
                                                              (exp(
                                                                alpha * mu
                                                              ) - 1))))
    
    cdf   <- (1 - exp(-lambda * (exp(alpha * y) - 1) / alpha)) ^ theta
    
    if (lower.tail == TRUE)
      cdf <- cdf
    else
      cdf <- 1 - cdf
    if (log.p == FALSE)
      cdf <- cdf
    else
      cdf <- log(cdf)
    cdf
  }

####################### Likelihood function #################################################
likGG <- function(par, q, data, linear_pred1, linear_pred2, linear_pred3, d, log = FALSE) {
  J1 <- dim(stats::model.matrix(linear_pred1, data))[2] 
  J2 <- dim(stats::model.matrix(linear_pred2, data))[2]
  J3 <- dim(stats::model.matrix(linear_pred3, data))[2] 
  
  beta <- par[1:J1]
  gama <- par[(1 + J1):(J1 + J2)]
  vi <- par[(J1 + J2 + 1):(J1 + J2 + J3)]
  
  cens <- get(d, data)
  t <-
    as.numeric(stats::model.extract(stats::model.frame(linear_pred1, data), 'response'))
  X <- stats::model.matrix(linear_pred1, data)
  mu <- exp(tcrossprod(X, t(beta)))
  
  X2 <- stats::model.matrix(linear_pred2, data)
  alpha <- -exp(tcrossprod(X2, t(gama)))
  
  X3 <- stats::model.matrix(linear_pred3, data)
  lambda <- exp(tcrossprod(X3, t(vi)))
  
  out <-
    sum(
      cens * dGG(
        t,
        alpha = alpha,
        lambda = lambda,
        mu = mu,
        q = q,
        log = TRUE
      ) +
        (1 - cens) * pGG(
          t,
          alpha = alpha,
          lambda = lambda,
          mu = mu,
          q = q,
          lower.tail = FALSE,
          log.p = TRUE
        )
    )
  
  if (log == TRUE)
    lik <- out
  else
    lik <- exp(out)
  return(lik)

}

#########################################################################################

################################ Bayesian ##############################################
#J is the number of covariates (do not include the intercept)
#guess is the initial values for the parameters - alpha,lambda,beta
bayesGG  <- function(linear_pred1,
                     linear_pred2,
                     linear_pred3,
                     data,
                     q,
                     d,
                     iter = 1000,
                     burn,
                     jump,
                     guess) {
  mon.names <- c("LP")
  J1 <- dim(stats::model.matrix(linear_pred1, data))[2] 
  J2 <- dim(stats::model.matrix(linear_pred2, data))[2] 
  J3 <- dim(stats::model.matrix(linear_pred3, data))[2] 
  parm.names <-
    as.parm.names(list(
      beta = rep(0, J1),
      gama = rep(0, J2), 
      vi = rep(0, J3)
    ))
  pos.beta <- grep("beta", parm.names)
  pos.gama <- grep("gama", parm.names)
  pos.vi <- grep("vi", parm.names)
  
  MyData <- list(
    linear_pred1 = linear_pred1,
    linear_pred2 = linear_pred2,
    linear_pred3 = linear_pred3,
    data = data,
    d = d,
    q = q,
    log = TRUE,
    mon.names = mon.names,
    parm.names = parm.names, 
    pos.beta = pos.beta,
    pos.gama = pos.gama,
    pos.vi = pos.vi,
    N = 1
  )
  
  Model <- function(parm, Data) {
    ### Parameters
    beta <-  parm[Data$pos.beta]
    gama <-  parm[Data$pos.gama]
    vi <-  parm[Data$pos.vi]
    ### Log(Prior Densities)
    beta.prior <- sum(dnormv(
      beta,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    gama.prior <- sum(dnormv(
      gama,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    vi.prior <- sum(dnormv(
      vi,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    ### Log-Likelihood
    LL <-
      likGG(
        par = c(beta, gama, vi),
        q = Data$q,
        d = Data$d,
        data = Data$data,
        linear_pred1 = Data$linear_pred1,
        linear_pred2 = Data$linear_pred2,
        linear_pred3 = Data$linear_pred3,
        log = Data$log
      ) 
    ### Log-Posterior
    LP <- LL + beta.prior + gama.prior + vi.prior
    Modelout <-
      list(
        LP = LP,
        Dev = -2 * LL,
        Monitor = LP,
        yhat = 1,
        parm = parm
      )
    
    return(Modelout)
  }
  
  Fit <- LaplacesDemon(
    Model = Model,
    Data = MyData,
    Initial.Values = guess,
    Covar = NULL,
    Iterations = burn + jump * iter,
    Status = 20000,
    Thinning = jump,
    Algorithm = "AM",
    Specs = list(Adaptive = 500, Periodicity = 100)
  )
  
  #######################  Metropolis-within-Gibbs  #########################
  #Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
  #     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
  #     Algorithm="MWG", Specs=list(B=NULL))
  ###################  Adaptive Metropolis-within-Gibbs  ####################
  #Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
  #     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
  #     Algorithm="AMWG", Specs=list(B=NULL, n=0, Periodicity=50))
  ###################  Adaptive Metropolis ####################
  # Fit <- LaplacesDemon(Model=Model,
  #                      Data=MyData,
  #                      Initial.Values=guess,
  #                      Covar=NULL,
  #                      Iterations=burn+jump*n.size,
  #                      Status=20000,
  #                      Thinning=jump,
  #                      Algorithm="AM",
  #                      Specs=list(Adaptive=500, Periodicity=100))
  
  ######################## AHMC #######################################
  # Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values=guess,
  #                      Covar=NULL,
  #                      Iterations=burn+jump*n.size,
  #                      Status=20000, Thinning=1,
  #                      Algorithm="AHMC",
  #                      Specs=list(epsilon=0.02, L=2, m=NULL,L=2, Periodicity=10))
  #
  
  Posterior <-
    Fit$Posterior1[(length(seq(jump, burn, jump)) + 1):length(Fit$Posterior1[, 1]), ]
  burn_rec <- Fit$Rec.BurnIn.Thinned
  jump_rec <- Fit$Rec.Thinning
  AR <- Fit$Acceptance.Rate
  DIC <- Fit$DIC1[1]
  
  output <-
    list(
      "post" = Posterior,
      "AR" = AR,
      "DIC" = DIC,
      "rec_burnin" = burn_rec,
      "rec_jump" = jump_rec, 
      "q" = q
    )
  return(output)
}

################################ Bayesian ##############################################
#J is the number of covariates (do not include the intercept)
#guess is the initial values for the parameters - alpha,lambda,beta
bayesGG_ergodicas  <- function(linear_pred1,
                               linear_pred2,
                               linear_pred3,
                               data,
                               q,
                               d,
                               iter = 1000,
                               burn,
                               jump,
                               guess) {
  mon.names <- c("LP")
  J1 <- dim(stats::model.matrix(linear_pred1, data))[2] 
  J2 <- dim(stats::model.matrix(linear_pred2, data))[2] 
  J3 <- dim(stats::model.matrix(linear_pred3, data))[2] 
  parm.names <-
    as.parm.names(list(
      beta = rep(0, J1),
      gama = rep(0, J2), 
      vi = rep(0, J3)
    ))
  pos.beta <- grep("beta", parm.names)
  pos.gama <- grep("gama", parm.names)
  pos.vi <- grep("vi", parm.names)
  
  MyData <- list(
    linear_pred1 = linear_pred1,
    linear_pred2 = linear_pred2,
    linear_pred3 = linear_pred3,
    data = data,
    d = d,
    q = q,
    log = TRUE,
    mon.names = mon.names,
    parm.names = parm.names, 
    pos.beta = pos.beta,
    pos.gama = pos.gama,
    pos.vi = pos.vi,
    N = 1
  )
  
  Model <- function(parm, Data) {
    ### Parameters
    beta <-  parm[Data$pos.beta]
    gama <-  parm[Data$pos.gama]
    vi <-  parm[Data$pos.vi]
    ### Log(Prior Densities)
    beta.prior <- sum(dnormv(
      beta,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    gama.prior <- sum(dnormv(
      gama,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    vi.prior <- sum(dnormv(
      vi,
      mean = 0,
      var = 100,
      log = TRUE
    ))
    ### Log-Likelihood
    LL <-
      likGG(
        par = c(beta, gama, vi),
        q = Data$q,
        d = Data$d,
        data = Data$data,
        linear_pred1 = Data$linear_pred1,
        linear_pred2 = Data$linear_pred2,
        linear_pred3 = Data$linear_pred3,
        log = Data$log
      ) 
    ### Log-Posterior
    LP <- LL + beta.prior + gama.prior + vi.prior
    Modelout <-
      list(
        LP = LP,
        Dev = -2 * LL,
        Monitor = LP,
        yhat = 1,
        parm = parm
      )
    
    return(Modelout)
  }
  
  Fit <- LaplacesDemon(
    Model = Model,
    Data = MyData,
    Initial.Values = guess,
    Covar = NULL,
    Iterations = burn + jump * iter,
    Status = 20000,
    Thinning = jump,
    Algorithm = "AM",
    Specs = list(Adaptive = 500, Periodicity = 100)
  )
  
  #######################  Metropolis-within-Gibbs  #########################
  #Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
  #     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
  #     Algorithm="MWG", Specs=list(B=NULL))
  ###################  Adaptive Metropolis-within-Gibbs  ####################
  #Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
  #     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
  #     Algorithm="AMWG", Specs=list(B=NULL, n=0, Periodicity=50))
  ###################  Adaptive Metropolis ####################
  # Fit <- LaplacesDemon(Model=Model,
  #                      Data=MyData,
  #                      Initial.Values=guess,
  #                      Covar=NULL,
  #                      Iterations=burn+jump*n.size,
  #                      Status=20000,
  #                      Thinning=jump,
  #                      Algorithm="AM",
  #                      Specs=list(Adaptive=500, Periodicity=100))
  
  ######################## AHMC #######################################
  # Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values=guess,
  #                      Covar=NULL,
  #                      Iterations=burn+jump*n.size,
  #                      Status=20000, Thinning=1,
  #                      Algorithm="AHMC",
  #                      Specs=list(epsilon=0.02, L=2, m=NULL,L=2, Periodicity=10))
  #
  
  Posterior <-
    Fit$Posterior1#[(length(seq(jump,burn,jump))+1):length(Fit$Posterior1[,1]),]
  burn_rec <- Fit$Rec.BurnIn.Thinned
  jump_rec <- Fit$Rec.Thinning
  AR <- Fit$Acceptance.Rate
  DIC <- Fit$DIC1[1]
  
  output <-
    list(
      "post" = Posterior,
      "AR" = AR,
      "DIC" = DIC,
      "rec_burnin" = burn_rec,
      "rec_jump" = jump_rec, 
      "q" = q
    )
  return(output)
}


##################### Funcoes saidas ####################################################

#### funcoes para saidas com summary ###################################################
##### so um vetor
saidas <- function(cadeia) {
  out <- cbind(t(summary(cadeia)),
               as.matrix(sd(cadeia)),
               t(quantile(cadeia, probs = c(0.025, 0.975))),
               t(emp.hpd(cadeia, conf = 0.95)))
  return(out)
}

#######################################################################################
##### para matriz
saidas.matriz <- function(cadeia) {
  out <- cbind(t(apply(cadeia, 2, summary)),
               as.matrix(apply(cadeia, 2, sd)),
               t(apply(cadeia, 2, function(x)
                 quantile(x, probs = c(0.025, 0.975)))),
               t(apply(cadeia, 2, function(x)
                 emp.hpd(x, conf = 0.95))))
  out <- as.data.frame(out)
  colnames(out) <-
    c(
      "min",
      "1qt",
      "med",
      "media",
      "3qt",
      "max",
      "dp",
      "emp_2.5",
      "emp_97.5",
      "HPD_2.5",
      "HPD_97.5"
    )
  return(out)
}

#################### CPO ##########################################
cpo <- function(Posterior, linear_pred1, linear_pred2, linear_pred3, dados, d, q) {
  #média do log do CPO
  J1 <- dim(stats::model.matrix(linear_pred1, dados))[2] 
  J2 <- dim(stats::model.matrix(linear_pred2, dados))[2]
  J3 <- dim(stats::model.matrix(linear_pred3, dados))[2] 
  
  beta <- Posterior[,1:J1]
  gama <- Posterior[,(1 + J1):(J1 + J2)]
  vi <- Posterior[,(J1 + J2 + 1):(J1 + J2 + J3)]
  
  cens <- get(d, dados)
  temp <-
    as.numeric(stats::model.extract(stats::model.frame(linear_pred1, dados), 'response'))
  X <- stats::model.matrix(linear_pred1, dados)
  X2 <- stats::model.matrix(linear_pred2, dados)
  X3 <- stats::model.matrix(linear_pred3, dados)
  
  aux <- 0
  for (i in 1:dim(dados)[1]) {
    d <- cens[i]
    t <- temp[i]
    if(dim(X)[2] > 1) {
      x_beta <- do.call(rbind, lapply(apply(beta, 1, list), function(b) {
        tcrossprod(X[i,], t(b[[1]]))
      }))
    } else {
      x_beta <- X[i,]*beta
    }
    mu <- exp(x_beta)
    
    if(dim(X2)[2] > 1) {
      x_alpha <- do.call(rbind, lapply(apply(gama, 1, list), function(g) {
        tcrossprod(X2[i,], t(g[[1]]))
      }))
    } else {
      x_alpha <- X2[i,]*gama
    }
    alpha <- -exp(x_alpha)
    
    if(dim(X3)[2] > 1) {
      x_lambda <- do.call(rbind, lapply(apply(vi, 1, list), function(v) {
        tcrossprod(X3[i,], t(v[[1]]))
      }))
    } else {
      x_lambda <- X3[i,]*vi
    }
    lambda <- exp(x_lambda)
    
    out <-
      sum(1 / ((
        dGG(
          t,
          alpha = alpha,
          lambda = lambda,
          mu = mu,
          q = q,
          log = FALSE
        ) ^ d
      ) *
        (
          pGG(
            t,
            alpha = alpha,
            lambda = lambda,
            mu = mu,
            q = q,
            lower.tail = FALSE,
            log.p = FALSE
          ) ^ (1 - d)
        )))
    
    CPO <- (out / dim(Posterior)[1]) ^ (-1)
    
    aux <- aux + log(CPO)
    
  }
  B <- aux #/dim(dados)[1]
  return(B)
}