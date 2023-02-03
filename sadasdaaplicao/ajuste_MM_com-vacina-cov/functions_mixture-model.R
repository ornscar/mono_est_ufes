############################## Functions ###########################################

### Function to return the theta parameter

### Function to return the value of mu_q  of the susceptibles as a function of quantile 0<q<1, alpha in R, lambda>0 and theta>0
quantil_mu_MM <- function(q = 0.5,
                          alpha = 0.5,
                          lambda = 2,
                          theta = 2) {
  out <-
    (1 / alpha) * (log(1 - (alpha / lambda) * log(1 - q ^ (1 / theta))))
  return(out)
}

### Function to return the theta parameter, as a function of mu_q (mu) >0, alpha in R, lambda>0 and 0<q<1
theta_par_MM <- function(mu = 0.5,
                         alpha = 0.5,
                         lambda = 2,
                         q = 0.5) {
  out <-
    log(q) / (log(1 - exp((lambda / alpha) * (1 - exp(
      alpha * mu)))))
  return(out)
}

### Function to return the lambda parameter, as a function of mu_q (mu) >0, alpha in R, lambda>0 and 0<q<1
lambda_par_MM <- function(mu = 0.5,
                          alpha = 0.5,
                          theta = 2,
                          q = 0.5) {
  out <- (alpha * log(1 - q ^ (1 / theta))) / (1 - exp(alpha * mu))
  return(out)
}

## GG density function
dGG_MM <- function(y,
                   alpha,
                   lambda,
                   theta, 
                   log = FALSE) {
  logfy  <-
    log(theta) + log(lambda) + (theta - 1) * log(1 - exp(-lambda *
                                                           (exp(alpha * y) - 1) / alpha)) - (-alpha ^ 2 * y + lambda * exp(alpha *
                                                                                                                             y) - lambda) / alpha
  ifelse(log, return(logfy), return(exp(logfy)))
}

#----------------------------------------------------------------------------------------
## GG distribution function
pGG_MM <-
  function(y,
           alpha,
           lambda,
           theta, 
           lower.tail = TRUE,
           log.p = FALSE) {
    cdf   <-
      (1 - exp(-lambda * (exp(alpha * y) - 1) / alpha)) ^ theta
    
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

#

################################ Bayesian ##############################################
#J is the number of covariates (do not include the intercept)
#guess is the initial values for the parameters - alpha,lambda,beta
bayesGE_MM  <-
  function(linear_pred1,
           linear_pred2,
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
    parm.names <- LaplacesDemon::as.parm.names(list(
      lambda = 0,
      beta = rep(0, J1),
      gama = rep(0, J2)
    )) 
    pos.lambda   <- grep("lambda", parm.names)
    pos.beta <- grep("beta", parm.names)
    pos.gama <- grep("gama", parm.names)
    
    MyData <- list(
      linear_pred1 = linear_pred1,
      linear_pred2 = linear_pred2,
      data = data, 
      d = d,
      q = q,
      log = TRUE,
      mon.names = mon.names,
      parm.names = parm.names, 
      pos.lambda = pos.lambda,
      pos.beta = pos.beta,
      pos.gama = pos.gama,
      N = 1
    )
    
    Model <- function(parm, Data) {
      ### Parameters
      lambda  <- parm[Data$pos.lambda]
      # parm[Data$pos.lambda] <- lambda
      beta <-  parm[Data$pos.beta]
      gama <-  parm[Data$pos.gama]
      ### Log(Prior Densities)
      lambda.prior <-
        stats::dnorm(lambda,
                      mean = 0,
                     sd = sqrt(100),
                      log = TRUE) #mean 1 and var 100
      beta.prior <- sum(LaplacesDemon::dnormv(
        beta,
        mean = 0,
        var = 100,
        log = TRUE
      ))
      gama.prior <- sum(LaplacesDemon::dnormv(
        gama,
        mean = 0,
        var = 100,
        log = TRUE
      ))
      
      ### Log-Likelihood
      LL <- likGG_MM_GE(  
        par = c(lambda, beta, gama),
        linear_pred1 = Data$linear_pred1,
        linear_pred2 = Data$linear_pred2,
        q = Data$q,
        d = Data$d,
        data = Data$data, 
        log = Data$log
      )

      ### Log-Posterior
      LP <- LL + lambda.prior + beta.prior + gama.prior
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
    
    # Fit <- LaplacesDemon::LaplacesDemon(
    #   Model = Model,
    #   Data = MyData,
    #   Initial.Values = guess,
    #   Covar = NULL,
    #   Iterations = burn + jump * iter,
    #   Status = 20000,
    #   Thinning = jump,
    #   Algorithm = "AM",
    #   Specs = list(Adaptive = 500, Periodicity = 100),
    #   Debug = list(
    #     DB.chol = FALSE,
    #     DB.eigen = FALSE,
    #     DB.MCSE = FALSE,
    #     DB.Model = FALSE
    #   )
    # )
    
    ######################  Robust Adaptive Metropolis  #######################
    #Fit <- LaplacesDemon(
    # Model,
    # Data = MyData,
    # Initial.Values = guess,
    # Covar = NULL,
    # Iterations = burn,
    # Status = burn,
    # Thinning = 2, 
    # Algorithm="RAM",
    # Specs=list(alpha.star=0.234, B=NULL, Dist="N",
    #     gamma=0.66, n=0))
    # initial <- Fit$Posterior1[dim(Fit$Posterior1)[1], ]
    
    #Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values = initial,
    #     Covar = NULL, Iterations = jump * iter, 
    # Status = jump * iter, 
    # Thinning = jump,
    #     Algorithm="RWM", Specs=NULL)
    
    
    Fit <- LaplacesDemon(
      Model,
      Data = MyData,
      Initial.Values = guess,
      Covar = NULL,
      Iterations = burn + jump * iter,
      Status = burn + jump * iter,
      Thinning = jump,
      Algorithm = "AMWG",
      Specs = list(B = NULL, n = 0, Periodicity = 200)
    )
    
    Posterior <-
      Fit$Posterior1[(length(seq(jump, burn, jump)) + 1):length(Fit$Posterior1[, 1]), ]
    
    Posterior[,1] <- exp(Posterior[,1]) 
    
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
        "rec_jump" = jump_rec
      )
    return(output)
  }


################################ Bayesian ##############################################
#J is the number of covariates (do not include the intercept)
#guess is the initial values for the parameters - alpha,lambda,beta
bayesGE_MM_ergodicas  <-
  function(linear_pred1,
           linear_pred2,
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
    parm.names <- LaplacesDemon::as.parm.names(list(
      lambda = 0,
      beta = rep(0, J1),
      gama = rep(0, J2)
    )) 
    pos.lambda   <- grep("lambda", parm.names)
    pos.beta <- grep("beta", parm.names)
    pos.gama <- grep("gama", parm.names)
    
    MyData <- list(
      linear_pred1 = linear_pred1,
      linear_pred2 = linear_pred2,
      data = data, 
      d = d,
      q = q,
      log = TRUE,
      mon.names = mon.names,
      parm.names = parm.names, 
      pos.lambda = pos.lambda,
      pos.beta = pos.beta,
      pos.gama = pos.gama,
      N = 1
    )
    
    Model <- function(parm, Data) {
      ### Parameters
      lambda  <- parm[Data$pos.lambda]
      # parm[Data$pos.lambda] <- lambda
      beta <-  parm[Data$pos.beta]
      gama <-  parm[Data$pos.gama]
      ### Log(Prior Densities)
      lambda.prior <-
        stats::dnorm(lambda,
                     mean = 0,
                     sd = sqrt(100),
                     log = TRUE) #mean 1 and var 100
      beta.prior <- sum(LaplacesDemon::dnormv(
        beta,
        mean = 0,
        var = 100,
        log = TRUE
      ))
      gama.prior <- sum(LaplacesDemon::dnormv(
        gama,
        mean = 0,
        var = 100,
        log = TRUE
      ))
      
      ### Log-Likelihood
      LL <- likGG_MM_GE(  
        par = c(lambda, beta, gama),
        linear_pred1 = Data$linear_pred1,
        linear_pred2 = Data$linear_pred2,
        q = Data$q,
        d = Data$d,
        data = Data$data, 
        log = Data$log
      )
      
      ### Log-Posterior
      LP <- LL + lambda.prior + beta.prior + gama.prior
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
    
    # Fit <- LaplacesDemon::LaplacesDemon(
    #   Model = Model,
    #   Data = MyData,
    #   Initial.Values = guess,
    #   Covar = NULL,
    #   Iterations = burn + jump * iter,
    #   Status = 20000,
    #   Thinning = jump,
    #   Algorithm = "AM",
    #   Specs = list(Adaptive = 500, Periodicity = 100),
    #   Debug = list(
    #     DB.chol = FALSE,
    #     DB.eigen = FALSE,
    #     DB.MCSE = FALSE,
    #     DB.Model = FALSE
    #   )
    # )
    
    ######################  Robust Adaptive Metropolis  #######################
    #Fit <- LaplacesDemon(
    # Model,
    # Data = MyData,
    # Initial.Values = guess,
    # Covar = NULL,
    # Iterations = burn,
    # Status = burn,
    # Thinning = 2, 
    # Algorithm="RAM",
    # Specs=list(alpha.star=0.234, B=NULL, Dist="N",
    #     gamma=0.66, n=0))
    # initial <- Fit$Posterior1[dim(Fit$Posterior1)[1], ]
    
    #Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values = initial,
    #     Covar = NULL, Iterations = jump * iter, 
    # Status = jump * iter, 
    # Thinning = jump,
    #     Algorithm="RWM", Specs=NULL)
    
    
    Fit <- LaplacesDemon(
      Model,
      Data = MyData,
      Initial.Values = guess,
      Covar = NULL,
      Iterations = burn + jump * iter,
      Status = burn + jump * iter,
      Thinning = jump,
      Algorithm = "AMWG",
      Specs = list(B = NULL, n = 0, Periodicity = 200)
    )
    
    Posterior <- Fit$Posterior1
    Posterior[,1] <- exp(Fit$Posterior1[,1]) 
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
        "rec_jump" = jump_rec
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

#################### CPO MM ##########################################
cpo_mm <- function(Posterior, linear_pred1, linear_pred2, dados, d, q) {
  #média do log do CPO
  J1 <- dim(stats::model.matrix(linear_pred1, dados))[2] 
  J2 <- dim(stats::model.matrix(linear_pred2, dados))[2] 
  
  alpha <- Posterior[,1]
  lambda <- Posterior[,2]
  beta <- Posterior[,3:(2+J1)]
  gama <- Posterior[,(3 + J1):(2 + J1 + J2)] 
  
  cens <- get(d, dados)
  temp <-
    as.numeric(stats::model.extract(stats::model.frame(linear_pred1, dados), 'response'))
  X <- stats::model.matrix(linear_pred1, dados)
  X2 <- stats::model.matrix(linear_pred2, dados) 
  
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
    
    theta <- theta_par_MM(mu = mu,
                          alpha = alpha,
                          lambda = lambda,
                          q = q)
    
    if(dim(X2)[2] > 1) {
      x_p <- do.call(rbind, lapply(apply(gama, 1, list), function(g) {
        tcrossprod(X2[i,], t(g[[1]]))
      }))
    } else {
      x_p <- X2[i,]*gama
    }
    p0 <- exp(x_p)/(1+exp(x_p))
    
    out <-
      sum(1 / ((
        ((1 - p0)*dGG_MM(
          t,
          alpha = alpha,
          lambda = lambda,
          theta = theta,
          log = FALSE
        ) )^ d
      ) *
        (
          (p0 + (1 - p0) * pGG_MM(
            t,
            alpha = alpha,
            lambda = lambda,
            theta = theta,
            lower.tail = FALSE,
            log.p = FALSE
          )) ^ (1 - d)
        )))
    
    CPO <- (out / dim(Posterior)[1]) ^ (-1)
    
    aux <- aux + log(CPO)
    
  }
  B <- aux #/dim(dados)[1]
  return(B)
}

