#### analises 29/12/2022

## pacotes
library(readxl)
library(janitor)
library(dplyr)
library(lubridate)
library(forcats)
library(pander)
library(skimr)
require(survival)
require(truncnorm)
require(LaplacesDemon)
require(TeachingDemos)
require(coda)
require(survival)
require(survMisc)
require(ggfortify)
library(foreign)
library(survminer)
library(maxstat)
library(ggcrqr)

############ Calling the functions #######################
source("functions_mixture-model.R")
source("diag_conver.R")
#function to get the convergence graphs

#####
#medias ergódicas média dos parâmetros
med.erg2 <- function(x)
{
  aux <- rep(0, length(x[, 1]))
  for (i in 1:length(x[, 1]))                              {
    aux[i] <-
      mean(x[i, ])
  }
  out <- rep(0, length(aux))
  out[1] <- aux[1]
  for (i in 2:length(aux))
  {
    out[i] <- (out[i - 1] * (i - 1) + aux[i]) / i
  }
  return(out)
}


########################

dados <-  readRDS("dados_modelo.rds")

dados$faixa2 <- as.factor(dados$faixa2)
dados$grupos <- as.factor(dados$grupos)
dados$obesidade <- as.factor(dados$obesidade)
dados$perd_pala <- as.factor(dados$perd_pala)
dados$saturacao <- as.factor(dados$saturacao)
dados$vacina_cov <- as.factor(dados$vacina_cov)
# dados$vacina_cov2 <- as.factor(dados$vacina_cov2)


## Ajustando modelos
q <- 1:9 / 10


J1 <- 6
J2 <- 6


burn <- 100000
jump <- 100


guess <- c(-10, -20, rep(0, J1 + 1), rep(0, J2 + 1))
 
models3 <- lapply(q, function(a) {
  aux <- ggcrqr::bayesGG_mm(
    linear_pred_mu = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov, 
    linear_pred_alpha = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
    data = dados,
    q = a,
    d = "ind_obito",
    iter = 1000,
    burn = burn,
    jump = jump,
    guess = guess
  )
  saveRDS(aux, paste("models_mixture_q_", a,".rds", sep = ""))
})







# 
# 
outputs <- function(chains) {
  out <- data.frame(
    names_vars = colnames(chains),
    post_mean = apply(chains, 2, mean),
    sd =  apply(chains, 2, stats::sd),
    emp_lq = apply(chains, 2, stats::quantile, probs = 0.025),
    emp_uq = apply(chains, 2, stats::quantile, probs = 0.925),
    hpd = t(apply(
      chains, 2, TeachingDemos::emp.hpd,
      conf = 0.95
    ))
  )
  rownames(out) <- NULL
  out
}

