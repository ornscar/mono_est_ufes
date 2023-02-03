#### analises 23/01/2023

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
library(foreign)
library(survminer)
library(ggcrqr)

############ Calling the functions #######################
source("functions.R")
source("diag_conver.R") 
#function to get the convergence graphs

##### 


#medias ergódicas média dos parâmetros  
med.erg2 <- function(x)
{
  aux <- rep(0,length(x[,1]))
  for (i in 1:length(x[,1]))                              {
    aux[i] <- mean(x[i,])                                                   
  }
  out <- rep(0,length(aux))
  out[1] <- aux[1]
  for (i in 2:length(aux))
  {
    out[i] <- (out[i-1]*(i-1)+aux[i])/i
  }
  return(out)
}



##########################
dados <-  readRDS("dados_modelo.rds")

dados$faixa2 <- as.factor(dados$faixa2)
dados$grupos <- as.factor(dados$grupos)
dados$obesidade <- as.factor(dados$obesidade)
dados$perd_pala <- as.factor(dados$perd_pala)
dados$saturacao <- as.factor(dados$saturacao)
dados$vacina_cov <- as.factor(dados$vacina_cov)
# dados$vacina_cov2 <- as.factor(dados$vacina_cov2)

###################################################################
############# todas covars in mu ##################################
###################################################################
## faixa2, grupos, obesidade, perd_pala, saturacao, vacina_cov
## 1, 1, 1, 1, 1, 1

#### convergence avaluation ########################
qe <- 0.5

J1 <- 6
J2 <- 0
J3 <- 0

##vamos primeiro encontrar valor de burnin
guess1 <- c(rep(0, J1 + 1), rep(0, J2 + 1), rep(0, J3 + 1))
guess2 <- c(rep(1, J1 + 1), rep(0, J2 + 1), rep(0, J3 + 1))
guess3 <- c(rep(0.7, J1 + 1), rep(0.1, J2 + 1), rep(0.1, J3 + 1))
guess4 <- c(rep(-0.5, J1 + 1), rep(-0.5, J2 + 1), rep(-0.5, J3 + 1))

burn <- 200000
jump <- 1

aj1 <- bayesGG_ergodicas(
  linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
  linear_pred2 = tempo_hosp_evolucao ~ 1,
  linear_pred3 = tempo_hosp_evolucao ~ 1, 
  data = dados,
  q = qe,
  d = "ind_obito",
  iter = 1000,
  burn = burn,
  jump = jump,
  guess = guess1
)

aj2 <-
  bayesGG_ergodicas(
    linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
    linear_pred2 = tempo_hosp_evolucao ~ 1,
    linear_pred3 = tempo_hosp_evolucao ~ 1, 
    data = dados,
    q = qe,
    d = "ind_obito",
    iter = 1000,
    burn = burn,
    jump = jump,
    guess = guess2
  )

aj3 <-
  bayesGG_ergodicas(
    linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
    linear_pred2 = tempo_hosp_evolucao ~ 1,
    linear_pred3 = tempo_hosp_evolucao ~ 1, 
    data = dados,
    q = qe,
    d = "ind_obito",
    iter = 1000,
    burn = burn,
    jump = jump,
    guess = guess3
  )

aj4 <-
  bayesGG_ergodicas(
    linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
    linear_pred2 = tempo_hosp_evolucao ~ 1,
    linear_pred3 = tempo_hosp_evolucao ~ 1, 
    data = dados,
    q = qe,
    d = "ind_obito",
    iter = 1000,
    burn = burn,
    jump = jump,
    guess = guess4
  )

#gráfico de médias ergódicas
png("medias_ergodicas_todas-todas-mu.png")
par(mfrow = c(1, 1))
plot(ts(med.erg2(aj1$post)), type = "l", col = "1")
lines(ts(med.erg2(aj2$post)), type = "l", col = "2")
lines(ts(med.erg2(aj3$post)), type = "l", col = "3")
lines(ts(med.erg2(aj4$post)), type = "l", col = "4")
dev.off()


##### burnin de 100000
burn <- 100000
jump <- 150

guess <- c(rep(0, J1 + 1), rep(0, J2 + 1), rep(0, J3 + 1))

aj <- ggcrqr::bayesGG(tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov, 
                ~ 1, 
                ~ 1,
                data = dados,
                q = qe,
                d = "ind_obito", burn = burn,
                jump = jump,
                guess = guess
)

saveRDS(aj, "ajuste_final_todas-mu.rds")

parm.names <- LaplacesDemon::as.parm.names(list(
  beta = rep(0, (J1 + 1)),
  gama = rep(0, (J2 + 1)), 
  vi = rep(0, (J3 + 1))
))
arquivo <- paste("Ajuste_todas-mu_q", qe, "_burn", burn, "_jump", jump, sep = "")
graphs_diag(aj$post, par.names = parm.names, arquivo)


#para diferentes q's
q <- 1:9 / 10


models3 <- lapply(q, function(a) {
  aux <- ggcrqr::bayesGG(tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov, 
                  ~ 1, 
                  ~ 1,
                  data = dados,
                  q = a,
                  d = "ind_obito", burn = burn,
                  jump = jump,
                  guess = guess
  )
  saveRDS(aux, paste("models_todas-mu_q_", a,".rds", sep = ""))
})

# saveRDS(models3, "models_todas-mu.rds")
# 
# models <- readRDS("models_todas-mu.rds")



#####################################################################
#### com grupo em alpha ###################
#####################################################################

qe <- 0.5

J1 <- 6
J2 <- 1
J3 <- 0

##vamos primeiro encontrar valor de burnin
guess1 <- c(rep(0, J1 + 1), rep(0, J2 + 1), rep(0, J3 + 1))
guess2 <- c(rep(1, J1 + 1), rep(0, J2 + 1), rep(0, J3 + 1))
guess3 <- c(rep(0.7, J1 + 1), rep(0.1, J2 + 1), rep(0.1, J3 + 1))
guess4 <- c(rep(-0.5, J1 + 1), rep(-0.5, J2 + 1), rep(-0.5, J3 + 1))

burn <- 200000
jump <- 1

aj1 <- bayesGG_ergodicas(
  linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
  linear_pred2 = tempo_hosp_evolucao ~ grupos,
  linear_pred3 = tempo_hosp_evolucao ~ 1, 
  data = dados,
  q = qe,
  d = "ind_obito",
  iter = 1000,
  burn = burn,
  jump = jump,
  guess = guess1
)

aj2 <-
  bayesGG_ergodicas(
    linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
    linear_pred2 = tempo_hosp_evolucao ~ grupos,
    linear_pred3 = tempo_hosp_evolucao ~ 1, 
    data = dados,
    q = qe,
    d = "ind_obito",
    iter = 1000,
    burn = burn,
    jump = jump,
    guess = guess2
  )

aj3 <-
  bayesGG_ergodicas(
    linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
    linear_pred2 = tempo_hosp_evolucao ~ grupos,
    linear_pred3 = tempo_hosp_evolucao ~ 1, 
    data = dados,
    q = qe,
    d = "ind_obito",
    iter = 1000,
    burn = burn,
    jump = jump,
    guess = guess3
  )

aj4 <-
  bayesGG_ergodicas(
    linear_pred1 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov,
    linear_pred2 = tempo_hosp_evolucao ~ grupos,
    linear_pred3 = tempo_hosp_evolucao ~ 1, 
    data = dados,
    q = qe,
    d = "ind_obito",
    iter = 1000,
    burn = burn,
    jump = jump,
    guess = guess4
  )

#gráfico de médias ergódicas
png("medias_ergodicas_grupo-alpha.png")
par(mfrow = c(1, 1))
plot(ts(med.erg2(aj1$post)), type = "l", col = "1")
lines(ts(med.erg2(aj2$post)), type = "l", col = "2")
lines(ts(med.erg2(aj3$post)), type = "l", col = "3")
lines(ts(med.erg2(aj4$post)), type = "l", col = "4")
dev.off()

##### burnin de 100000
burn <- 100000
jump <- 150

guess <- c(rep(0, J1 + 1), rep(0, J2 + 1), rep(0, J3 + 1))

aj <- ggcrqr::bayesGG(tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov, 
                      ~ grupos, 
                      ~ 1,
                      data = dados,
                      q = qe,
                      d = "ind_obito", burn = burn,
                      jump = jump,
                      guess = guess
)


saveRDS(aj, "ajuste_final_grupo-alpha.rds")

parm.names <- LaplacesDemon::as.parm.names(list(
  beta = rep(0, (J1 + 1)),
  gama = rep(0, (J2 + 1)),
  vi = rep(0, (J3 + 1))
))
arquivo <- paste("Ajuste_grupo-alpha_q", qe, "_burn", burn, "_jump", jump, sep = "")
graphs_diag(aj$post, par.names = parm.names, arquivo)



#para diferentes q's
q <- 1:9 / 10


models3 <- lapply(q, function(a) {
  aux <- ggcrqr::bayesGG(tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + perd_pala + saturacao + vacina_cov, 
                         ~ grupos, 
                         ~ 1,
                         data = dados,
                         q = a,
                         d = "ind_obito", burn = burn,
                         jump = jump,
                         guess = guess
  )
  saveRDS(aux, paste("models_grupo-alpha_q_", a,".rds", sep = ""))
})

 




 