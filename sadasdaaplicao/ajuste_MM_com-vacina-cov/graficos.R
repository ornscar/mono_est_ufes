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
library("patchwork")

############ Calling the functions #######################
source("functions_mixture-model.R")
source("diag_conver.R") 
#function to get the convergence graphs

##########################
dados <-  readRDS("dados_modelo.rds")

dados$faixa2 <- as.factor(dados$faixa2)
dados$grupos <- as.factor(dados$grupos)
dados$obesidade <- as.factor(dados$obesidade)
dados$perd_pala <- as.factor(dados$perd_pala)
dados$saturacao <- as.factor(dados$saturacao)
dados$vacina_cov <- as.factor(dados$vacina_cov)
# dados$vacina_cov2 <- as.factor(dados$vacina_cov2)

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

#para diferentes q's
q <- 1:9 / 10

####### todos mu ##############################
m1 <- readRDS("models_mixture_q_0.1.rds")
m2 <- readRDS("models_mixture_q_0.2.rds")
m3 <- readRDS("models_mixture_q_0.3.rds")
m4 <- readRDS("models_mixture_q_0.4.rds")
m5 <- readRDS("models_mixture_q_0.5.rds")
m6 <- readRDS("models_mixture_q_0.6.rds")
m7 <- readRDS("models_mixture_q_0.7.rds")
m8 <- readRDS("models_mixture_q_0.8.rds")
m9 <- readRDS("models_mixture_q_0.9.rds")

models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)

estimates <- lapply(models, function(a) {
  estimates <- outputs(a$post)
  # estimates$tau <- a$tau
  estimates
}) %>% do.call(rbind.data.frame, .)

estimates$tau <- rep(q, each = 16)

e1 <- estimates %>% 
  filter(tau == 0.5)

str(estimates)

d <- estimates %>%
  mutate(names_vars = as.factor(names_vars)) %>%
  filter(names_vars != "alpha" & names_vars != "lambda" &
           names_vars != "gama[1]" & names_vars != "gama[2]" &
           names_vars != "gama[3]" & names_vars != "gama[4]" &
           names_vars != "gama[5]" & names_vars != "gama[6]" &
           names_vars != "gama[7]") %>%
  mutate(nomes_vars = factor(rep(
    c("intercepto",
      "idade >= 30",
      "puérpera",
      "obesidade",
      "perda de paladar",
      "saturação",    
      "vacina"
      ),
    times = length(models)
  ),
  levels =         c("intercepto",
                     "idade >= 30",
                     "puérpera",
                     "obesidade",
                     "perda de paladar",
                     "saturação",    
                     "vacina"
  )))

g <- ggplot2::ggplot(d, aes(x = tau)) +
  facet_wrap( ~ nomes_vars, ncol = 2, scales = "free") +
  geom_line(aes(y = post_mean)) +
  geom_line(aes(y = hpd.1), linetype = 2) +
  geom_line(aes(y = hpd.2), linetype = 2) +
  # geom_line(aes(y = media), colour = 'red',  alpha = 0.5) +
  # geom_line(aes(y = HPD_2.5), colour = 'red', linetype = 2,  alpha = 0.5) +
  # geom_line(aes(y = HPD_97.5), colour = 'red', linetype = 2,  alpha = 0.5) +
  theme(text = element_text(size = 14)) +
  theme_minimal() +
  labs(y = "estimativas", x = "quantis") +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0, 1, 0.1)) +
  theme(
    text = element_text(size = 20),
    panel.background =  element_rect(),
    legend.background = element_rect()
    #axis.ticks.x = element_line(size = 0.5)
  ) +
  geom_hline(yintercept = 0,
             linetype = 3,
             color = "grey50")

g

ggsave(g,
       filename = "coeficientes_MM_todas.png",
       width = 14,
       height = 14)


####### final model ##############################
m1 <- readRDS("models_mixture_final_q_0.1.rds")
m2 <- readRDS("models_mixture_final_q_0.2.rds")
m3 <- readRDS("models_mixture_final_q_0.3.rds")
m4 <- readRDS("models_mixture_final_q_0.4.rds")
m5 <- readRDS("models_mixture_final_q_0.5.rds")
m6 <- readRDS("models_mixture_final_q_0.6.rds")
m7 <- readRDS("models_mixture_final_q_0.7.rds")
m8 <- readRDS("models_mixture_final_q_0.8.rds")
m9 <- readRDS("models_mixture_final_q_0.9.rds")

models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)

estimates <- lapply(models, function(a) {
  estimates <- outputs(a$post)
  # estimates$tau <- a$tau
  estimates
}) %>% do.call(rbind.data.frame, .)

estimates$tau <- rep(q, each = 12)

str(estimates)

e2 <- estimates %>% 
  filter(tau == 0.5)

d <- estimates %>%
  mutate(names_vars = as.factor(names_vars)) %>%
  filter(names_vars != "alpha" & names_vars != "lambda" &
           names_vars != "gama[1]" & names_vars != "gama[2]" &
           names_vars != "gama[3]" & names_vars != "gama[4]" &
           names_vars != "gama[5]" & names_vars != "gama[6]") %>%
  mutate(nomes_vars = factor(rep(
    c("intercepto",
      "puérpera",
      "saturação",    
      "vacina" 
    ),
    times = length(models)
  ),
  levels =         c("intercepto",
                     "puérpera",
                     "saturação",    
                     "vacina"
  )))

g <- ggplot2::ggplot(d, aes(x = tau)) +
  facet_wrap( ~ nomes_vars, ncol = 2, scales = "free") +
  geom_line(aes(y = post_mean)) +
  geom_line(aes(y = hpd.1), linetype = 2) +
  geom_line(aes(y = hpd.2), linetype = 2) +
  # geom_line(aes(y = media), colour = 'red',  alpha = 0.5) +
  # geom_line(aes(y = HPD_2.5), colour = 'red', linetype = 2,  alpha = 0.5) +
  # geom_line(aes(y = HPD_97.5), colour = 'red', linetype = 2,  alpha = 0.5) +
  theme(text = element_text(size = 14)) +
  theme_minimal() +
  labs(y = "estimativas", x = "quantis") +
  scale_x_continuous(limits = c(0.05, 0.95), breaks = seq(0, 1, 0.1)) +
  theme(
    text = element_text(size = 20),
    panel.background =  element_rect(),
    legend.background = element_rect()
    #axis.ticks.x = element_line(size = 0.5)
  ) +
  geom_hline(yintercept = 0,
             linetype = 3,
             color = "grey50")

g

ggsave(g,
       filename = "coeficientes_MM_modelo-final.png",
       width = 14,
       height = 14)


comb_covars <- expand.grid(inter = 1, faixa2 = 0:1, grupos = 0:1,
                           obesidade = 0:1,
                           saturacao = 0:1, vacina_cov = 0:1)
comb_covars  <- as.matrix(comb_covars)

q <- 1:9 / 10
for (j in 1:length(q)){
  models[[j]][["q"]] = q[j]
}

p0_estimates <- lapply(1:dim(comb_covars)[1], function(a) {
  est <- lapply(1:length(models), function(aaa) {
    est_p0 <- exp(comb_covars[a,] %*% t(models[[aaa]]$post[, 7:12]))/(1 + exp(comb_covars[a,] %*% t(models[[aaa]]$post[, 7:12])))
 
    est_p0
  })
   est[[7]]
})


faixa2_post <- c("<30", ">=30")
grupos_post <- c("gestante", "puérpera")
obesidade_post <- c("não", "sim") 
saturacao_post  <- c("não", "sim")
vacina_cov_post <- c("não", "sim")


size_psample <- dim(models[[5]]$post)[1]

data_plot2 <- data.frame(
  p0_estimates = unlist(p0_estimates),
  faixa2 = rep(
    c(rep(faixa2_post[1], size_psample),
    rep(faixa2_post[2], size_psample)), 16
  ),
  grupos = rep(
    c(rep(grupos_post[1], size_psample*2),
      rep(grupos_post[2], size_psample*2)), 8
  ),
  obesidade = rep(
    c(rep(obesidade_post[1], size_psample*4),
      rep(obesidade_post[2], size_psample*4)), 4
  ), 
  saturacao = rep(
    c(rep(saturacao_post[1], size_psample*8),
      rep(saturacao_post[2], size_psample*8)), 2
  ),
  vacina_cov = c(rep(vacina_cov_post[1], size_psample*16),
      rep(vacina_cov_post[2], size_psample*16)
      )
)


##### por grupo ##########################################
data_plot2$grupos <- forcats::fct_relevel(data_plot2$grupos,
                                          "gestante", "puérpera")



g2_faixa <- ggplot(data_plot2) +
  theme_minimal() +
  aes(x = p0_estimates, y = ..density.., colour = faixa2) +
  geom_density(adjust = 1.5) +
  facet_wrap( ~ grupos, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "probabilidade de cura",
       y = "densidade da posteriori - fração de cura",
       colour = "idade") +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 20),
        panel.background =  element_rect(),
        #axis.ticks.x = element_line(size = 0.5)
        )
g2_faixa


g2_obesidade <- ggplot(data_plot2) +
  theme_minimal() +
  aes(x = p0_estimates, y = ..density.., colour = obesidade) +
  geom_density(adjust = 1.5) +
  facet_wrap( ~ grupos, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "probabilidade de cura",
       y = "densidade da posteriori - fração de cura",
       colour = "obesidade") +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 20),
        panel.background =  element_rect(),
        #axis.ticks.x = element_line(size = 0.5)
        )
g2_obesidade


g2_saturacao <- ggplot(data_plot2) +
  theme_minimal() +
  aes(x = p0_estimates, y = ..density.., colour = saturacao) +
  geom_density(adjust = 1.5) +
  facet_wrap( ~ grupos, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "probabilidade de cura",
       y = "densidade da posteriori - fração de cura",
       colour = "saturação") +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 20),
        panel.background =  element_rect(),
        #axis.ticks.x = element_line(size = 0.5)
        )
g2_saturacao


g2_vacina_cov <- ggplot(data_plot2) +
  theme_minimal() +
  aes(x = p0_estimates, y = ..density.., colour = vacina_cov) +
  geom_density(adjust = 1.5) +
  facet_wrap( ~ grupos, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "probabilidade de cura",
       y = "densidade da posteriori - fração de cura",
       colour = "vacina") +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 20),
        panel.background =  element_rect(),
        #axis.ticks.x = element_line(size = 0.5)
        )
g2_vacina_cov


p2 <- (g2_vacina_cov + g2_saturacao)/(g2_faixa + g2_obesidade)
ggsave('graficos_probs-modelo-final_por-grupo.png', p2,
       width = 14,
       height = 14)
 

### Table for CPO ####
cpo <- rep(0, length(q))
for (i in 1:length(q)) {
  cpo[i] <- cpo_mm(Posterior = models[[i]]$post, 
                   linear_pred1 = tempo_hosp_evolucao ~  grupos + saturacao + vacina_cov, 
                   linear_pred2 = tempo_hosp_evolucao ~ faixa2 + grupos + obesidade + saturacao + vacina_cov,
                   dados = dados, d = "ind_obito",
                   q = q[i]) 
}

write.table(cpo, "cpo_mm.txt")


  
 
