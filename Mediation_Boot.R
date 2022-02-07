
##################
#PARTE 2
##################
# https://stats.stackexchange.com/questions/184654/simulating-data-for-mediation-model

#test for GIT

# Let say you have the mediation model : x→m→y and you know the correlations: 
#  a=ρ(x,m),
#  b=ρ(m,y|x) and 
#  c=ρ(x,y). 
#  You then know the indirect effect ab=a∗b and the direct effect c′=c−ab.

#var(X + Y) = var(X) + var(Y) + 2*cov(X,Y)

# https://stats.stackexchange.com/questions/72468/simulating-data-to-fit-a-mediation-model

# 1) Generating data STEP BY STEP

n  <-  100
a <-  .6
b <-  .4
c <-  .3
ab <-  a*b # regra multiplicativa a*b = c-cp
cp <-  c-ab #cp = c'

# (1.1) generating x
set.seed(123)
x <- rnorm(n) # generantig x with normal distribution
hist(x)
summary(x)
library(psych)
describe(x)

# (1.2) generating m
em <- sqrt(1-a^2) # variância residual de m - takes a (standardized) beta, squares it to calculate the proportion of variance explained (R2), then subtracts that from 1 to get the remaining variance, then takes the square root to get back to a beta coefficient.

m = a*x + em*rnorm(n) # Variância aleatória com média 0 e SD 1
hist(m)

# (1.3) generating y
# Y follows a normal distribution, Y ∼ N (0; 1)
ey = 1-(cp^2 + b^2 + 2*a*cp*b) # var(x+m)= var(x) + var(m)+ 2cov(xm)

# the error of y (ey) is the residual variance of var(x)+var(m)+2cov(x,m).
# see wright's tracing rules

y = cp*x + b*m + ey*rnorm(n)

hist(y)



df1 <-  data.frame(x,m,y)
library(psych)
describe(df1)

round(cor(df1),2)


## 1.2) Using builtin function from caron & Volois, 2018

GenerateMediationData <- function(n = 10000, a = .60, b =.40, c =.30, mean.x = 0,
                                  sd.x = 1, mean.m = 0, sd.m = 1, mean.y = 0, sd.y = 1) {
  # a is a_xm
  # b is b_my|x
  # c is c_xy
  # mean.x, sd.x , mean.m, sd.m, mean.y and sd.y will create unstandardized data according to the specified
  # means and standard deviations
  if(missing(a) | missing(b) | missing(c)){
    stop("One or more arguments are missing")
  }
  ab <- a*b
  cp <- c-ab # cp = c’ = c_xy|m
  ey <- 1-(cp^2 +b^2 + 2*a*cp*b)
  if ((ey < 0) | (ey > 1)){print("WARNING : Sum of square of coefficients is too
high to generate standardized data")}
  # Generate data
  x <- rnorm(n, mean = 0, sd = 1)
  em <- sqrt(1-a^2)
  m <- a*x + em*rnorm(n, mean = 0, sd = 1)
  ey2 <- sqrt(ey)
  y <- cp*x + b*m + ey2*rnorm(n, mean = 0, sd = 1)
  x <- x * sd.x + mean.x
  m <- m * sd.m + mean.m
  y <- y * sd.y + mean.y
  data <- as.data.frame(cbind(x, m, y))
  return(data)
}

set.seed(123)
df2 <- GenerateMediationData(n = 100, a = .6, b =.4, c =.3)
describe(df2)
cor(df2)
cor(df1)
#################################
#################################
#Simulação de dados
# Mediation models:
# (1) Y=b11+b12X+e1 
# (2) med=b21+b22X+e2
# (3) Y=b31+b32X+b32med+e3
set.seed(1234)
df3 <- GenerateMediationData(n = 200, a = .6, b =.4, c =.3)
describe(df3)
round(cor(df3),2)

hist(df3$x)
hist(df3$y)
hist(df3$m)

# Check the relationships between x, med, and y
model1 <- lm(y ~ x, data=df3)
summary(model1)

# Check the relationship between x and med
model2 <- lm(m ~ x, data=df3)
summary(model2)


model3 <- lm(y ~ x + m, data=df3)
summary(model3)

########################
#Monte carlo (using simglm package)
########################
install.packages("simglm")
library(simglm)

n  <-  200
a <-  .6
b <-  .4
c <-  .3
ab <-  a*b # regra multiplicativa a*b = c-cp
cp <-  c-ab #cp = c'

###################
#simulating x
###################
sim_arguments <- list(
  formula = y ~ 1 + x,
  fixed = list(x = list(var_type = 'continuous', 
                         mean = 0, sd = 1)),
  sample_size = 200
)

dfsim_x <- simulate_fixed(data = NULL, sim_arguments)

data.frame(dfsim_x$x)
mean(dfsim_x$x)
var(dfsim_x$x)
hist(dfsim_x$x)
###################
# simulate error for m
###################
sim_arguments <- list(
  error = list(variance = sqrt(1-a^2)),
  sample_size = 200
)

dfsim_em <- simulate_error(data = NULL, sim_arguments)

mean(dfsim_em$error) # close to 0
var(dfsim_em$error) # 

###################
# simulating m
###################
sim_arguments <- list(
  formula = m ~ 1 + x,
  fixed = list(x = list(var_type = 'continuous', mean = 0, sd = 1)),
  error = list(variance = sqrt(1-a^2)),
  sample_size = 200,
  reg_weights = c(0,a)
)

dfsim_m <- simulate_fixed(data = NULL, sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments)

hist(dfsim_m$m)
names(dfsim_m)
summary(lm(m~x, data=dfsim_m))

#Fitting the model
set.seed(321) 
sim_arguments <- list(
  formula =  m ~ 1 + x,
  fixed = list(x = list(var_type = 'continuous', mean = 0, sd = 1)),
  error = list(variance = sqrt(1-a^2)),
  sample_size = 200,
  reg_weights = c(0,a)
)
#These were estimated using the lm function based on the same formula defined in sim_arguments. 

simulate_fixed(data = NULL, sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments) %>% 
  model_fit(sim_arguments) %>%
  extract_coefficients()

data_m <- simulate_fixed(data = NULL, sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments) 

names(data_m)

###################
# Simulating y
###################
# Y follows a normal distribution, Y ∼ N (0; 1)
# ey = 1-(cp^2 + b^2 + 2*a*cp*b) # var(x+m)= var(x) + var(m)+ 2cov(xm)

# the error of y (ey) is the residual variance of var(x)+var(m)+2cov(x,m).
# see wright's tracing rules
# y = cp*x + b*m + ey*rnorm(n)
#
##Criando novos dados (ver opção usando dados de x e m prévios)
# sim_arguments <- list(
#   formula =  y ~ 1 + x + m,
#   fixed = list(x = list(var_type = 'continuous', mean = 0, sd = 1),
#                m = list(var_type = 'continuous', mean = 0, sd = 1)),
#   error = list(variance = 1-(cp^2 + b^2 + 2*a*cp*b)),
#   sample_size = 200,
#   reg_weights = c(0,cp, b)
# )
# simulate_fixed(data = NULL, sim_arguments) %>%
#   simulate_error(sim_arguments) %>%
#   generate_response(sim_arguments) %>% 
#   model_fit(sim_arguments) %>%
#   extract_coefficients()
# 
# data_y <- simulate_fixed(data = NULL, sim_arguments) %>%
#   simulate_error(sim_arguments) %>%
#   generate_response(sim_arguments) 

# using m data

sim_arguments <- list(
  formula =  y ~ 1 + x + m,
  error = list(variance = 1-(cp^2 + b^2 + 2*a*cp*b)),
  sample_size = 200,
  reg_weights = c(cp, b)
)

simulate_error(data = data_m[,c(2,7)], sim_arguments) %>%
  generate_response(sim_arguments) %>% 
  model_fit(sim_arguments) %>%
  extract_coefficients()


df2 <- simulate_error(data = data_m[,c(2,7)], sim_arguments) %>%
  generate_response(sim_arguments)

cor(df2[, c("x", "m", "y")])
# simulate_error(data = data_m[,c(2,7)], sim_arguments) %>%
#   generate_response(sim_arguments)
# 
# names(data_m)
# data_test1 <- simulate_error(data = data_m[,c(2,7)], list(
#   formula =  y ~ 1 + x + m,
#   error = list(variance = 1-(cp^2 + b^2 + 2*a*cp*b)),
#   sample_size = 200)
#   ) 
# names(data_test1)

# data_test1$y <- 0 + cp*data_test1$x + b*data_test1$m + data_test1$error
# data_test2 <- generate_response(data_test1, sim_args=list(
#   formula =  y ~ 1 + x + m + error,
#   reg_weights = c(0, cp, b)))
#   
# summary(lm(y ~ 1 + x + m, data=data_test2))
# 
# names(data_test)
# head(data_test)
# generate_response(data_m, sim_args=sim_arguments) %>% 
#   model_fit(sim_arguments) 
# 
# ?generate_response
# 
# generate_response(sim_arguments) %>% 
#   model_fit(sim_arguments) %>%
#   extract_coefficients()
# 
# simulate_error(data_m, sim_arguments) %>%
# generate_response(sim_arguments) %>% 
#     model_fit(sim_arguments) %>%
#     extract_coefficients()
# 
# 
# head(data_y)
# 
# head(data_m)

# Determining power using monte carlo approach (replication) for x -> m 
set.seed(321) 

sim_arguments_1 <- list(
  formula =  m ~ 1 + x,
  fixed = list(x = list(var_type = 'continuous', mean = 0, sd = 1)),
  error = list(variance = sqrt(1-a^2)),
  sample_size = 200,
  reg_weights = c(0,a),
  model_fit = list(formula = m ~ 1 + x,
                   model_function = 'lm'),
  reg_weights_model = c(0,a),
  replications = 5000,
  extract_coefficients = TRUE
)

dfsim_m_rep <- replicate_simulation(sim_arguments_1) 

head(dfsim_m_rep)

simcoef <- function(listsim,
                    parameter = "x",
                    value = "estimate") {
  list_coef <- lapply(listsim, as.data.frame)
  coefs <- as.numeric()
  for (i in 1:length(list_coef)) {
    coefs[i] <- list_coef[[i]][list_coef[[i]]$term == parameter, value]
  }
  return(coefs)
}

mean(simcoef(dfsim_m_rep))
hist(simcoef(dfsim_m_rep))

compute_statistics(dfsim_m_rep, sim_arguments_1)


# Determining power using monte carlo approach (replication) for x + m -> y 
set.seed(321) 

sim_arguments_2 <- list(
  formula =  y ~ 1 + x + m,
  fixed = list(x = list(var_type = 'continuous', mean = 0, sd = 1),
               m = list(var_type = 'continuous', mean = 0, sd = 1)),
  error = list(variance = 1-(cp^2 + b^2 + 2*a*cp*b)),
  sample_size = 200,
  reg_weights = c(0,cp, b),
  model_fit = list(formula = y ~ 1 + x + m,
                   model_function = 'lm'),
  reg_weights_model = c(0,cp, b),
  replications = 5,
  extract_coefficients = TRUE
)

dfsim_y_rep <- replicate_simulation(sim_arguments_2) 

head(dfsim_y_rep)




mean(simcoef(dfsim_y_rep, parameter="x"))
hist(simcoef(dfsim_y_rep, parameter="x"))

mean(simcoef(dfsim_y_rep, parameter="m"))
hist(simcoef(dfsim_y_rep, parameter="m"))

compute_statistics(dfsim_y_rep, sim_arguments_2)




#################################
#################################
# testing inderect effect
#################################
#################################

# https://m-clark.github.io/posts/2019-03-12-mediation-models/#data


install.packages('bda') 
library(bda) 
head(df1)

########
#Sobel test
########

names(df1)

library(bda)
library(MBESS)
mediation.test
mediation.test(mv=df1$m,iv=df1$x,dv=df1$y)
?mediate


sobeltest <- function(data, x, m, y, alpha=.05){
#Passo 1: extrair coeficientes
model1 <- lm(formula = m ~ x, data=data)
summary(model1)
a <- model1$coefficient[2]

model2 <- lm(formula = y ~ x + m, data=data)
summary(model2)
b <- model2$coefficient[3]

# Passo 2: Extrair desvio padrão (SE) das estimativas
SEa <- coef(summary(model1))[2, 2]
SEb <- coef(summary(model2))[3, 2]

# Passo 3: devio padrão do efeito indireto e estatística do teste
SE <- sqrt(a^2*SEb^2 + b^2*SEa^2)
z <- a*b/SE

# Passo 4: Encontrar probabilidade de observar o valor do teste tomando por referência um distribuição z. Ou então verificar se o valor absoluto da estatistaca do teste é superior ao valor crítico para um determinado alpha
p <- pnorm(-abs(z)) * 2
# and not p <- 1-pnorm(z) as in paper

sig <- qnorm(1-alpha/2) < abs(z)
return(list(teststat=z, pvalue=p, sig=sig))
}

test <- sobeltest(x, m, y, data=df1)

#############
# Implementamos dois tipos de bootstrap o Percentile Bootstrap Method (bootstrap percentílico) e o Bias corrected e accelarated bootstrap
#############

hist(rnorm(30)*rnorm(30))


# Para os dois tipos de bootstrap é necessário criar uma função geradora das estimativas do para ab da mediação. Usamos cálculo matricial para encontrar o parametro ab - poderia ser usado a função lm mas assim o processo torna-se computacionamente mais eficiente

# ver apontamentos de estatística aplicada
# (mat <- matrix[,1:3]) # Construir matriz X (matrix do modelo)
# (mat_t <- t(mat)) # matriz transposta Xt
# (mat_XtX <- mat_t %*% mat) # Calcular Xt*X
# (mat_XtX_inv <- solve(mat_XtX)) # inverter matriz
# (betas <- mat_inv %*% mat_t %*% y)

d <- as.matrix(cbind(int=rep(1, length(df$x)), x=df$x, m=df$m, y=df$y))
b <- solve(t(d[,1:3])%*%d[,1:3])%*%t(d[,1:3])%*%d[,4]
a <- solve(t(d[,1:2])%*%d[,1:2])%*%t(d[,1:2])%*%d[,3]
ab <- a[2,1]*b[3,1]

# com base neste procedimento, criar função para estimar efeitos indiretos no contexto do bootstrap
stat_function <- function(data){
  d <- as.matrix(cbind(int=rep(1, dim(data)[1]), x=data[,"x"], m=data[,"m"], y=data[,"y"]))
  b <- solve(t(d[,1:3])%*%d[,1:3])%*%t(d[,1:3])%*%d[,4]
  a <- solve(t(d[,1:2])%*%d[,1:2])%*%t(d[,1:2])%*%d[,3]
  ab <- a[2,1]*b[3,1] 
  return(ab)
}
stat_function(df1)

#############
# Percentile Bootstrap Method (bootstrap percentílico)
#############
boot_perc <- function(data, R=5000, alpha = .05) {
  data <- as.matrix(data)
  est <- stat_function(data)
  n <- dim(data)[1]
  N <- 1:n
  res <- rep(NA, R)
  for (i in 1:R) {
    #Select a bootstrap sample
    id <- sample(n, replace = TRUE)
    #Estimate index
    res[i] <- stat_function(data[id, ])
  }
  limits <- quantile(res, c(alpha / 2, 1 - alpha / 2))
  CI <- c(limits[[1]], limits[[2]])
  sig <- 0 < prod(sign(CI)) # serve para ver se o intervalo de confiânça contem ou não zero
  return(list(
    estimate = est,
    bootdist=res,
    BootConfInt = limits,
    sig = sig
  ))
}

testboot1 <- boot_perc(df1)
str(testboot1)

#Plot bootstrap distribution
hist(testboot1[2], xlab='theta', main='Bootstrap distribution')
abline(v=testboot1[1], lwd=2, col='blue')
abline(v=testboot1[[3]][1], lwd=2, col='orange')
abline(v=testboot1[[3]][2], lwd=2, col='orange')


###########
#2. Bias corrected e accelarated bootstrap
###########
# https://stats.stackexchange.com/questions/437477/calculate-accelerated-bootstrap-interval-in-r

# Essentially, resampling methods allow us to estimate the bias of an estimator (see the Jackknife estimator for a simple example). The accelerated bootstrap also accounts for (and adjusts for) the possible skewness of the Bootstrap distribution.

# Although the method of the previous section is a straightforward and natural way to obtain endpoints for a confidence interval, there are several alternatives which have been shown to perform better in a variety of settings. The Accelerated Bootstrap is one such method.
# 
# The endpoints to the CI in this approach are found by considering the function
# 
# 
# Este bootstrap implica um constant de acelaração que é calculada através do método jackknife


#res <- BCa.boot(data = d, stat = indirect, R = R)
# indirect <- ab
# stat <- indirect

boot_bca <- function(data, R=5000, alpha = .05) {
  data <- as.matrix(data)
  n <- dim(data)[1]
  est <- stat_function(data) # cálculo do parametro amostral
  # o primero loop faz o bootstrap
  res <- rep(NA, R)
  for (i in 1:R) {
    #Select a bootstrap sample
    id <- sample(n, replace = TRUE)
    #Estimate index
    res[i] <- stat_function(data[id,])
  }
  #Calcular as constantes
  z0 <- qnorm(mean(res <= rep(est, R)))
  zc <- qnorm(c(alpha / 2, 1 - alpha / 2)) # Quantis pretendidos
  
  # este segundo faz o jackknife para estimar o parametro a (acelaração)
  I <- rep(NA, n)
  for (i in 1:n) {
    #Remove ith data point
    datanew <- data[-i, ]
    #Estimate theta - distribuição jackkife para todos os possíveis n-1
    theta_jack <- stat_function(datanew)
    I[i] <- (n - 1) * (est - theta_jack)
  }
  
  #calculo da constante de acelaração
  a <- 1 / 6 * (sum(I ^ 3) / sum(I ^ 2) ^ 1.5)
  
  # Calculo dos quantis ajustados
  adj.alpha <- pnorm(z0 + (z0 + zc) / (1 - a * (z0 + zc)))
  limits <- quantile(res, adj.alpha)
  CI <- c(limits[[1]], limits[[2]])
  sig <- 0 < prod(sign(CI)) # serve para ver se o intervalo de confiânça contem ou não zero
  return(list(
    estimate = est,
    bootdist = res,
    BootConfInt = limits,
    sig = sig
  ))
}


(testboot2 <- boot_bca(df1, R=5000))

#Plot bootstrap distribution
hist(testboot1[[2]], xlab='theta', main='Bootstrap distribution')
abline(v=testboot1[[1]], lwd=2, col='blue')
abline(v=testboot1[[3]][1], lwd=2, col='orange')
abline(v=testboot1[[3]][2], lwd=2, col='orange')
abline(v=testboot2[[3]][1], lwd=2, col='darkgreen')
abline(v=testboot2[[3]][2], lwd=2, col='darkgreen')


#################################
#################################
# Análise da potência do test
#################################
#################################
# Power refers to the probability to find a significant result when the null hypothesis is false (there is an indirect effect). Failure to find a significant result is a type II error. 


# The purpose of power analysis is to simulate an experiment with known and non-null population parameters, check whether the result is significant or not, and redo the above a tremendous amount of times. 

# There are two main components in the script: the generation of data (Listing 1) and the indirect effect test (Listings 2 to 4). The outcome of the function is the power of the mediation test given a sample size n


PowerMediation <- function(MediationTest, a = .25, b = .6, c = .0, n = 40, R =
                             100, alpha = 0.05){
  # Warning : This function can be excessively slow with high replication values and high sample sizes,
  # especially with bootstrap
  # MediationTest = SobelTest.R or BaronKenny.R or BootTest.R or any function returning an output labelled sig indicating if the result is significant (TRUE or FALSE)
  TestSIG <- rep(NA, R)
  for(j in 1:R){
    data <- GenerateMediationData(n=n, a=a, b=b, c=c)
    RES <- MediationTest(data=data, alpha=alpha)
    TestSIG[j] <- RES$sig
  }
  Power <- round(sum(TestSIG)/R,3)
  return(list(Results=TestSIG, Power=Power))
}

##############
#Cenário 1
##############
# n = 20, 30, 40, 50, 60, 70, 80, 90, 100, 150
# ba= .1
# cp = 0

a <-  .335
b <-  .30
c <-  .15
(ab <-  a*b)# regra multiplicativa a*b = c-cp
(cp <-  c-ab) #cp = c'

#SOBEL
(powersobel_20_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(20), R=100))
(powersobel_30_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(30), R=100))
(powersobel_40_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(40), R=100))
(powersobel_50_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(50), R=100))
(powersobel_60_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(60), R=100))
(powersobel_70_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(70), R=100))
(powersobel_80_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(80), R=100))
(powersobel_90_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(90), R=100))
(powersobel_100_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(100), R=100))
(powersobel_150_c1 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(150), R=100))

#boot_perc
(powerbot1_20_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(20), R=100))
(powerbot1_30_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(30), R=100))
(powerbot1_40_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(40), R=100))
(powerbot1_50_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(50), R=100))
(powerbot1_60_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(60), R=100))
(powerbot1_70_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(70), R=100))
(powerbot1_80_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(80), R=100))
(powerbot1_90_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(90), R=100))
(powerbot1_100_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(100), R=100))
(powerbot1_150_c1 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(150), R=100))

# boot_bca
(powerbot2_20_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(20), R=100))
(powerbot2_30_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(30), R=100))
(powerbot2_40_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(40), R=100))
(powerbot2_50_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(50), R=100))
(powerbot2_60_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(60), R=100))
(powerbot2_70_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(70), R=100))
(powerbot2_80_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(80), R=100))
(powerbot2_90_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(90), R=100))
(powerbot2_100_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(100), R=100))
(powerbot2_150_c1 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(150), R=100))


##############
#Cenário 2 - efeito indireto = .2
##############
# n = 20, 30, 40, 50, 60, 70, 80, 90, 100, 150
# ba= .5
# cp = 0

a <-  .5
b <-  .4
c <-  .25
(ab <-  a*b) # regra multiplicativa a*b = c-cp
(cp <-  c-ab) #cp = c'

#SOBEL
(powersobel_20_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(20), R=100))
(powersobel_30_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(30), R=100))
(powersobel_40_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(40), R=100))
(powersobel_50_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(50), R=100))
(powersobel_60_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(60), R=100))
(powersobel_70_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(70), R=100))
(powersobel_80_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(80), R=100))
(powersobel_90_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(90), R=100))
(powersobel_100_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(100), R=100))
(powersobel_150_c2 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(150), R=100))

#boot_perc
(powerbot1_20_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(20), R=100))
(powerbot1_30_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(30), R=100))
(powerbot1_40_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(40), R=100))
(powerbot1_50_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(50), R=100))
(powerbot1_60_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(60), R=100))
(powerbot1_70_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(70), R=100))
(powerbot1_80_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(80), R=100))
(powerbot1_90_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(90), R=100))
(powerbot1_100_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(100), R=100))
(powerbot1_150_c2 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(150), R=100))

# boot_bca
(powerbot2_20_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(20), R=100))
(powerbot2_30_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(30), R=100))
(powerbot2_40_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(40), R=100))
(powerbot2_50_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(50), R=100))
(powerbot2_60_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(60), R=100))
(powerbot2_70_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(70), R=100))
(powerbot2_80_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(80), R=100))
(powerbot2_90_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(90), R=100))
(powerbot2_100_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(100), R=100))
(powerbot2_150_c2 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(150), R=100))


##############
#Cenário 3 - efeito indireto = .3
##############
# n = 20, 30, 40, 50, 60, 70, 80, 90, 100, 150
# ba= .5
# cp = 0
a <-  .6
b <-  .5
c <-  .35
(ab <-  a*b)# regra multiplicativa a*b = c-cp
(cp <-  c-ab) #cp = c'


#SOBEL
(powersobel_20_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(20), R=100))
(powersobel_30_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(30), R=100))
(powersobel_40_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(40), R=100))
(powersobel_50_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(50), R=100))
(powersobel_60_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(60), R=100))
(powersobel_70_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(70), R=100))
(powersobel_80_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(80), R=100))
(powersobel_90_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(90), R=100))
(powersobel_100_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(100), R=100))
(powersobel_150_c3 <- PowerMediation(sobeltest, a = a, b = b, c = c, n = c(150), R=100))

#boot_perc
(powerbot1_20_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(20), R=100))
(powerbot1_30_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(30), R=100))
(powerbot1_40_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(40), R=100))
(powerbot1_50_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(50), R=100))
(powerbot1_60_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(60), R=100))
(powerbot1_70_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(70), R=100))
(powerbot1_80_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(80), R=100))
(powerbot1_90_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(90), R=100))
(powerbot1_100_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(100), R=100))
(powerbot1_150_c3 <- PowerMediation(boot_perc, a = a, b = b, c = c, n = c(150), R=100))

# boot_bca
(powerbot2_20_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(20), R=100))
(powerbot2_30_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(30), R=100))
(powerbot2_40_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(40), R=100))
(powerbot2_50_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(50), R=100))
(powerbot2_60_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(60), R=100))
(powerbot2_70_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(70), R=100))
(powerbot2_80_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(80), R=100))
(powerbot2_90_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(90), R=100))
(powerbot2_100_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(100), R=100))
(powerbot2_150_c3 <- PowerMediation(boot_bca, a = a, b = b, c = c, n = c(150), R=100))


#table and plot for power


porwer_df <- data.frame(
  condition = c(rep("Cenário 1: ba=.1", 10), rep("Cenário 2: ba=.2", 10), rep("Cenário 2: ba=.3", 10)),
  n = rep(c(seq(20, 100, 10), 150), 3),
  power_sobel =  c(
    c(powersobel_20_c1$Power, powersobel_30_c1$Power, powersobel_40_c1$Power, 
      powersobel_50_c1$Power, powersobel_60_c1$Power, powersobel_70_c1$Power,
      powersobel_80_c1$Power, powersobel_90_c1$Power, powersobel_100_c1$Power,
      powersobel_150_c1$Power),
    c(powersobel_20_c2$Power, powersobel_30_c2$Power, powersobel_40_c2$Power,
      powersobel_50_c2$Power, powersobel_60_c2$Power, powersobel_70_c2$Power,
      powersobel_80_c2$Power, powersobel_90_c2$Power, powersobel_100_c2$Power,
      powersobel_150_c2$Power),
    c(powersobel_20_c3$Power, powersobel_30_c3$Power, powersobel_40_c3$Power,
      powersobel_50_c3$Power, powersobel_60_c3$Power, powersobel_70_c3$Power,
      powersobel_80_c3$Power, powersobel_90_c3$Power, powersobel_100_c3$Power,
      powersobel_150_c3$Power)),
  power_bootperc =  c(
    c(powerbot1_20_c1$Power, powerbot1_30_c1$Power, powerbot1_40_c1$Power,
      powerbot1_50_c1$Power, powerbot1_60_c1$Power, powerbot1_70_c1$Power, 
      powerbot1_80_c1$Power, powerbot1_90_c1$Power, powerbot1_100_c1$Power,
      powerbot1_150_c1$Power),
    c(powerbot1_20_c2$Power, powerbot1_30_c2$Power, powerbot1_40_c2$Power,
      powerbot1_50_c2$Power, powerbot1_60_c2$Power, powerbot1_70_c2$Power,
      powerbot1_80_c2$Power, powerbot1_90_c2$Power, powerbot1_100_c2$Power,
      powerbot1_150_c2$Power), 
    c(powerbot1_20_c3$Power, powerbot1_30_c3$Power, powerbot1_40_c3$Power, 
      powerbot1_50_c3$Power, powerbot1_60_c3$Power, powerbot1_70_c3$Power, 
      powerbot1_80_c3$Power, powerbot1_90_c3$Power, powerbot1_100_c3$Power, 
      powerbot1_150_c3$Power)), 
  power_bootcba =  c(
    c(powerbot2_20_c1$Power, powerbot2_30_c1$Power, powerbot2_40_c1$Power,
      powerbot2_50_c1$Power, powerbot2_60_c1$Power, powerbot2_70_c1$Power, 
      powerbot2_80_c1$Power, powerbot2_90_c1$Power, powerbot2_100_c1$Power,
      powerbot2_150_c1$Power),
    c(powerbot2_20_c2$Power, powerbot2_30_c2$Power, powerbot2_40_c2$Power, 
      powerbot2_50_c2$Power, powerbot2_60_c2$Power, powerbot2_70_c2$Power, 
      powerbot2_80_c2$Power, powerbot2_90_c2$Power, powerbot2_100_c2$Power, 
      powerbot2_150_c2$Power), 
    c(powerbot2_20_c3$Power, powerbot2_30_c3$Power, powerbot2_40_c3$Power, 
      powerbot2_50_c3$Power, powerbot2_60_c3$Power, powerbot2_70_c3$Power, 
      powerbot2_80_c3$Power, powerbot2_90_c3$Power, powerbot2_100_c3$Power, 
      powerbot2_150_c3$Power))
)


power_sobel <- cbind(test=rep("sobel",30), porwer_df[, c(1,2,3)])
names(power_sobel)[4] <- "power"

power_bootperc <- cbind(test=rep("bootperc",30), porwer_df[, c(1,2,4)])
names(power_bootperc)[4] <- "power"

power_bootcba <- cbind(test=rep("bootbca",30), porwer_df[, c(1,2,5)])
names(power_bootcba)[4] <- "power"


power_df <- rbind(power_sobel, power_bootperc, power_bootcba)
head(power_df)

save(power_df, file="power.Rdata")
load(file="power.Rdata")
# smoothingSpline_sobel  <-  smooth.spline(power_df$n, power_df$sobel, spar=0.5)
# smoothingSpline_bot1  <-  smooth.spline(power_df$n, power_df$bootperc, spar=0.5)
# smoothingSpline_bot2  <-  smooth.spline(power_df$n, power_df$bootBCa, spar=0.5)
# plot(power_df$n, power_df$sobel, type="l", col="white")
# lines(smoothingSpline_sobel, col="darkred", lwd=2)
# lines(smoothingSpline_bot1, col="darkgreen", lwd=2)
# lines(smoothingSpline_bot2, col="darkblue", lwd=2)



# qplot(x,y, geom='smooth', span =0.5)

library(ggplot2)

ggplot(power_df, aes(x=n, y=power, color=test)) + 
  geom_point() +
  geom_line() + 
  facet_grid(cols = vars(condition)) + 
  ggtitle("Potencia dos testes para diferentes n tamanho do efeito indireto ") + 
  geom_hline(yintercept=.8, linetype="dashed", 
             color = "darkgrey", size=.8)


ggplot(power_df) + 
  geom_point(aes(x=n, y=power, color=test), size=1)+
  geom_smooth(ymax=1,aes(x=n, y=power, color=test), method = "glm", formula = y ~ poly(x, 2), se = FALSE) + 
  facet_grid(cols = vars(condition)) +  
  ggtitle("Potencia dos testes para diferentes n tamanho do efeito indireto \n (linha suavizada)") + 
  geom_hline(yintercept=.8, linetype="dashed", 
             color = "darkgrey", size=.8)




################################################
#Caron & Valois, 2018
################################################


GenerateMediationData <- function(n = 10000, a = .60, b =.40, c =.30, mean.x = 0,
                                  sd.x = 1, mean.m = 0, sd.m = 1, mean.y = 0, sd.y = 1) {
  # a is a_xm
  # b is b_my|x
  # c is c_xy
  # mean.x, sd.x , mean.m, sd.m, mean.y and sd.y will create unstandardized data according to the specified
  # means and standard deviations
  if(missing(a) | missing(b) | missing(c)){
    stop("One or more arguments are missing")
  }
  ab <- a*b
  cp <- c-ab # cp = c’ = c_xy|m
  ey <- 1-(cp^2 +b^2 + 2*a*cp*b)
  if ((ey < 0) | (ey > 1)){print("WARNING : Sum of square of coefficients is too
high to generate standardized data")}
  # Generate data
  x <- rnorm(n, mean = 0, sd = 1)
  em <- sqrt(1-a^2)
  m <- a*x + em*rnorm(n, mean = 0, sd = 1)
  ey2 <- sqrt(ey)
  y <- cp*x + b*m + ey2*rnorm(n, mean = 0, sd = 1)
  x <- x * sd.x + mean.x
  m <- m * sd.m + mean.m
  y <- y * sd.y + mean.y
  data <- as.data.frame(cbind(x, m, y))
  return(data)
}

df_test <- GenerateMediationData(100, a=a,b=b,c=c)

BaronKenny <- function(x, m, y, data, alpha = 0.05) {
  # x is the column name of the predictor in data
  # m is the column name of the mediator in data
  # y is the column name of the dependent variable in data
  # data is a data.frame or a matrix that contain columns with the names of x, y and m.
  # If x, m or y are missing , data [,1] will be used for x, data [,2] will be
  # used for m and data [,3] will be used for y.
  if(missing(data)) {
    stop("There’s no data")
  }
  if(is.data.frame(data) != TRUE & is.matrix(data)) {
    d <- as.data.frame(data)
  } else if (is.data.frame(data) != TRUE & is.matrix(data) != TRUE) {
    stop('"data" should be a matrix or a data.frame')
  }
  if(missing(x)){x <- data[,1]} else if (is.numeric(x) == TRUE) {x <- data[,x]}
  else {x <- data[,match(x, table = colnames(data))]}
  if(missing(m)){m <- data[,1]} else if (is.numeric(m) == TRUE) {m <- data[,m]}
  else {m <- data[,match(m, table = colnames(data))]}
  if(missing(y)){y <- data[,1]} else if (is.numeric(y) == TRUE) {y <- data[,y]}
  else {y <- data[,match(y, table = colnames(data))]}
  Sig <- FALSE
  out <- 'No Mediation'
  #regression 1
  step1 <- lm(formula = y ~ x)
  pC <- summary(step1)$coefficients[2,4]
  if (pC <= alpha){
    #regression 2
    step2 <- lm(formula = m ~ x)
    pA <- summary(step2)$coefficients[2,4]
    if (pA <= alpha){
      #regression 3
      step3 <- lm(formula = y ~ x + m)
      pB <- summary(step3)$coefficients[3,4]
      Sig <- (pB <= alpha)
      if (Sig){
        if (summary(step3)$coefficients[2,4] <= alpha){
          out <- 'This is a partial mediation'} else {
            out <- 'This is a complete mediation'}
      }
    }
  }
  return(list(sig=Sig, conclusion=out))
}

BaronKenny(x="x", y="y", m="m",data=df_test)

SobelTest <- function(x, y, m, data, alpha = 0.05) {
  # x is the column name or number of the predictor in data
  # m is the column name or number of the mediator in data
  # y is the column name or number of the dependent variable in data
  # data is a data.frame or a matrix that contain columns with the names of x, y and m.
  # If x, m or y are missing , data [,1] will be used for x, data [,2] will be
  # used for m and data [,3] will be used for y.
  if(missing(data)) {
    stop("There’s no data")
  }
  if(is.data.frame(data) != TRUE & is.matrix(data)) {
    d <- as.data.frame(data)
  } else if (is.data.frame(data) != TRUE & is.matrix(data) != TRUE) {
    stop('"data" should be a matrix or a data.frame')
  }
  if(missing(x)){x <- data[,1]} else if (is.numeric(x) == TRUE) {x <- data[,x]}
  else {x <- data[,match(x, table = colnames(data))]}
  if(missing(m)){m <- data[,1]} else if (is.numeric(m) == TRUE) {m <- data[,m]}
  else {m <- data[,match(m, table = colnames(data))]}
  if(missing(y)){y <- data[,1]} else if (is.numeric(y) == TRUE) {y <- data[,y]}
  else {y <- data[,match(y, table = colnames(data))]}
  step1 <- lm(formula = m ~ x)
  step2 <- lm(formula = y ~ x + m)
  a <- step1$coefficient[2]
  SEa <- coef(summary(step1))[2, 2]
  b <- step2$coefficient[3]
  SEb <- coef(summary(step2))[3, 2]
  SE <- sqrt(a^2*SEb^2 + b^2*SEa^2)
  z <- a*b/SE
  p <- 1-pnorm(z)
  sig <- qnorm(1-alpha/2) < abs(z)
  return(list(z = z, p = p, sig = sig))
}

names(df_test)
SobelTest(x="x", y="y", m="m",data=df_test)

BootTest <- function(x, y, m, data, alpha = 0.05, R = 5000) {
  # Warning : This function can be excessively slow with high replication values and high sample sizes
  # x is the column name or number of the predictor in data
  # m is the column name or number of the mediator in data
  # y is the column name or number of the dependent variable in data
  # data is a data.frame or a matrix that contain columns with the names pf x, y and m.
  # If x, m or y are missing , data [,1] will be used for x, data [,2] will be
  # used for m and data [,3] will be used for y.
  # R is the number of replication
  if(missing(data)) {
    stop("There’s no data")
  }
  if(is.data.frame(data) != TRUE & is.matrix(data)) {
    d <- as.data.frame(data)
  } else if (is.data.frame(data) != TRUE & is.matrix(data) != TRUE) {
    stop('"data" should be a matrix or a data.frame')
  }
  if(missing(x)){x <- data[,1]} else if (is.numeric(x) == TRUE) {x <- data[,x]}
  else {x <- data[,match(x, table = colnames(data))]}
  if(missing(m)){m <- data[,1]} else if (is.numeric(m) == TRUE) {m <- data[,m]}
  else {m <- data[,match(m, table = colnames(data))]}
  if(missing(y)){y <- data[,1]} else if (is.numeric(y) == TRUE) {y <- data[,y]}
  else {y <- data[,match(y, table = colnames(data))]}
  d <- as.matrix(cbind(x, m, y))
  # Compute the indirect effect for the Bca.boot function
  indirect <- function(data, indice) {
    d <- data[indice,]
    b <- solve(t(d[,1:2])%*%d[,1:2])%*%t(d[,1:2])%*%d[,3]
    a <- solve(t(d[,1])%*%d[,1])%*%t(d[,1])%*%d[,2]
    ab <- a*b[2]
    return(ab)
  }
  res <- BCa.boot(data = d, stat = indirect, R = R)
  sig <- 0 < prod(sign(res$BCaCI))
  return(list(ab = round(res$estimate,3), CI = round(res$BCaCI, 3), sig = sig))
}

# Bootstraping function with the bias corrected and accelerated boostrap interval (BCa)
BCa.boot = function(data, stat, R = 5000, alpha=0.05){
  # data is the data to bootstrap # stat is the function to bootstrap
  # R is the number of replication # alpha is significance threshold
  data <- as.matrix(data)
  n <- dim(data)[1]
  N <- 1:n
  res <- rep(0,R)
  zj <- rep(0,n)
  est <- stat(data,indice=N)
  M <- max(R,n)
  for (i in 1:M){
    if(i<=R){
      id <- sample(n, replace = TRUE)
      res[i] <- stat(data=data,indice=id)
    }
    if(i<=n){
      J <- N[1:(n-1)]
      zj[i] <- stat(data[-i,],J)
    }
  }
  z0 <- qnorm(sum(res < rep(est,R))/R)
  zc <- qnorm(c(alpha/2,1-alpha/2))
  L <- mean(zj)-zj
  a <- sum(L^3)/(6*sum(L^2)^1.5)
  adj.alpha <- pnorm(z0 + (z0 + zc) / (1 - a * (z0 + zc)))
  limits <- quantile(res,adj.alpha)
  CI <- c(limits[[1]], limits[[2]])
  return(list(estimate = est, BCa=limits, BCaCI = CI))
}

df_test <- GenerateMediationData(100, a=a,b=b,c=c)
BootTest(x="x", y="y", m="m",data=df_test)

PowerMediation <- function(MediationTest, a = .25, b = .6, c = .0, n = 40, R =
                             5000, alpha = 0.05){
  # Warning : This function can be excessively slow with high replication values and high sample sizes ,
  # especially with bootstrap
  # MediationTest = SobelTest.R or BaronKenny.R or BootTest.R or any function returning an output
  # labelled sig indicating if the result is significant (TRUE or FALSE)
  SIG <- 0
  for(j in 1:R){
    data <- GenerateMediationData(n=n, a=a, b=b, c=c)
    RES <- MediationTest(data=data, x="x", y="y", m="m", alpha=alpha)
    SIG <- RES$sig + SIG
  }
  Power <- round(SIG/R,3)
  return(list(Power=Power))
}

a <-  .3
b <-  .4
c <-  .2
a*b
PowerMediation(BaronKenny, a = a, b = b, c = c, R=50, n=40)
PowerMediation(SobelTest, a = a, b = b, c = c, R=50, n=40)
PowerMediation(BootTest, a = a, b = b, c = c, R=50, n=40)

###############################################
###############################################
###############################################
# other boootstrap approaches
install.packages("robmed")
library("robmed")
BootTest(x="x", m="m", y="y", data=df1)

standard <- test_mediation(test, 
                           x = "x",
                           y = "y",
                           m = "m",
                           robust = FALSE)
summary(standard)
robust <- test_mediation(test, 
                           x = "x",
                           y = "y",
                           m = "m",
                           robust = TRUE)
summary(robust)




##############################
##############################
##############################
# Adding skew to a normal distribution
install.packages("fGarch")
library(fGarch)
library(tidyverse)
N <- 10000
x <- rnbinom(N, 10, .5)
dsnorm(x, mean = 0, sd = 1, xi = 1.5, log = FALSE)
# psnorm(x, mean = 0, sd = 1, xi = 1.5) %>% plot()
# qsnorm(x, mean = 0, sd = 1, xi = 1.5) %>% plot()
# rsnorm(x, mean = 0, sd = 1, xi = 1.5) %>% plot()

hist(rsnorm(N, mean = 0, sd = 1, xi = 1.5))

?dsnorm

x <- rsnorm(10000, mean = 0, sd = 1, xi = 1.5)
##############
##############
##############
##############
##############
##############
##############
##############

