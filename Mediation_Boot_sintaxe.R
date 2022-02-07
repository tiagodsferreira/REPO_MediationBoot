#################################
#################################
#PARTE 2: simulação de dados e monte carlo
#################################
#################################
n  <-  100
a <-  .6
b <-  .4
c <-  .3
ab <-  a*b # regra multiplicativa a*b = c-cp
cp <-  c-ab # cp = c'

# (1.1) gerar x
set.seed(123)
x <- rnorm(n)
hist(x)
summary(x)
library(psych)
describe(x)

# (1.2) gerar m
em <- sqrt(1-a^2) # variância residual de m 
m = a*x + em*rnorm(n) # Variância aleatória com média 0 e SD 1
hist(m)

# (1.3) gerar y
ey = 1-(cp^2 + b^2 + 2*a*cp*b) # var(x+m)= var(x) + var(m)+ 2cov(xm)
y = cp*x + b*m + ey*rnorm(n) # erro de y (var(x)+var(m)+2cov(x,m)
hist(y)

df1 <-  data.frame(x,m,y)
library(psych)
describe(df1)
round(cor(df1),2)

# gerar dados com caron & Volois, 2018
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

set.seed(1234)
df3 <- GenerateMediationData(n = 200, a = .6, b =.4, c =.3)
describe(df3)
round(cor(df3),2)

hist(df3$x)
hist(df3$y)
hist(df3$m)

model1 <- lm(y ~ x, data=df3)
summary(model1)

model2 <- lm(m ~ x, data=df3)
summary(model2)

model3 <- lm(y ~ x + m, data=df3)
summary(model3)


#Monte carlo (simglm package)
install.packages("simglm")
library(simglm)

n  <-  200
a <-  .6
b <-  .4
c <-  .3
ab <-  a*b # regra multiplicativa a*b = c-cp
cp <-  c-ab #cp = c'

#simular x
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

# simulate error for m
sim_arguments <- list(
  error = list(variance = sqrt(1-a^2)),
  sample_size = 200
)

dfsim_em <- simulate_error(data = NULL, sim_arguments)

mean(dfsim_em$error) # close to 0
var(dfsim_em$error) # 

# simular m
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

#ajustar o modelo usando lm
set.seed(321) 
sim_arguments <- list(
  formula =  m ~ 1 + x,
  fixed = list(x = list(var_type = 'continuous', mean = 0, sd = 1)),
  error = list(variance = sqrt(1-a^2)),
  sample_size = 200,
  reg_weights = c(0,a)
)

simulate_fixed(data = NULL, sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments) %>% 
  model_fit(sim_arguments) %>%
  extract_coefficients()

data_m <- simulate_fixed(data = NULL, sim_arguments) %>%
  simulate_error(sim_arguments) %>%
  generate_response(sim_arguments) 

names(data_m)

# Simular y
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


# Determinar potência com monte carlo 
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
# PARTE 4: Teste de efeitos indiretos
#################################
#################################
head(df1)

########
#Sobel test
########
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
# Análise da potência do test
#################################

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

