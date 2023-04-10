library(lmtest)
set.seed(7)

dane <- read.csv('MB_2_dane.csv', sep = ';')
dane[,1:3] <- log(dane[,1:3])

dane0 <- read.csv('MB_2_dane0.csv', sep = ';')
dane0[,1:3] <- log(dane0[,1:3])

index = sort(sample(nrow(dane), 5))
dane_testowe <- dane[index,]
dane <- dane[-index,]

# ------------------------ MNK - wcześniejsze badania ------------------
model_mnk0 <- lm(cena ~ ., data = dane0)
summary(model_mnk0)

# ------------------------ MNK -----------------------------------------
model_mnk <- lm(cena ~ ., data = dane)
summary(model_mnk)
shapiro.test(model_mnk$residuals)
bptest(model_mnk)

# ------------------------ a priori # ----------------------------------
dane0_X <- as.matrix(cbind(rep(1, nrow(dane0)),  dane0[,2:4]))

alfa0 <- model_mnk0$df.residual
beta0 <- as.matrix(model_mnk0$coefficients)
delta0 <- sum(model_mnk0$residuals^2)
sigma0 <- solve(t(dane0_X) %*% dane0_X)

# ------------------------ a posteriori --------------------------------
dane_X <- as.matrix(cbind(rep(1, nrow(dane)),  dane[,2:4]))
dane_Y <- as.matrix(dane[,1])

sigma1 <- solve(t(dane_X)%*%dane_X + solve(sigma0))
beta1 <- sigma1%*%(t(dane_X)%*%dane_Y + solve(sigma0)%*%beta0)
alfa1 <- alfa0 + nrow(dane)
delta1 <- delta0 + t(dane_Y)%*%dane_Y - t(beta1)%*%solve(sigma1)%*% beta1 + t(beta0)%*%solve(sigma0)%*%beta0

# ------------------------ rozkłady brzegowe ---------------------------
library(invgamma)
library(MCMCpack)
library(metRology)
n <- 1000000 # liczba losowań z danego rozkładu

# rozkłady brzegowe rozkładu a priori
density(invgamma::rinvgamma(n, shape = alfa0/2, rate = delta0/2)) # parametr sigma
plot(density(invgamma::rinvgamma(n, shape = alfa0/2, rate = delta0/2)))

rt.scaled(n, mean = beta0[1], sd = as.vector((delta0/alfa0) * sigma0[1,1]), df = alfa0) # stała
rt.scaled(n, mean = beta0[2], sd = as.vector((delta0/alfa0) * sigma0[2,2]), df = alfa0) # B1
rt.scaled(n, mean = beta0[3], sd = as.vector((delta0/alfa0) * sigma0[3,3]), df = alfa0) # B2
rt.scaled(n, mean = beta0[4], sd = as.vector((delta0/alfa0) * sigma0[4,4]), df = alfa0) # B3

# rozkłady brzegowe rozkładu a posteriori
density(rinvgamma(n, shape = alfa1/2, rate = delta1/2)) # parametr sigma
plot(density(rinvgamma(n, shape = alfa1/2, rate = delta1/2)))

a_posteriori <- as.data.frame(matrix(NA, nrow = n, ncol = ncol(dane)))
for(i in 1:(ncol(a_posteriori))){
  a_posteriori[,i] <- rt.scaled(n, mean = beta1[i], sd = as.vector((delta1/alfa1) * sigma1[i,i]), df = alfa1)
}

# wykresy (przykład):
library(ggplot2)
library(plyr)
df1 <- data.frame(
  Rozkład=factor(rep(c("a priori", "a posteriori"), each=n)),
  wartość=c(rt.scaled(n, mean = beta0[1], sd = as.vector((delta0/alfa0) * sigma0[1,1]), df = alfa0),
                 rt.scaled(n, mean = beta1[1], sd = as.vector((delta1/alfa1) * sigma1[1,1]), df = alfa1)))
mu1 <- ddply(df1, "Rozkład", summarise, grp.mean=mean(wartość))
p1 <- ggplot(df1, aes(x=wartość, fill=Rozkład)) +
  geom_density(alpha=0.4) + 
  geom_vline(data=mu1, aes(xintercept=grp.mean, color=Rozkład), linetype="dashed") +
  theme_bw() +
  ylab('') +
  ggtitle(expression("Rozkłady parametru" ~ beta[0]))
p1

# ------------------------ bayesowskie estymatory parametrów modelu ---------------------------
estymator_bayesowski <- colMeans(a_posteriori)
apply(a_posteriori, 2, mean)
apply(a_posteriori, 2, sd)

# HPDI
library(bayestestR)
apply(a_posteriori, 2, hdi)

# ------------------------ rozkłady predykcyjne dla zbioru testowego ---------------------------
library(Metrics)
x_tau <- as.matrix(cbind(rep(1, nrow(dane_testowe)), dane_testowe[,2:4]))
  
a_posteriori_prognoza <- as.data.frame(matrix(NA, nrow = n, ncol = nrow(dane_testowe)))
p1 <- diag(1, nrow(dane_testowe)) +  x_tau%*%sigma1%*%t(x_tau) 
for(i in 1:(nrow(dane_testowe))){
  a_posteriori_prognoza[,i] <- rt.scaled(n, 
                                         mean = (x_tau %*% as.matrix(beta1))[i], 
                                         sd = as.vector((delta1/alfa1) * p1[i,i]), 
                                         df = alfa1)
}

colMeans(a_posteriori_prognoza)
predict(model_mnk, dane_testowe[,2:4], interval = "confidence")

# hdi
apply(a_posteriori_prognoza, 2, hdi)

# błędy prognozy 
mae(dane_testowe$cena, c(predict(model_mnk, dane_testowe[,2:4])))
mae(dane_testowe$cena, colMeans(a_posteriori_prognoza))

mse(dane_testowe$cena, c(predict(model_mnk, dane_testowe[,2:4])))
mse(dane_testowe$cena, colMeans(a_posteriori_prognoza))

rmse(dane_testowe$cena, c(predict(model_mnk, dane_testowe[,2:4])))
rmse(dane_testowe$cena, colMeans(a_posteriori_prognoza))

# ------------------------ rozkłady predykcyjne dla zbioru testowego - II wersja ---------------------
library(fCopulae)

a_posteriori_prognoza_2 <- rmvst(n,
      dim = nrow(dane_testowe),
      mu = c(x_tau %*% as.matrix(beta1)),
      Omega = round(as.vector((delta1/alfa1)) * (diag(1,nrow(dane_testowe)) + x_tau%*%sigma1%*%t(x_tau)), 10),
      df = alfa1
)

colMeans(a_posteriori_prognoza_2)

# ------------------------ istotność jednej zmiennej X1 (pamiec) ---------------------------
# model a priori dla modelu z restrykcjami R1
model_mnk0_R1 <- lm(cena ~ PPI + lcd, data = dane0)
summary(model_mnk0_R1)

dane0_X_R1 <- as.matrix(cbind(rep(1, nrow(dane0)),  dane0[,3:4]))
alfa0_R1 <- model_mnk0_R1$df.residual
beta0_R1 <- as.matrix(model_mnk0_R1$coefficients)
delta0_R1 <- sum(model_mnk0_R1$residuals^2)
sigma0_R1 <- solve(t(dane0_X_R1) %*% dane0_X_R1)

model_mnk_R1 <- lm(cena ~ PPI + lcd, data = dane)
summary(model_mnk_R1)

dane_X_R1 <- as.matrix(cbind(rep(1, nrow(dane)),  dane[,3:4]))
dane_Y_R1 <- as.matrix(dane[,1])
sigma1_R1 <- solve(t(dane_X_R1)%*%dane_X_R1 + solve(sigma0_R1))
beta1_R1 <- sigma1_R1%*%(t(dane_X_R1)%*%dane_Y_R1 + solve(sigma0_R1)%*%beta0_R1)
alfa1_R1 <- alfa0_R1 + nrow(dane)
delta1_R1 <- delta0_R1 + t(dane_Y_R1)%*%dane_Y_R1 - t(beta1_R1)%*%solve(sigma1_R1)%*% beta1_R1 + 
  t(beta0_R1)%*%solve(sigma0_R1)%*%beta0_R1

statystyka_R1 = sqrt((det(sigma1)*det(sigma0_R1)) / (det(sigma1_R1)*det(sigma0))) * 
  (delta0^(alfa0/2) / delta1^(alfa1/2)) / (delta0_R1^(alfa0_R1/2) / delta1_R1^(alfa1_R1/2)) *
  (gamma(alfa1/2) / gamma(alfa0/2)) / (gamma(alfa1_R1/2) / gamma(alfa0_R1/2))
     
log(statystyka_R1, 10) # model bez restrykcji jest dużo lepszy 

# ------------------------ istotność dwóch zmiennych X1 i X2(pamiec i PPI) ---------------------------
# model a priori dla modelu z restrykcjami R2
model_mnk0_R2 <- lm(cena ~ lcd, data = dane0)
summary(model_mnk0_R2)

dane0_X_R2 <- as.matrix(cbind(rep(1, nrow(dane0)),  dane0[,4]))
alfa0_R2 <- model_mnk0_R2$df.residual
beta0_R2 <- as.matrix(model_mnk0_R2$coefficients)
delta0_R2 <- sum(model_mnk0_R2$residuals^2)
sigma0_R2 <- solve(t(dane0_X_R2) %*% dane0_X_R2)

model_mnk_R2 <- lm(cena ~ lcd, data = dane)
summary(model_mnk_R2)

dane_X_R2 <- as.matrix(cbind(rep(1, nrow(dane)),  dane[,4]))
dane_Y_R2 <- as.matrix(dane[,1])
sigma1_R2 <- solve(t(dane_X_R2)%*%dane_X_R2 + solve(sigma0_R2))
beta1_R2 <- sigma1_R2%*%(t(dane_X_R2)%*%dane_Y_R2 + solve(sigma0_R2)%*%beta0_R2)
alfa1_R2 <- alfa0_R2 + nrow(dane)
delta1_R2 <- delta0_R2 + t(dane_Y_R2)%*%dane_Y_R2 - t(beta1_R2)%*%solve(sigma1_R2)%*% beta1_R2 + 
  t(beta0_R2)%*%solve(sigma0_R2)%*%beta0_R2

statystyka_R2 = sqrt((det(sigma1)*det(sigma0_R2)) / (det(sigma1_R2)*det(sigma0))) * 
  (delta0^(alfa0/2) / delta1^(alfa1/2)) / (delta0_R2^(alfa0_R2/2) / delta1_R2^(alfa1_R2/2)) *
  (gamma(alfa1/2) / gamma(alfa0/2)) / (gamma(alfa1_R2/2) / gamma(alfa0_R2/2))

log(statystyka_R2, 10) # model bez restrykcji jest dużo lepszy 

library(car)
summary(model_mnk)
summary(model_mnk_R2)
linearHypothesis(model_mnk, c("pamiec=0", "PPI=0"))

RSSS = sum(model_mnk_R2$residuals^2)
URSS = sum(model_mnk$residuals^2)
F_statystyka = ((RSSS-URSS)/URSS)*((nrow(dane)-3-1)/2)
F_statystyka
# odrzucamy H0, restrykcje są zbędne
