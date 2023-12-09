## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(1)
data1 <- rnorm(100)
samples <- summary(data1)
population <- c(-Inf, qnorm(0.25), qnorm(0.5), 0, qnorm(0.75), +Inf)
table1 <- cbind(samples, population)
knitr::kable(table1)

## -----------------------------------------------------------------------------
x <- seq(-3, 3, 0.1)
#plot(x, dexp(x), ty="l", main="指数分布", xlab="x", ylab="f(x)")
hist(data1, freq = FALSE, main = 'Histogram of samples from N(0,1)', xlim = c(-3,3), xlab = '')
lines(x, dnorm(x), col="red")

## -----------------------------------------------------------------------------
d <- iris
knitr::kable(d[1:6,])

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
hist(d[1:50, "Sepal.Length"], freq = FALSE, main = 'Histogram of setosa', xlim = c(4,6), xlab = '')
lines(density(d[1:50, "Sepal.Length"]), col = "red")
hist(d[51:100, "Sepal.Length"], freq = FALSE, main = 'Histogram of versicolor', xlim = c(4,8), xlab = '')
lines(density(d[51:100, "Sepal.Length"]), col = "red")
hist(d[101:150, "Sepal.Length"], freq = FALSE, main = 'Histogram of virginica', xlim = c(4,9), xlab = '')
lines(density(d[101:150, "Sepal.Length"]), col = "red")

## -----------------------------------------------------------------------------
my.sample <- function(x, n, p = 0){
#  cat(p,'\n')
  if(length(p) == 1){
    cp <- seq(1/length(x), 1, 1/length(x))
  }else{
    cp <- cumsum(p)
  }
#  cat(cp,'\n')
  U <- runif(n)
  r <- x[findInterval(U,cp)+1]
  return(r)
}

## -----------------------------------------------------------------------------
set.seed(1)
example1.sample <- my.sample(x = 1:3, n = 2000)
example1.ct <- as.vector(table(example1.sample))
example1.ct/sum(example1.ct)/(1/3)

## -----------------------------------------------------------------------------
set.seed(1)
example2.sample <- my.sample(x = 1:4, n = 5000, p = seq(0.1,0.4,0.1))
example2.ct <- as.vector(table(example2.sample))
example2.ct/sum(example2.ct)/(seq(0.1,0.4,0.1))

## -----------------------------------------------------------------------------
set.seed(1)
example3.sample <- my.sample(x = c(1,2,10), n = 1000, p = c(0.1,0.7,0.2))
example3.ct <- as.vector(table(example3.sample))
example3.ct/sum(example3.ct)/(c(0.1,0.7,0.2))

## -----------------------------------------------------------------------------
# the inverse function of distribution function F
F.inverse <- function(x){
  if(0 <= x && x <= 0.5){
    y <- log(2*x)
  }else{
    y <- -log(2*(1-x))
  }
  return(y)
}

## -----------------------------------------------------------------------------
# inverse transform method
set.seed(1)
n <- 1000
u <- runif(n)
x.sample.1 <- sapply(u, F.inverse)

## -----------------------------------------------------------------------------
x.1 <- seq(-10,10,0.1)
hist(x.sample.1, prob = TRUE, xlim = c(-10,10), ylim = c(0,0.5), xlab = "x", main =  "Histogram of x")
lines(x.1, exp(-abs(x.1))/2)

## -----------------------------------------------------------------------------
# preparation
a <- 3
b <- 2
n <- 1000
accept <- 0
x.sample.2 <- numeric(n)

# acceptance-rejection algorithm
set.seed(1)
while(accept < n){
  u <- runif(1)
  x <- runif(1) # random variable from g(.)
  if(x^(a-1)*(1-x)^(b-1) > u){
    # we accept x
    accept <- accept + 1
    x.sample.2[accept] <- x
  }
}

## -----------------------------------------------------------------------------
x.2 <- seq(0,1,0.02)
hist(x.sample.2, prob = TRUE, xlim = c(0,1), ylim = c(0,2), xlab = "x", main = "Histogram of x")
lines(x.2, x.2^(a-1)*(1-x.2)^(b-1)*gamma(a+b)/(gamma(a)*gamma(b)))

## -----------------------------------------------------------------------------
fe <- function(n){
  x <- numeric(n)
  u1 <- runif(n, -1, 1)
  u2 <- runif(n, -1, 1)
  u3 <- runif(n, -1, 1)
  for(i in 1:n){
    if((abs(u3[i]) >= abs(u2[i])) && (abs(u3[i]) >= abs(u1[i]))){
      x[i] <- u2[i]
    }else x[i] <- u3[i]
  }
  return(x)
}

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1000
x.sample.3 <- fe(n)

x.3 <- seq(-1,1,0.05)
hist(x.sample.3, prob = TRUE, xlim = c(-1,1), ylim = c(0, 0.8), xlab = "x", main = "Histogram of x")
lines(x.3, (1-x.3^2)^2*3/4)

## -----------------------------------------------------------------------------
pihat <- function(rho){
  d <- 1
  l <- rho
  n <- 1e6
  X <- runif(n,0,d/2)
  Y <- runif(n,0,pi/2)
  pi.hat <- 2*l/d/mean(l/2*sin(Y)>X)
  return(pi.hat)
}

## -----------------------------------------------------------------------------
K <- 100
rho <- c(0.1, 0.5, 1)
pihat.sample <- matrix(0, nrow = length(rho), ncol = K)

# obtain pihat
set.seed(1)
pihat.sample[1,] <- sapply(rep(rho[1], times = K), pihat)
pihat.sample[2,] <- sapply(rep(rho[2], times = K), pihat)
pihat.sample[3,] <- sapply(rep(rho[3], times = K), pihat)

## -----------------------------------------------------------------------------
# calculate the standard deviation
pihat.var <- apply(pihat.sample, 1, var)

## -----------------------------------------------------------------------------
names(pihat.var) <- c("rho=0.1", "rho=0.5", "rho=1")
#knitr::kable(variance)
pander::pander(pihat.var)

## -----------------------------------------------------------------------------
-exp(2)+3*exp(1)-1
(-3*exp(2)+10*exp(1)-5)/4

## -----------------------------------------------------------------------------
reduction.theoretical <- 2*(exp(2)-3*exp(1)+1)/((exp(1)-1)*(3-exp(1)))
cat("The theoretical percent reduction is", reduction.theoretical, "." ,'\n')

## -----------------------------------------------------------------------------
MC_s <- function(n){
  return(mean(exp(runif(n))))
}
MC_a <- function(n){
  u <- runif(n/2)
  v <- 1-u
  return(mean(exp(u)+exp(v))/2)
}

## -----------------------------------------------------------------------------
n <- 10000
K <- 1000

set.seed(1)
theta_s <- sapply(rep(n,times = K), MC_s)
theta_a <- sapply(rep(n,times = K), MC_a)

## -----------------------------------------------------------------------------
var_s <- var(theta_s)
var_a <- var(theta_a)
reduction.empirical <- (var_s-var_a)/var_s
cat("The empirical percent reduction is",reduction.empirical, "." ,'\n')

## -----------------------------------------------------------------------------
reduction <- c(reduction.theoretical, reduction.empirical)
names(reduction) <- c("theoretical", "empirical")
pander::pander(reduction)

## -----------------------------------------------------------------------------
g <- function(x) x^2 * exp(-x^2/2)/sqrt(2 * pi)*(x>1)
f1 <- function(x) dnorm(x, mean = 0, sd = 1)
f2 <- function(x) dgamma(x, shape = 3, rate = 1)

## -----------------------------------------------------------------------------
x <- seq(1, 8, 0.01)
plot(x, g(x), type = "l", ylim = c(0, 0.4), ylab = "", main = "Figure of g(x), f1(x), f2(x)")
lines(x, f1(x), lty = 2)
lines(x, f2(x), lty = 3)
legend("topright", inset = 0.02, legend = c(expression(g(x)==x^2*e^{-x^2/2}/sqrt(2*pi)), expression(f[1](x)==e^{-x^2/2}/sqrt(2*pi)), expression(f[2](x)==x^2*e^{-x}/2)), lty = 1:3)

## -----------------------------------------------------------------------------
x <- seq(1,5,0.01)
plot(x, g(x)/f1(x), type = "l", lty = 2, ylim = c(0,3), ylab = "", main = "Figure of g(x)/f1(x) and g(x)/f2(x)")
lines(x, g(x)/f2(x), lty = 3)
legend("topright", inset = 0.02, legend = c(expression(f[1](x)==e^{-x^2/2}/sqrt(2*pi)), expression(f[2](x)==x^2*e^{-x}/2)), lty = 2:3)

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- sd.theta.hat <- numeric(2)

# the importance function is f1
set.seed(1)
x <- rnorm(m)
theta.hat[1] <- mean(g(x)/f1(x))
sd.theta.hat[1] <- sd(g(x)/f1(x))

# the importance function is f2
set.seed(1)
x <- rgamma(m,shape = 3,rate = 1)
theta.hat[2] <- mean(g(x)/f2(x))
sd.theta.hat[2] <- sd(g(x)/f2(x))

# result reporting
rbind(theta.hat, sd.theta.hat)

## -----------------------------------------------------------------------------
g <- function(x) exp(-x)/(1+x^2)*(x>0)*(x<1)
f <- function(x) exp(-x)/(1-exp(-1))*(x>0)*(x<1)
M <- 10000
theta.hat <- sd.theta.hat <- numeric(2)

## -----------------------------------------------------------------------------
# data generation and analysis
set.seed(1)
u <- runif(m)
x <- -log(1-u*(1-exp(-1)))
fg <- g(x)/f(x)
theta.hat[1] <- mean(fg)
sd.theta.hat[1] <-sd(fg)

# result reporting
c(theta.hat[1], sd.theta.hat[1])

## -----------------------------------------------------------------------------
k <- 5 # number of subintervals
r <- M/k # replicates per stratum
theta.hat_s <- var_s <- numeric(k)

# data generation and analysis
set.seed(1)
for (j in 1:k) {
  u <- runif(r,(j-1)/5,j/5)
  x <- -log(1-(1-exp(-1))*u)
  fg <- g(x)/k/f(x)
  theta.hat_s[j] <- mean(fg)
  var_s[j] <- var(fg)
}
theta.hat[2] <- sum(theta.hat_s)
sd.theta.hat[2] <- sqrt(sum(var_s))

# result reporting
c(theta.hat[2], sd.theta.hat[2])

## -----------------------------------------------------------------------------
rbind(theta.hat, sd.theta.hat)

## -----------------------------------------------------------------------------
mc <- function(n){
  t0 <- qt(c(0.025, 0.975), df = n - 1)
  x <-  rchisq(n, df = 2)
  ci <- mean(x) + t0 * sd(x)/sqrt(n)
  return(ci)
}

## -----------------------------------------------------------------------------
n <- 20
set.seed(1)
CI <- replicate(10000, expr = mc(n))
LCL <- CI[1, ]
UCL <- CI[2, ]
c(mean(LCL), mean(UCL))
mean(LCL < 2 & UCL > 2)

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05

set.seed(1)
UCL <- replicate(1000, expr = {
  x <- rchisq(n, df = 2)
  (n-1) * var(x) / qchisq(alpha, df = n-1)
} )
mean(UCL)
mean(UCL>4)

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05
m <- 10000 #number of replicates
p <- matrix(0, nrow = 3, ncol = m) #storage for p-values
p.hat <- numeric(3)
se.p.hat <- numeric(3)

## -----------------------------------------------------------------------------
mu0 <- 1

set.seed(1)
for (j in 1:m) {
  x <- rchisq(n, df = 1)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[1,j] <- ttest$p.value
}
p.hat[1] <- mean(p[1,] < alpha)
se.p.hat[1] <- sqrt(p.hat[1] * (1 - p.hat[1]) / m)
print(c(p.hat[1], se.p.hat[1]))

## -----------------------------------------------------------------------------
mu0 <- 1

set.seed(1)
for (j in 1:m) {
  x <- runif(n, min = 0, max = 2)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[2,j] <- ttest$p.value
}
p.hat[2] <- mean(p[2,] < alpha)
se.p.hat[2] <- sqrt(p.hat[2] * (1 - p.hat[2]) / m)
print(c(p.hat[2], se.p.hat[2]))

## -----------------------------------------------------------------------------
mu0 <- 1

set.seed(1)
for (j in 1:m) {
  x <- rexp(n, rate = 1)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[3,j] <- ttest$p.value
}
p.hat[3] <- mean(p[3,] < alpha)
se.p.hat[3] <- sqrt(p.hat[3] * (1 - p.hat[3]) / m)
print(c(p.hat[3], se.p.hat[3]))

## -----------------------------------------------------------------------------
rbind(p.hat, se.p.hat)

## -----------------------------------------------------------------------------
x <- seq(-4, 4, 0.01)
plot(x, dnorm(x, mean = 0, sd = 1), type = "l", ylim = c(0, 2), ylab = "", col = 1)
lines(x, dchisq(x, df = 1), lty = 2, col = 2)
lines(x, dunif(x, min = 0, max = 2), lty = 3, col = 3)
lines(x, dexp(x, rate = 1), lty = 4, col = 4)
legend("topleft", inset = 0.02, legend = c("N(0,1)", "Chisq(1)", "Uniform(0,2)", "Exponential(1)"), lty = 1:4, col = 1:4)

## -----------------------------------------------------------------------------
adjust <- function(m, alpha){
  m0 <- m*0.95
  m1 <- m - m0
  
  # generate p-values
  p0 <- runif(m0)
  p1 <- rbeta(m1, 0.1, 1)
  p <- c(p0, p1)
  
  # adjust p-values
  p.adj_bonf = p.adjust(p, method = 'bonferroni')
  p.adj_BH = p.adjust(p, method = 'BH')
  
  # calculate FWER, FDR and TPR
  FWER_bonf <- 1 - prod(p.adj_bonf[1:m0] >= alpha)
  FWER_BH <- 1 - prod(p.adj_BH[1:m0] >= alpha)
  FDR_bonf <- sum(p.adj_bonf[1:m0] < alpha) / sum(p.adj_bonf < alpha)
  FDR_BH <- sum(p.adj_BH[1:m0] < alpha) / sum(p.adj_BH < alpha)
  TPR_bonf <- sum(p.adj_bonf[m0:m] < alpha) / m1
  TPR_BH <- sum(p.adj_BH[m0:m] < alpha) / m1
  
  return(c(FWER_bonf, FWER_BH, FDR_bonf, FDR_BH, TPR_bonf, TPR_BH))
}

## -----------------------------------------------------------------------------
m <- 1000
M <- 1000
alpha <- 0.1
rate <- matrix(0, nrow = M, ncol = 6)

set.seed(1)
for(i in 1:M){
  rate[i,] <- adjust(m, alpha)
}
estimate <- matrix(apply(rate, 2, mean), nrow = 2, ncol = 3)

## -----------------------------------------------------------------------------
estimate <- round(estimate, 3)
colnames(estimate) <- c("FWER", "FDR", "TPR")
rownames(estimate) <- c("Bonf", "B-H")
knitr::kable(estimate)

## -----------------------------------------------------------------------------
lambda <- 2 #true value
n <- c(5, 10, 20) #sample size
B <- 1000 #bootstrap replicates
m <- 1000 #simulations
lambda.star <- numeric(B)
bias.estimate <- matrix(0, nrow = m, ncol = 3)
se.estimate <- matrix(0, nrow = m, ncol = 3)

## -----------------------------------------------------------------------------
set.seed(1)
for(i in 1:3){
  for(j in 1:m){
    x <- rexp(n[i], rate = 2)
    lambda.hat <- 1/mean(x)
    for(b in 1:B){
      xstar <- sample(x, replace=TRUE)
      lambda.star[b] <- 1/mean(xstar)
    }
    bias.estimate[j,i] <- mean(lambda.star)-lambda.hat
    se.estimate[j,i] <- sd(lambda.star)
  }
}

## -----------------------------------------------------------------------------
bias.theta.hat <- rbind(apply(bias.estimate, 2, mean), lambda/(n-1))
se.theta.hat <- rbind(apply(se.estimate, 2, mean), lambda*n/(n-1)/sqrt(n-2))
colnames(bias.theta.hat) <- colnames(se.theta.hat) <- c("n=5", "n=10", "n=20")
rownames(bias.theta.hat) <- rownames(se.theta.hat) <- c("estimate", "theoretical")
pander::pander(bias.theta.hat)
pander::pander(se.theta.hat)

## -----------------------------------------------------------------------------
d <- bootstrap::law
d <- as.matrix(d)
n <- nrow(d)

B <-  1000 #bootstrap replicates for t
R <-  25 #bootstrap replicates for se
level <-  0.95
stat.star <- numeric(B)
se.star <- numeric(B)

## -----------------------------------------------------------------------------
b.cor <- function(x,i = 1:nrow(x)) cor(x[i,1],x[i,2])
boot.se <- function(x, R, statistic) {
  m <- nrow(x)
  th <- replicate(R, expr = {
    i <- sample(1:m, size = m, replace = TRUE)
    statistic(x[i, ])
  })
  return(sd(th))
}

## -----------------------------------------------------------------------------
set.seed(1)
for (b in 1:B) {
  j <- sample(1:n, size = n, replace = TRUE)
  xstar <- d[j, ]
  stat.star[b] <- b.cor(xstar)
  se.star[b] <- boot.se(xstar, R = R, statistic = b.cor)
#  se[b] <- sd(boot(xstar, R = R, statistic = b.cor)$t)
}
stat.hat <- b.cor(d)
t.stats <- (stat.star - stat.hat)/se.star
se.hat <- sd(stat.star)
alpha <- 1 - level
Qt <- quantile(t.stats, c(alpha/2, 1 - alpha/2), type = 1)
names(Qt) <- rev(names(Qt))
CI <- rev(stat.hat - Qt * se.hat)
CI

## -----------------------------------------------------------------------------
library(boot)
x <- aircondit[1]
mean.x <- function(x, i) return(mean(as.matrix(x[i, ])))

set.seed(1)
boot.x <- boot(x, statistic = mean.x, R = 2000)
boot.x
boot.ci(boot.x, type = c("norm", "perc", "basic", "bca"))

## -----------------------------------------------------------------------------
hist(boot.x$t, prob = TRUE, main = "", xlab = "")
points(boot.x$t0, 0, cex = 1, pch = 16)

## -----------------------------------------------------------------------------
detach(package:boot)
rm(list = ls())

## -----------------------------------------------------------------------------
x <- as.matrix(bootstrap::scor)
n <- nrow(x)
theta.jack <- numeric(n)

# original theta.hat
lambda <- eigen(cov(x))$values
theta.hat <- max(lambda/sum(lambda))

# n jackknife estimates theta.jack
for (i in 1:n) {
  x.jack <- x[-i, ]
  cov.jack <- cov(x.jack)
  lambda <- eigen(cov.jack)$values
  theta.jack[i] <- max(lambda/sum(lambda))
}
bias.jack <- (n - 1) * (mean(theta.jack) - theta.hat)
se.jack <- sqrt((n - 1)/n * sum((theta.jack - mean(theta.jack))^2))

# result reporting
result <- c(theta.hat, bias.jack, se.jack)
names(result) <- c("est", "bias", "se")
pander::pander(result)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
d <- DAAG::ironslag
d <- as.matrix(d)
n <- nrow(d)
N <- choose(n, 2)

e_lin <- e_quad <- e_exp <- e_loglog <- numeric(N)
ij <- 1
for (i in 1:(n - 1)) for (j in (i + 1):n) {
  k <- c(i, j)
  y <- d[-k, "magnetic"]
  x <- d[-k, "chemical"]
  
  # linear model
  J_lin <- lm(y ~ x)
  yhat_lin <- J_lin$coef[1] + J_lin$coef[2] * d[k, "chemical"]
  e_lin[ij] <- sum((d[k, "magnetic"] - yhat_lin)^2)
  
  # quadratic model
  J_quad <- lm(y ~ x + I(x^2))
  yhat_quad <- J_quad$coef[1] + J_quad$coef[2] * d[k, "chemical"] + J_quad$coef[3] * d[k, "chemical"]^2
  e_quad[ij] <- sum((d[k, "magnetic"] - yhat_quad)^2)
  
  # exponential model
  J_exp <- lm(log(y) ~ x)
  logyhat_exp <- J_exp$coef[1] + J_exp$coef[2] * d[k, "chemical"]
  yhat_exp <- exp(logyhat_exp)
  e_exp[ij] <- sum((d[k, "magnetic"] - yhat_exp)^2)
  
  # log-log model
  J_loglog <- lm(log(y) ~ log(x))
  logyhat_loglog <- J_loglog$coef[1] + J_loglog$coef[2] * log(d[k, "chemical"])
  yhat_loglog <- exp(logyhat_loglog)
  e_loglog[ij] <- sum((d[k, "magnetic"] - yhat_loglog)^2)

  ij <- ij + 1
}

# result reporting
spe <- cbind(e_lin, e_quad, e_exp, e_loglog)
aspe <- apply(spe, 2, mean)
names(aspe) <- c("linear model", "quadratic model", "exponential model", "log-log model")
pander::pander(aspe)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
cvm.test <- function(x, y, R = 199) {
  n <- length(x)
  m <- length(y)
  z <- c(x, y)
  N <- n + m
  
  # Fn: ecdf of the sample x1,...,xn; Gm: ecdf of the sample y1,...,ym.
  Fn <- numeric(N)
  Gm <- numeric(N)
  for(i in 1:N){
    Fn[i] <- mean(as.integer(z[i] <= x))
    Gm[i] <- mean(as.integer(z[i] <= y))
  }
  
  # statistic W2
  cvm0 <- ((n * m)/N) * sum((Fn - Gm)^2)
  
  # p-value of the test
  cvm <- numeric(R)
  for(j in 1:R){
    index <- sample(1:N)
    Z <- z[index]
    X <- Z[1:n]
    Y <- Z[(n + 1):N]
    for (i in 1:N) {
      Fn[i] <- mean(as.integer(Z[i] <= X))
      Gm[i] <- mean(as.integer(Z[i] <= Y))
    }
    cvm[j] <- ((n * m)/N) * sum((Fn - Gm)^2)
  }
  
  # result reporting
  cvm1 <- c(cvm, cvm0)
  return(c(cvm0, mean(cvm1 >= cvm0)))
}

## -----------------------------------------------------------------------------
attach(chickwts)
d <- chickwts
detach(chickwts)
d <- as.matrix(d)
d[d[,2]=="horsebean", 2] <- 1
d[d[,2]=="linseed", 2] <- 2
d[d[,2]=="soybean", 2] <- 3
d[d[,2]=="sunflower", 2] <- 4
d[d[,2]=="meatmeal", 2] <- 5
d[d[,2]=="casein", 2] <- 6

## -----------------------------------------------------------------------------
set.seed(100)
M2 <- matrix(NA, nrow = 6, ncol = 6)
p_value <- matrix(NA, nrow = 6, ncol = 6)
for(i in 1:5){
  for(j in (i+1):6){
    result <- cvm.test(d[d[,2]==i, 1], d[d[,2]==j, 1])
    M2[i,j] <- result[1]
    p_value[i,j] <- result[2]
  }
}

## -----------------------------------------------------------------------------
rownames(M2) <- colnames(M2) <- rownames(p_value) <- colnames(p_value) <- c("horsebean", "linseed", "soybean", "sunflower", "meatmeal", "casein")
M2
p_value

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
maxoutliers <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}

## -----------------------------------------------------------------------------
maxout <- function(x, y, R = 199) {
  z <- c(x, y)
  n <- length(x)
  N <- length(z)
  
  # statistic
  stat0 <- maxoutliers(x, y)
  
  # p-value of the test
  stat <- replicate(R, expr = {
    index <- sample(1:N)
    maxoutliers(z[index[1:n]], z[index[(n+1):N]])
  })
  
  stats <- c(stat0, stat)
  tab <- table(stats)/(R + 1)
  return(list(estimate = stat0, p = mean(stats >= stat0), freq = tab, cdf = cumsum(tab)))
}

## -----------------------------------------------------------------------------
set.seed(100)
n1 <- n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
maxout(x, y)

## -----------------------------------------------------------------------------
set.seed(100)
n1 <- 20
n2 <- 40
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
maxout(x, y)

## -----------------------------------------------------------------------------
set.seed(100)
n1 <- n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 2
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
maxout(x, y)

## -----------------------------------------------------------------------------
set.seed(100)
n1 <- 20
n2 <- 40
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 2
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
maxout(x, y)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
output.a <- function(N, b, f0){
  #生成N组样本(x1i, x2i, x3i), i=1,...,N.
  x <- rbind(rpois(N, lambda = 1), rexp(N, rate = 1), sample(0:1, N, replace=TRUE))
  
  #构造函数P(Y=1)-f0用来求解a
  fa_f0 <- function(a){
    prob <- 1/(1 + exp(-a-b[1]*x[1,]-b[2]*x[2,]-b[3]*x[3,]))
    log(mean(prob)) - log(f0)
  }
  
  #利用uniroot求解f(a)-f_0=0
  solution <- uniroot(fa_f0, c(-20,0))
  return(solution$root)
}

## -----------------------------------------------------------------------------
#按照题目要求指定 f0, b1, b2, b3, N (large enough)
N <- 1e6
b <- c(0, 1, -1)
f0 <- c(0.1, 0.01, 0.001, 0.0001)

#用前面定义的函数output.a得到a的估计
set.seed(1)
a <- sapply(f0, output.a, N = N, b = b)
names(a) <- c("f0=0.1", "f0=0.01", "f0=0.001", "f0=0.0001")
a

## -----------------------------------------------------------------------------
estimate.logf0 <- function(a, b){
  #生成N组样本(x1i, x2i, x3i), i=1,...,N.
  set.seed(1)
  x <- rbind(rpois(N, lambda = 1), rexp(N, rate = 1), sample(0:1, N, replace=TRUE))
  
  #P(Y=1)的概率
  prob <- 1/(1 + exp(-a-b[1]*x[1,]-b[2]*x[2,]-b[3]*x[3,]))
  
  #return(log(mean(p)))
  #用二项分布抽样估计概率
  return(log(mean(rbinom(N, 1, prob))))
}

## -----------------------------------------------------------------------------
#-logf0关于a的函数图像
a_axis <- seq(-12, 0, 0.1)
logf0 <- sapply(a_axis, estimate.logf0, b = b)

plot(a_axis, -logf0, type = 'l', xlab = "a", ylab = "-log f0")
for(i in 1:4){
  points(a[i], -log(f0[i]), cex = 1, pch = 16)
}

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
#Laplace分布密度函数
d.Laplace <- function(x){
  exp(-abs(x))/2
}

#Metropolis抽样方法，其中提议分布为正态分布
mcmc.Metropolis.normal <- function(sigma, x.initial, N, f) {
  x <- numeric(N) #M-H方法生成的样本链
  x[1] <- x.initial #初始值
  k <- 0  #记录接受Y的次数，计算接受率
  
  u <- runif(N)
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y)/f(x[i-1]))){x[i] <- y
      k <- k + 1
    } #接受
    else{x[i] <- x[i-1]} #拒绝
  }
  return(c(x, k/N))
}

## -----------------------------------------------------------------------------
N <- 5000 #样本链的长度
b <- 1000 #预烧期的长度
sigma <- c(0.5, 1, 4, 16) #给定的σ值
k <- length(sigma) #样本链的数量

## -----------------------------------------------------------------------------
#M-H方法生成4条样本链和对应的接受率
set.seed(1)
x.initial <- rnorm(1) #初始值
mcmc.result <- sapply(sigma, mcmc.Metropolis.normal, x.initial = x.initial, N = N, f = d.Laplace)
X.chain <- t(mcmc.result[1:N,]) #样本链
accept.prob <- as.vector(mcmc.result[(N+1),]) #接受率

## -----------------------------------------------------------------------------
x_axis <- seq(-6, 6, 0.01)
x_axis.Laplace <- d.Laplace(x_axis)
qq <- seq(0.005, 1, 0.005)
qq.population <- c(-rev(qexp(qq, 1)), qexp(qq, 1))
qq.chain <- apply(X.chain[,-(1:b)], 1, function(x) quantile(x, qq))

for(i in c(0,2)){
  par(mfrow = c(2,3))
  for(j in 1:2){
    #样本路径图
    plot(X.chain[i+j,], type = "l", main = bquote(sigma == .(sigma[i+j])), xlab = "t", ylab = "X", ylim = c(-8, 8))
    
    #样本频数分布直方图
    hist(X.chain[i+j, -(1:b)], breaks = "Scott", freq = FALSE, main = bquote(sigma == .(sigma[i+j])), xlim = c(-6, 6), ylim = c(0, 0.5), xlab = "x")
    lines(x_axis, x_axis.Laplace, col = 2)
    
    #样本Q-Q图
    qqplot(qq.population, qq.chain[,i+j], cex = 0.5, main = bquote(sigma == .(sigma[i+j])), xlab = "Sample quantile", ylab = "Population quantile")
    abline(a = 0, b = 1, col = 2)
  }
}

## -----------------------------------------------------------------------------
#生成服从拉普拉斯分布的样本，得到总体分位数的估计
set.seed(1)
x.Laplace <- rexp(200, 1)
x.Laplace <- x.Laplace * sample(c(-1,1), 200, replace = TRUE)
q <- c(0.05, seq(0.1, 0.9, 0.1), 0.95)
q.true <- quantile(x.Laplace, q)

#每条样本链的样本分位数
q.chain <- apply(X.chain[,-(1:b)], 1, function(x) quantile(x, q))

#比较总体分位数和样本链的样本分位数
q.table <- cbind(q.true, q.chain)
colnames(q.table) <- c('True', 'sigma=0.5', 'sigma=1', 'sigma=4', 'sigma=16')
knitr::kable(q.table)

## -----------------------------------------------------------------------------
#接受率
names(accept.prob) <- c("sigma=0.5", "sigma=1", "sigma=4", "sigma=16")
accept.prob

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
#Gibbs方法估计二元正态分布
mcmc.Gibbs.binormal <- function(mu, sigma, rho, initial, N) {
  X <- matrix(0, nrow = 2, ncol = N) #样本链
  X[,1] <- initial #初始值
  s <- sqrt(1 - rho^2) * sigma #条件正态分布的标准差
  
  for (t in 2:N) {
    mu1 <- mu[1] + rho*(X[2, t-1] - mu[2])*sigma[1]/sigma[2]
    X[1, t] <- rnorm(1, mu1, s[1])#从Y_{t-1}得到X_{t}
    mu2 <- mu[2] + rho*(X[1, t] - mu[1])*sigma[2]/sigma[1]
    X[2, t] <- rnorm(1, mu2, s[2])#从X_{t}得到Y_{t}
  }
  return(X)
}

## -----------------------------------------------------------------------------
N <- 5000 #length of chain
b <- 1000 #burn-in length
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)

#用Gibbs方法得到样本链
set.seed(1)
mcmc.result <- mcmc.Gibbs.binormal(mu, sigma, rho, mu, N)
#样本链X和Y的协方差
cov(t(mcmc.result))

## -----------------------------------------------------------------------------
#样本链X和Y的散点图
X.chain <- mcmc.result[1,]
Y.chain <- mcmc.result[2,]
plot(X.chain[-(1:b)], Y.chain[-(1:b)], xlab = "X", ylab = "Y", main = "Scatter diagram of (X,Y)", cex = 0.5, ylim = c(-5,5))
abline(h = 0, v = 0, lty = 2)

## -----------------------------------------------------------------------------
#Y关于X的线性回归模型
lin.reg <- lm(Y.chain[-(1:b)] ~ X.chain[-(1:b)])
lin.reg
summary(lin.reg)

## -----------------------------------------------------------------------------
#残差的期望和方差的理论值和估计值对比
mean.var <- c(0, mean(lin.reg$residuals), 0.19, var(lin.reg$residuals))
names(mean.var) <- c("theoretical mean", "estimated mean", "theoretical var", "estimated var")
pander::pander(mean.var)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
#残差的频数分布直方图
x_axis <- seq(-2, 2, 0.01)
hist(lin.reg$residuals, breaks = "Scott", freq = FALSE, main = "Density of residuals", xlim = c(-2, 2), ylim = c(0, 1), xlab = "Residual", ylab = "Density")
lines(x_axis, dnorm(x_axis, 0, sqrt(0.19)), col = 2)

#残差的Q-Q图
qqnorm(lin.reg$residuals, cex = 0.5)
qqline(lin.reg$residuals, col = 2)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
#Raleigh分布的密度
d.Raleigh <- function(x, sigma){
  if (any(x >= 0)){
    return((x/sigma^2) * exp(-x^2/(2*sigma^2)))
  }else{return(0)}
}

#M-H方法，其中提议分布为卡方分布
mcmc.MH.chisq <- function(sigma, N, x.initial, f){
  x <- numeric(N) #M-H方法生成的样本链
  x[1] <- x.initial #初始值
  u <- runif(N)
  
  for (i in 2:N) {
    y <- rchisq(1, df = x[i-1])
    numerator <- f(y, sigma) * dchisq(x[i-1], df = y)
    denominator <- f(x[i-1], sigma) * dchisq(y, df = x[i-1])
    if (u[i] <= numerator/denominator){
      x[i] <- y #接受
    }else{x[i] <- x[i-1]} #拒绝
  }
  return(x)
}

## -----------------------------------------------------------------------------
#Gelman-Rubin方法
Gelman.Rubin <- function(phi) {
  n <- ncol(phi)
  
  #计算Bn, Wn
  phi.i.mean <- apply(phi, 1, mean)
  B <- n * var(phi.i.mean)
  phi.i.var <- apply(phi, 1, var)
  W <- mean(phi.i.var)
  
  #计算Var.hat, R.hat
  var.hat <- W*(n - 1)/n + B/n
  r.hat <- var.hat/W
  return(r.hat)
}

## -----------------------------------------------------------------------------
sigma <- 4
x.initial <- c(0.5, 1, 4, 16) #给定初始值
k <- length(x.initial) #样本链的数量
N <- 2000 #样本链的长度
b <- 500 #预烧期的长度

#用M-H方法生成样本链
set.seed(1)
X.chain <- sapply(x.initial, mcmc.MH.chisq, sigma = sigma, N = N, f = d.Raleigh)
#处理样本链得到Φ_{ij}为第i个链前j项的均值
phi.sum <- t(apply(X.chain, 2, cumsum))
phi.mean <- t(apply(phi.sum, 1, function(x) x/(1:N)))

## -----------------------------------------------------------------------------
#用Gelman-Rubin方法监测收敛
r.hat.chain <- rep(0, N)
for (j in (b+1):N){
  r.hat.chain[j] <- Gelman.Rubin(phi.mean[, 1:j])
}

#R.hat的图像
plot(r.hat.chain[(b+1):N], type = "l", xlab = "t", ylab = "R.hat")
abline(h = 1.2, lty = 2)
#N=2000时的R.hat
r.hat.final <- Gelman.Rubin(phi.mean)
r.hat.final

## -----------------------------------------------------------------------------
library(coda)

## -----------------------------------------------------------------------------
X.chain <- t(X.chain)
#为了使用coda中的Gelman-Rubin诊断函数，我们将链的数据结构转换为mcmc型
mcmc.chain <- mcmc.list(as.mcmc(X.chain[1, ]), as.mcmc(X.chain[2, ]), as.mcmc(X.chain[3, ]), as.mcmc(X.chain[4, ]))
#Gelman-Rubin诊断
gelman.diag(mcmc.chain)
gelman.plot(mcmc.chain)

## -----------------------------------------------------------------------------
detach(package:coda)
rm(list = ls())

## -----------------------------------------------------------------------------
#观测数据对数似然的一阶导函数
d.log.likelihood_o <- function(lambda, u, v){
  numerator <- v*exp(-lambda*v) - u*exp(-lambda*u)
  denominator <- exp(-lambda*u) - exp(-lambda*v)
  return(sum(numerator/denominator))
}

#观测数据：x的下界
u <- c(11,8,27,13,16,0,23,10,24,2)
#观测数据：x的下界
v <- c(12,9,28,14,17,1,24,11,25,3)
#存储两种方法得到的估计值
lambda.estimate <- numeric(2)

## -----------------------------------------------------------------------------
#方法1：关于观测数据似然函数的MLE
#求解对数似然方程
result <- uniroot(d.log.likelihood_o, interval = c(0.01,0.1), u=u, v=v)
result
lambda.estimate[1] <- as.numeric(result[1])

## -----------------------------------------------------------------------------
#方法2：EM算法
n <- length(u) #样本容量
N <- 1e3 #迭代次数上界
lambda.seq <- numeric(N+1) #λ数列
lambda.seq[1] <- 0.01 #初始化
rate <- .Machine$double.eps #当数列中元素变化足够小判断收敛
#迭代
for(i in 1:N){
  lambda.seq[i+1] <- n/(n/lambda.seq[i] - d.log.likelihood_o(lambda.seq[i],u,v))
  if ((abs(lambda.seq[i+1] - lambda.seq[i])/lambda.seq[i]) < rate){
    cat("The sequence converges after", i, "iterations.", '\n')
    break
  }
}
lambda.estimate[2] <- lambda.seq[i+1]

## -----------------------------------------------------------------------------
#数据整理
names(lambda.estimate) <- c("MLE","EM")
lambda.estimate

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
Morra.game <- function(A) {
  #用单纯形法求解二人零和博弈
  
  #对A进行线性变换使得A中元素属于[0,1]，且目标函数v >= 0
  original.A <- A
  A_min <- min(A)
  A <- A - A_min #由此可得A中元素 >= 0，于是v >= 0
  A_max <- max(A)
  A <- A / A_max #由此可得A中元素属于[0,1]
  m.x <- nrow(A)
  n.y <- ncol(A)
  iteration <- n.y^3 #在单纯形法的每个阶段要进行的最大迭代次数
  
  #先优化玩家1，再优化玩家2
  
  #设玩家1的策略为x^*=(x_1,...,x_m)^T，把v作为额外变量，即要求最小值的目标函数，令x=(x_1,..., x_m, x_{m+1}=v)^T。最大化v=x_{m+1}=e_{m+1}^T x，受制于A1x <= 0, A3x = 1。
  object.x <- c(rep(0, m.x), 1) #目标函数：收益矩阵的值v=x_{m+1}=object.x^T x
  A1.x.leq <- -cbind(t(A), rep(-1, n.y)) #A1是由-A的转置与元素全为1的列合并得到的n*(m+1)维矩阵，即A1=-(A^T -1_n)
  b1.x.leq <- rep(0, n.y) #约束A_1x <= 0
  A3.x.eq <- t(as.matrix(c(rep(1, m.x), 0))) #A3是前m项为1，第m+1项为0的矩阵，即A3=(1_m^T 0)
  b3.x.eq <- 1 #约束A3x = 1，即sum(x^*)=1
  #利用函数simplex求解，解得x=(x_1,..., x_m, x_{m+1}=v)
  simplex.x <- simplex(a=object.x, A1=A1.x.leq, b1=b1.x.leq, A3=A3.x.eq, b3=b3.x.eq, maxi=TRUE, n.iter=iteration)
  
  #设玩家2的策略为y^*=(y_1,...,y_n)，把v作为额外变量，即要求最小值的目标函数，令y=(y_1,..., y_n, y_{n+1}=v)^T。最大化v=y_{n+1}=e_{n+1}^T y，受制于A1y <= 0, A3y = 1。
  object.y <- c(rep(0, n.y), 1) #目标函数：收益矩阵的值v=y_{n+1}=object.y^T y
  A1.y.leq <- cbind(A, rep(-1, m.x)) #A1是由A与元素全为-1的列合并得到的m*(n+1)维矩阵，即A1=(A 1_m)
  b1.y.leq <- rep(0, m.x) #约束A_1y <= 0
  A3.y.eq <- t(as.matrix(c(rep(1, n.y), 0))) #A3是前n项为1，第n+1项为0的矩阵，即A3=(1_n^T 0)
  b3.y.eq <- 1 #约束A3y = 1，即sum(y^*)=1
  #利用函数simplex求解，解得y=(x_1,..., y_n, y_{n+1}=v)
  simplex.y <- simplex(a=object.y, A1=A1.y.leq, b1=b1.y.leq, A3=A3.y.eq, b3=b3.y.eq, maxi=FALSE, n.iter=iteration)
  
  #结果整理，以列表形式输出收益矩阵，最优策略，以及Morra游戏的值
  result <- list("payoff matrix" = original.A, "strategy.x" = simplex.x$soln[1:m.x], "strategy.y" = simplex.y$soln[1:n.y], "value" = simplex.x$soln[m.x+1] * A_max + A_min)
  return(result)
}

## -----------------------------------------------------------------------------
A <- matrix(c(0,-2,-2,3,0,0,4,0,0,
              2,0,0,0,-3,-3,4,0,0,
              2,0,0,3,0,0,0,-4,-4,
              -3,0,-3,0,4,0,0,5,0,
              0,3,0,-4,0,-4,0,5,0,
              0,3,0,0,4,0,-5,0,-5,
              -4,-4,0,0,0,5,0,0,6,
              0,0,4,-5,-5,0,0,0,6,
              0,0,4,0,0,5,-6,-6,0), nrow = 9, ncol = 9) #收益矩阵
B <- A + 2 #收益矩阵每个元素都增加一个常数
library(boot) #包含函数simplex的包
#对A和B分别使用函数Morra.game，得到最优策略及Morra游戏的值
result.A <- Morra.game(A)
result.B <- Morra.game(B)

## -----------------------------------------------------------------------------
#结果整理
#比较游戏A和B的最优策略
compare.strategy <- round(cbind(result.A$strategy.x, result.A$strategy.y, result.B$strategy.x, result.B$strategy.y), 8)
colnames(compare.strategy) <- c("strategy 1 for A", "strategy 2 for A", "strategy 1 for B", "strategy 2 for B")
knitr::kable(compare.strategy)

#观察最优策略是式(11.12)-(11.15)中的哪一个
extreme.pionts <- cbind(rep(0,4), rep(0,4), c(5/12,16/37,20/47,25/61), rep(0,4), c(4/12,12/37,15/47,20/61), rep(0,4), c(3/12,9/37,12/47,16/61), rep(0,4), rep(0,4))
difference <- extreme.pionts - matrix(rep(compare.strategy[,1],4), nrow = 4, ncol = 9, byrow = TRUE)
#apply(extreme.pionts, 1, function(x) prod(all.equal(x,compare.strategy[,1])))
#apply(extreme.pionts, 1, function(x) identical(x,compare.strategy[,1]))
as.logical(apply(abs(difference)<1e-8, 1, prod))

#比较游戏A和B的值
compare.value <- c(result.A$value, result.B$value)
names(compare.value) <- c("value of game A", "value of game B")
pander::pander(compare.value)

## -----------------------------------------------------------------------------
detach(package:boot)
rm(list = ls())

## -----------------------------------------------------------------------------
test.list <- list(1:3)
test.list; str(test.list)
as.vector(test.list); str(as.vector(test.list))
1:3; str(1:3)
unlist(test.list); str(unlist(test.list))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
test.vector <- 1:3; str(test.vector)
dim(test.vector)
dim(matrix(test.vector, nrow = 1, ncol = 3))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
test.matrix <- matrix(1:3, nrow = 1, ncol = 3)
is.matrix(test.matrix)
is.array(test.matrix)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
test.matrix <- matrix(1:8, nrow = 4, ncol = 2)
colnames(test.matrix) <- c("col1", "col2") 
#只要行数与数据框匹配，数据框的列也可以是矩阵或数组。
#但是这里需要命名矩阵的列，否则运行as.matrix()会由于列数和列名长度不同而报错。

test.dataframe <- data.frame(a = 1:4, b = c(1.0,2.0,3.0,4.0), c = c(TRUE, FALSE, TRUE, FALSE), d = c("a", "b", "c", "d"), e = test.matrix, f = I(list(1:2, 1:3, c(1.0,2.0), c(1:2,3.0))), g = I(list(1:2, c("a", "b", "c"), c(TRUE, FALSE), test.matrix)))
test.dataframe; str(test.dataframe)
as.matrix(test.dataframe); str(as.matrix(test.dataframe))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
#使用函数data.frame()时只设置col.names无法生成0行的数据框。
row0.dataframe <- data.frame(col.names = c('a','b'))
row0.dataframe; dim(row0.dataframe)

#使用函数data.frame()时只设置row.names可以生成0列的数据框。
col0.dataframe <- data.frame(row.names = c('a','b'))
col0.dataframe; dim(col0.dataframe)

## -----------------------------------------------------------------------------
#通过将一个0行的矩阵转换为数据框，生成一个0行的数据框。
row0.matrix <- matrix(nrow = 0, ncol = 2)
row0.dataframe <- data.frame(row0.matrix)
row0.dataframe; dim(row0.dataframe)

#通过将一个0列的矩阵转换为数据框，生成一个0列的数据框。
col0.dataframe <- data.frame(t(row0.matrix))
col0.dataframe; dim(col0.dataframe)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
#生成一个每列都是int或double型的数据框
set.seed(1)
numeric.dataframe <- data.frame(a = 1:5, b = c(5.0, 4.0, 3.0, 2.0, 1.0), c = rnorm(5))
numeric.dataframe
#对数据框的每一列使用函数scale01()
data.frame(lapply(numeric.dataframe, scale01))

## -----------------------------------------------------------------------------
#生成列为int, double, logical, character或list的数据框
mixed.dataframe <- data.frame(a = 1:4, b = c(1.0,2.0,3.0,4.0), c = c(TRUE, FALSE, TRUE, FALSE), d = c("a", "b", "c", "d"), e = I(list(1:2, 1:3, c(1.0,2.0), c(1:2,3.0))), f = I(list(1:2, c("a", "b", "c"), c(TRUE, FALSE), matrix(1:4, nrow = 2, ncol = 2))))
mixed.dataframe
#对数据框的每一个numeric型的列使用函数scale01()
scale.dataframe <- data.frame(lapply(mixed.dataframe, function(x) if (is.numeric(x)) scale01(x) else x))
scale.dataframe; str(scale.dataframe)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
#生成一个每个元素都是numeric型的数据框
set.seed(1)
numeric.dataframe <- data.frame(a = 1:4, b = c(4.0, 3.0, 2.0, 1.0), c = rnorm(4))
numeric.dataframe
#对数据框的每一列使用函数sd()
vapply(numeric.dataframe, sd, FUN.VALUE = 1)

## -----------------------------------------------------------------------------
#生成列为int, double, logical, character或list的数据框
mixed.dataframe <- data.frame(a = 1:4, b = c(4.0,3.0,2.0,1.0), c = rnorm(4), d = c(TRUE, FALSE, TRUE, FALSE), e = c("a", "b", "c", "d"), f = I(list(1:2, 1:3, c(1.0,2.0), c(1:2,3.0))), g = I(list(1:2, c("a", "b", "c"), c(TRUE, FALSE), matrix(1:4, nrow = 2, ncol = 2))))
mixed.dataframe
#先判断数据框的numeric型的列有哪些，然后对每一个numeric型的列使用函数sd()
numeric.col <- vapply(mixed.dataframe, is.numeric, FUN.VALUE = TRUE)
vapply(mixed.dataframe[,numeric.col], sd, FUN.VALUE = 1)

## -----------------------------------------------------------------------------
#Gibbs方法估计二元分布
mcmc.Gibbs.R <- function(a, b, n, initial, N) {
  #a, b, n: 二元分布参数
  #initial: 初始值
  #N: 需要生成的样本链长度
  X <- matrix(0, nrow = 2, ncol = N) #样本链
  X[,1] <- initial #赋初值
  
  for (t in 2:N){
    #从Y_{t-1}得到X_{t}，(X_t|Y_{t-1})服从Binomial(n,Y_{t-1})分布
    X[1, t] <- rbinom(1, size = n, prob = X[2, t-1])
    #从X_{t}得到Y_{t}，(Y_t|X_t)服从Beta(X_t+a,n-X_t+b)分布
    X[2, t] <- rbeta(1, shape1 = X[1, t] + a, shape2 =  n - X[1, t] + b)
  }
  return(X)
}

## -----------------------------------------------------------------------------
library(Rcpp)

## -----------------------------------------------------------------------------
N <- 5000 #样本链长度
B <- 1000 #预烧期长度
a <- 2; b <- 5; n <- 9 #二元分布参数
initial <- c(1, 0.5)

#用R函数实现Gibbs方法得到样本链
set.seed(1)
mcmc.result.R <- mcmc.Gibbs.R(a, b, n, initial, N)
X.chain.R <- mcmc.result.R[1, -(1:B)]
Y.chain.R <- mcmc.result.R[2, -(1:B)]

#用Rcpp函数实现Gibbs方法得到样本链
set.seed(1)
sourceCpp('../src/mcmc.Gibbs.Rcpp.cpp')
mcmc.result.Rcpp <- mcmc_Gibbs_Rcpp(a, b, n, initial, N)
X.chain.Rcpp <- mcmc.result.Rcpp[1, -(1:B)]
Y.chain.Rcpp <- mcmc.result.Rcpp[2, -(1:B)]

## -----------------------------------------------------------------------------
#列表比较关于X的样本链的分布和X的理论边缘分布
#X样本链分布
X.chain.R.freq <- table(X.chain.R) / (N-B)
X.chain.Rcpp.freq <- table(X.chain.Rcpp) / (N-B)
#X理论分布
X.value <- 0:9
X.population <- choose(n, X.value) * beta(X.value+a, n-X.value+b) / beta(a, b)

X.compare <- rbind(X.chain.R.freq, X.chain.Rcpp.freq, X.population)
rownames(X.compare) <- c("Samples from R function", "Samples from Rcpp function", "Theoretical distribution")
round(X.compare,4)

## -----------------------------------------------------------------------------
par(mfrow = c(1,2))
#作图比较关于由R函数生成的X的样本链的分布和X的理论边缘分布
barplot(X.chain.R.freq, space = 0, main = "Density of X from R function", xlim = c(0, 9), ylim = c(0, 0.25), xlab = "X", ylab = "Density") #X样本链分布
points(0:n + 0.5, X.population, col = 2, pch = 20) #X理论分布

#作图比较关于由R函数生成的Y的样本链的分布和Y的理论边缘分布
Y_axis <- seq(0,1,0.01)
hist(Y.chain.R, breaks = "Scott", freq = FALSE, main = "Density of Y from R function", xlim = c(0, 1), ylim = c(0, 2.5), xlab = "Y", ylab = "Density") #Y样本链分布
lines(Y_axis, Y_axis^(a-1) * (1-Y_axis)^(b-1) / beta(a,b), col = 2) #Y理论分布

par(mfrow = c(1,2))
#作图比较关于由Rcpp函数生成的X的样本链的分布和X的理论边缘分布
barplot(X.chain.Rcpp.freq, space = 0, main = "Density of X from Rcpp function", xlim = c(0, 9), ylim = c(0, 0.25), xlab = "X", ylab = "Density") #X样本链分布
points(0:n + 0.5, X.population, col = 2, pch = 20) #X理论分布

#作图比较关于由Rcpp函数生成的Y的样本链的分布和Y的理论边缘分布
Y_axis <- seq(0,1,0.01)
hist(Y.chain.Rcpp, breaks = "Scott", freq = FALSE, main = "Density of Y from Rcpp function", xlim = c(0, 1), ylim = c(0, 2.5), xlab = "Y", ylab = "Density") #Y样本链分布
lines(Y_axis, Y_axis^(a-1) * (1-Y_axis)^(b-1) / beta(a,b), col = 2) #Y理论分布

## -----------------------------------------------------------------------------
library(microbenchmark)
compare.time <- microbenchmark(R.function = mcmc.Gibbs.R(a, b, n, initial, N), Rcpp.function = mcmc_Gibbs_Rcpp(a, b, n, initial, N))
summary(compare.time)[,c(1,3,5,6)]

## -----------------------------------------------------------------------------
detach(package:microbenchmark)
rm(list = ls())

