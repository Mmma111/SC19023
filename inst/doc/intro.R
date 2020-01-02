## ----eval=FALSE---------------------------------------------------------------
#  linear_reg <- function(formula, data){
#    x <- scale(data)
#    result <- lm(formula, x)
#    return(result)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  rwMetropolisR <- function(sigma, x0, N) {
#    x <- numeric(N)
#    x[1] <- x0
#    u <- runif(N)
#    k <- 0
#    for (i in 2:N) {
#      y <- rnorm(1, x[i-1], sigma)
#      if (u[i] <= exp(abs(x[i-1]) - abs(y)))
#        x[i] <- y else {
#          x[i] <- x[i-1]
#          k <- k + 1
#        }
#    }
#    return(list(x=x))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  library(Rcpp)
#  cppFunction('NumericVector rwMetropolisC(double sigma, double x0, double N) {
#    NumericVector x(N+1);
#    x[0] = x0;
#    NumericVector u = runif(N);
#    double k = 0;
#    double y = 0;
#    for (int i = 2; i < N+1; i++) {
#      y = rnorm(1, x[i-2], sigma)[0];
#      if (u[i-2] <= exp(-((abs(y))-(abs(x[i-2])))) ) {
#        x[i-1] = y;
#      }
#      else {
#        x[i-1] = x[i-2];
#        k++;
#      }
#    }
#    x[N] = k;
#    return(x);
#  }')

## ----eval=TRUE----------------------------------------------------------------
library(SC19023)
library(microbenchmark)

tm <- microbenchmark(
  vR = rwMetropolisR( .5, 25, 2000),
  vC = rwMetropolisC( .5, 25, 2000))
knitr::kable(summary(tm)[, c(1,3,5,6)])

## -----------------------------------------------------------------------------
a <- "Hello,USTC!"
print(a)

## -----------------------------------------------------------------------------
x <- rnorm(10)
y <- rnorm(10)
plot(x,y)

## -----------------------------------------------------------------------------
m <- 1:4;n <- 8
data.frame(m,n)

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
sigma <- c(1)
x <- {-2*sigma*log(1-u)}^(1/2) #F(x)=1-exp{-x^2/2*sigma^2}
hist(x,prob=TRUE)
y <- seq(0,4, .01)
lines(y,exp(-1*y^2/2*sigma^2)*y/sigma^2)

## -----------------------------------------------------------------------------
n <- 1000
x1 <- rnorm(n,0,1)
x2 <- rnorm(n,3,1)
u <- runif(n)
k <- as.integer(u>0.25)
x <- k * x1 + (1-k) * x2
hist(x, prob=TRUE)
y <- seq(-4,6, .01)
lines(y,0.75/(2*pi)^1/2*exp(-y^2/2)+0.25/(2*pi)^1/2*exp(-(y-3)^2/2))

## -----------------------------------------------------------------------------
n <- 1000


## -----------------------------------------------------------------------------
m <- 10000
x <- runif(m, min=0, max=(pi/3))
theta.hat <- mean(sin(x)) * (pi/3)
print(theta.hat)
print(1 - cos(pi/3))

## -----------------------------------------------------------------------------

R <- 10000
u <- runif(R/2)
v <- 1 - u
g1 <- exp(-u)/(1+u^2)
g2 <- exp(-v)/(1+v^2)
theta.hat1 <- mean((g1+g2)/2) 
theta.hat1


## -----------------------------------------------------------------------------

M <- 10000
x <- runif(M)
f <- exp(-x)/(1+x^2)
theta.hat2 <- mean(f)
theta.hat2


## -----------------------------------------------------------------------------
m <- 1000
theta.hat1 <- theta.hat2 <- numeric(m)
var_theta.hat1 <- (sd(g1)+sd(g2)+cov(g1,g2))/4
var_theta.hat2 <- sd(f)
print (var_theta.hat1)
print (var_theta.hat2)
print((var_theta.hat2 - var_theta.hat1)/var_theta.hat2)

## -----------------------------------------------------------------------------
M <- 5000 #number of replicates
k <- 5 #number of strata
r <- M / k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
estimates <- matrix(0, N, 2)

g <- function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

fg <- function(x) {
g(x)/(exp(-x) / (1 - exp(-1)) * (x > 0) * (x < 1))
}

for (i in 1:N) {
  estimates[i, 1] <- mean(fg(runif(M)))
  for (j in 1:k)
    T2[j] <- mean(fg(runif(M/k, (j-1)/k, j/k)))
  estimates[i, 2] <- mean(T2)
}

## -----------------------------------------------------------------------------
apply(estimates, 2, mean)
apply(estimates, 2, var)

## -----------------------------------------------------------------------------


n <- 20
alpha <-  .05
ULC <- replicate(1000,expr={
  x <- rt(n,df=2)
  (n-1)*var(x)/qchisq(alpha,df=n-1)
})
sum(ULC > 4)
mean(ULC > 4)



## -----------------------------------------------------------------------------


sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

m <- 10000
n <- 1000
sqrt.b1 <- numeric(m)
for (j in 1:m) {
x <- rnorm(n)
sqrt.b1[j] <- sk(x)
}

q <- c(0.025,0.05,0.95,0.975)
quantile(sqrt.b1, probs = q)




## -----------------------------------------------------------------------------


n <- 1000
q <- c(0.025,0.05,0.95,0.975)
cv <- qnorm(q, 0, sqrt(6/n))
cv



## -----------------------------------------------------------------------------
xq <- as.numeric(quantile(sqrt.b1, probs = q))
sd_sqrt.b1 <- sqrt(6*(n-2)/((n+1)*(n+3)))

sd.hat <- numeric(4)
for(i in 1:4){
f <- dnorm(xq[i], mean = 0, sd = sd_sqrt.b1)
sd.hat[i] <- q[i]*(1-q[i])/(n*f^2)
}

sd.hat

## -----------------------------------------------------------------------------

sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

b <- .1
n <- 30
m <- 2500
alpha <- c(1:31)
N <- length(alpha)
pwr <- numeric(N)
cv <- qnorm(1-b/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  a <- alpha[j]
  sktests <- numeric(m)
  for (i in 1:m) { 
    x <- rbeta(n, a, a)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
    }
pwr[j] <- mean(sktests)
}

pwr

plot(alpha, pwr, type = "b",
xlab = bquote(alpha), ylim = c(0, .1))
abline(h = .01, lty = 3)


## -----------------------------------------------------------------------------

sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

b <- .1
n <- 30
m <- 2500
v <- c(1:31)
N <- length(v)
pwr <- numeric(N)
cv <- qnorm(1-b/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  a <- as.integer(v[j])
  sktests <- numeric(m)
  for (i in 1:m) { 
    x <- rt(n, a)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
    }
pwr[j] <- mean(sktests)
}

pwr

plot(v, pwr, type = "b",
xlab = bquote(v), ylim = c(0,1))
abline(h = .01, lty = 3)


## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
I <- numeric(m) 
cv1 <- qchisq(1-alpha/2, 1)
cv2 <- qchisq(alpha/2, 1)

for (j in 1:m) {
  chi <- rchisq(n, 1)
  I[j] <- as.integer(chi >= cv1||chi <= cv2)
}
p.hat <- mean(I)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)

print(c(p.hat,se.hat))



## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
I <- numeric(m) 
cv1 <- qunif(1-alpha/2, 0, 2)
cv2 <- qunif(alpha/2, 0, 2)

for (j in 1:m) {
  chi <- runif(n, 0, 2)
  I[j] <- as.integer(chi >= cv1||chi <= cv2)
}
p.hat <- mean(I)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)

print(c(p.hat,se.hat))


## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
I <- numeric(m) 
cv1 <- qexp(1-alpha/2, 1)
cv2 <- qexp(alpha/2, 1)

for (j in 1:m) {
  chi <- rexp(n, 1)
  I[j] <- as.integer(chi >= cv1||chi <= cv2)
}
p.hat <- mean(I)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)

print(c(p.hat,se.hat))


## -----------------------------------------------------------------------------

library(bootstrap)
library(ggplot2)

ggplot(data = scor)+
  geom_point(mapping = aes(x = 1:88, y = scor$mec, color="scor$mec"))+
  geom_point(mapping = aes(x = 1:88, y = scor$vec, color="scor$vec"))+
  geom_point(mapping = aes(x = 1:88, y = scor$alg, color="scor$alg"))+
  geom_point(mapping = aes(x = 1:88, y = scor$ana, color="scor$ana"))+
  geom_point(mapping = aes(x = 1:88, y = scor$sta, color="scor$sta"))



## -----------------------------------------------------------------------------
r12 <- function(x, i) {cor(x[i,1], x[i,2])}
r34 <- function(x, i) {cor(x[i,3], x[i,4])}
r35 <- function(x, i) {cor(x[i,3], x[i,5])}
r45 <- function(x, i) {cor(x[i,4], x[i,5])}


library(boot) 
obj12 <- boot(data = scor, statistic = r12, R = 2000)
obj34 <- boot(data = scor, statistic = r34, R = 2000)
obj35 <- boot(data = scor, statistic = r35, R = 2000)
obj45 <- boot(data = scor, statistic = r45, R = 2000)


cor12 <- mean(obj12$t)
cor34 <- mean(obj34$t)
cor35 <- mean(obj35$t)
cor45 <- mean(obj45$t)

c(cor12,cor34,cor35,cor45)


## -----------------------------------------------------------------------------

library(boot)

n <- 100
x <- rnorm(n, 0, 1)

theta.boot <- function(x,i) {
xbar <- mean(x[i])
m3 <- mean((x[i] - xbar)^3)
m2 <- mean((x[i] - xbar)^2)
return( m3 / m2^1.5 )
}

boot.obj1 <- boot(x, R = 2000, statistic = theta.boot)

print(boot.obj1)
print(boot.ci(boot.obj1, type=c("basic","norm","perc")))


## -----------------------------------------------------------------------------

library(boot)

n <- 100
x <- rchisq(n, 5)

theta.boot <- function(x,i) {
xbar <- mean(x[i])
m3 <- mean((x[i] - xbar)^3)
m2 <- mean((x[i] - xbar)^2)
return( m3 / m2^1.5 )
}

boot.obj2 <- boot(x, R = 2000, statistic = theta.boot)

print(boot.obj2)
print(boot.ci(boot.obj2, type=c("basic","norm","perc")))


## -----------------------------------------------------------------------------
library(bootstrap)
data(scor)
scor_pca <- princomp(scor, cor = F)
scor_pca_summary <- summary(scor_pca)

theta.hat <- 0.619115


n <- nrow(scor)
theta.hat.jackknife <- numeric(n)
for (i in 1:n) {
scor_jackknife <- scor[-i, ]
scor_jackknife_cov <- cov(scor_jackknife)
theta.hat.jackknife[i] <-
eigen(scor_jackknife_cov)$value[1] / sum(eigen(scor_jackknife_cov)$value)
}


bias_theta <- (n - 1) * (mean(theta.hat.jackknife) - theta.hat)
 

se_theta <- sqrt((n - 1) * mean((theta.hat.jackknife - mean(theta.hat.jackknife)) ^ 2))

output = matrix(0,1,2)
colnames(output) <- c("bias_theta", "se_theta")
output[1,] <- c(bias_theta, se_theta)

output


## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
#data(ironslag, package = "DAAG")

a <- seq(10, 40, .1) 

par(mar = c(4,3,3,1)+0.1,mfrow=c(2,2))

L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
plot(chemical, magnetic, main="Cubic Polynomial", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)


n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- numeric(n)
ey1 <- ey2 <- ey3 <- ey4 <- numeric(n)
ybar <- mean(magnetic)


for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]

J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1
ey1[k] <- magnetic[k] - ybar

J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
ey2[k] <- magnetic[k] - ybar

J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
ey3[k] <- magnetic[k] - ybar

J4 <- lm(y ~ x + I(x^2) + I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * (chemical[k]^2) + J4$coef[4] * (chemical[k]^3)
e4[k] <- magnetic[k] - yhat4
ey4[k] <- magnetic[k] - ybar
}

sigmahat2 <- c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

r2 <- c((1-sum(e1^2)/sum(ey1^2)),(1-sum(e2^2)/sum(ey2^2)),(1-sum(e3^2)/sum(ey3^2)),(1-sum(e4^2)/sum(ey4^2)))

a <- matrix(0,2,4)
a[1,] <- sigmahat2
a[2,] <- r2
colnames(a) <- c("Linear","Quadratic","Exponential","Cubic Polynomial")
rownames(a) <- c("sigmahat^2","r^2")
a

L2


## -----------------------------------------------------------------------------

#Count Seven test statistics
maxout <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}


count7test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(as.integer(max(c(outx, outy)) > 7))
}


R <- 999 
set.seed(12345)
x <- rnorm(20, 0, sd = 1)
y <- rnorm(30, 0, sd = 1)
z <- c(x, y) 
K <- 1:50

test <- numeric(R)
stat <- numeric(R)

for (i in 1:R) {
  k <- sample(K, size = 20, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  stat[i] <- maxout(x1,y1)
  x1 <- z[k]-mean(x1)
  y1 <- z[-k]-mean(y1)
  test[i] <- count7test(x1,y1)
}


#The empirica quantiles
print(quantile(stat,c( .8, .9, .95)))

#empirical Type I error rate 
alphahat <- mean(test)
alphahat

## -----------------------------------------------------------------------------
library(boot)
library(Ball)
library(MASS)
library(ggplot2)

dCov <- function(x, y) { 
  x <- as.matrix(x) 
  y <- as.matrix(y) 
  n <- nrow(x) 
  m <- nrow(y) 
  if (n != m || n < 2) stop("Sample sizes must agree") 
  if (! (all(is.finite(c(x, y))))) stop("Data contains missing or infinite values")
  
  Akl <- function(x) {
    d <- as.matrix(dist(x)) 
    m <- rowMeans(d) 
    M <- mean(d) 
    a <- sweep(d, 1, m) 
    b <- sweep(a, 2, m) 
    return(b + M) 
    } 
  A <- Akl(x) 
  B <- Akl(y) 
  dCov <- sqrt(mean(A * B)) 
  dCov
}


ndCov2 <- function(z, ix, dims) { #dims contains dimensions of x and y 
  p <- dims[1] 
  q <- dims[2] 
  d <- p + q 
  x <- z[ , 1:p] #leave x as is 
  y <- z[ix, -(1:p)] #permute rows of y 
  return(nrow(z) * dCov(x, y)^2) 
  } 

alpha <- 0.01
n <- c(seq(50,120,5))
m <- 50


set.seed(12345) 


pow_dCor_Model1<-pow_ball_Model1<-pow_dCor_Model2<-pow_ball_Model2<-numeric(length(n))
dcor1<-dcor2<-p_ball1<-p_ball2<-numeric(m)

for (j in 1:length(n)) {
  
  
  for (i in 1:m) {
    
    set.seed(i)
    
    X<-mvrnorm(n[j],rep(0,2),diag(2))
    e<-mvrnorm(n[j],rep(0,2),diag(2))
    Y1<-(X/4)+e
    Z1<-cbind(X,Y1)
    Y2<-(X/4)*e
    Z2<-cbind(X,Y2)
   
    #model 1
    t1<-bcov.test(X,Y1,R=100)
    p_ball1[i]<-t1$p.value
    boot.obj1<-boot(data=Z1,statistic=ndCov2,R=100,sim="permutation",dims=c(2, 2))
    temp1<-c(boot.obj1$t0, boot.obj1$t)
    dcor1[i]<-mean(temp1>=temp1[1])
    
    #model 2
    t2<-bcov.test(X,Y2,R=100)
    p_ball2[i]<-t2$p.value
    boot.obj2<-boot(data=Z2,statistic=ndCov2,R=100,sim="permutation",dims=c(2, 2))
    temp2<-c(boot.obj2$t0, boot.obj2$t)
    dcor2[i]<-mean(temp2>=temp2[1])
    
    }
  pow_dCor_Model1[j]<-mean(dcor1<alpha)
  pow_ball_Model1[j]<-mean(p_ball1<alpha)
  pow_dCor_Model2[j]<-mean(dcor2<alpha)
  pow_ball_Model2[j]<-mean(p_ball2<alpha)
}

dat1<-data.frame(pow_dCor_Model1,pow_ball_Model1)

ggplot(dat1,aes(n))+geom_point(y=pow_dCor_Model1,fill="white")+geom_line(y=pow_dCor_Model1,colour="red")+geom_point(y=pow_ball_Model1,fill="white")+geom_line(y=pow_ball_Model1,colour="green")



dat2<-data.frame(pow_dCor_Model2,pow_ball_Model2)

ggplot(dat2,aes(n))+geom_point(y=pow_dCor_Model2,fill="white")+geom_line(y=pow_dCor_Model2,colour="red")+geom_point(y=pow_ball_Model2,fill="white")+geom_line(y=pow_ball_Model2,colour="green")

## -----------------------------------------------------------------------------

dt <- function(t){
  exp(-abs(t))
}
 
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (dt(y) / dt(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
        }
    }
  return(list(x=x, k=k))
  }

N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

rw <- c(rw1$k, rw2$k, rw3$k, rw4$k)
rej <- rw/N
rec <- 1-rej
m <- matrix(0, 2, 4)
rownames(m) <- c("rejection","reception") 
colnames(m) <- NULL
m[1,] <- rej
m[2,] <- rec
m

a <- 1:2000
par(mar = c(4,4,4,2)+0.1,mfcol=c(2,2))
plot(a, rw1$x, type = "s", main = "", xlab = "sigma = 0.05", ylab = "x")
plot(a, rw2$x, type = "s", main = "", xlab = "sigma = 0.5", ylab = "x")
plot(a, rw3$x, type = "s", main = "", xlab = "sigma = 2", ylab = "x")
plot(a, rw4$x, type = "s", main = "", xlab = "sigma = 16", ylab = "x")




## -----------------------------------------------------------------------------
x <- numeric(100)

isTRUE(identical(log(exp(x)),exp(log(x))))

isTRUE(all.equal(log(exp(x)),exp(log(x))))


## -----------------------------------------------------------------------------

k <- 100

g <- function(v){
  return(exp(log(2)+lgamma(v/2)-(1/2)*(log(pi * (v-1)))-lgamma((v-1)/2)))
  }

f <- function(a)
  (g(k)*integrate(function(u)(1 + (u^2)/(k-1))^(-k/2),lower = 0,upper =sqrt((a^2)*(k-1)/(k-a^2)))$value)-(g(k+1)*integrate(function(u)(1+(u^2)/k)^(-(k+1)/2),lower = 0,upper = sqrt((a^2)*k/(k+1-a^2)))$value)

out <- uniroot(f,lower = 0.1,upper = 1+sqrt(k)/2)

out


## -----------------------------------------------------------------------------
n_A <- 28; n_B <- 24; n_OO <- 41; n_AB <- 70

p <- q <- r <- numeric(100)

p[1] <- 0.2
q[1] <- 0.2
r[1] <- (1- p[1]- q[1])
   
f1 <- function(a,b) {
   return((n_B*b/(2-b-2*a)+n_B+n_AB)/(n_A*a/(2-a-2*b)+n_A+n_AB))
}

f2 <- function(a,b) {
   return(((1-a/(2-a-2*b))*n_A+(1-b/(2-b-2*a))*n_B+2*n_OO)/((n_B*b/(2-b-2*a)+
                          n_B+n_AB)))
}
threshold <- 1e-5


for (k in 2:100) {
   p[k] <- 1/(1+f1(p[k-1],q[k-1])*(1+f2(p[k-1],q[k-1])))
   q[k] <- f1(p[k-1],q[k-1])/(1+f1(p[k-1],q[k-1])*(1+f2(p[k-1],q[k-1])))
   r[k] <- 1- p[k] - q[k]
   #The iterative relation is obtained by EM algorithm
   if((p[k]-p[k-1] <= threshold) & (q[k]-q[k-1] <= threshold) & (r[k]-r[k-1] <= threshold))
      {print(c(k, p[k], q[k],r[k]))
      break
      }
   }

x <- seq(1,k,1)

plot(x, p[1:k], "b", col = "red",ylim=c(0,0.6), xlab = "Iterations", ylab = "Iteration value", main = "The log-maximum likelihood values in M-steps")
lines(x, q[1:k], "b", col = "green")
lines(x, r[1:k], "b", col = "blue")
legend("topright", legend = c("p", "q", "r"),lty = 1, col = c("red", "green", "blue"))


## -----------------------------------------------------------------------------

formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

f1 <- lapply(formulas, lm, data = mtcars)
f2 <- lapply(formulas, function(x) lm(formula = x, data = mtcars))

#loop
lf1 <- vector("list", length(formulas))
for (i in seq_along(formulas)){
  lf1[[i]] <- lm(formulas[[i]], data = mtcars)
}



## -----------------------------------------------------------------------------

bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

#without an anonymous function
f <- lapply(bootstraps, lm, formula = mpg ~ disp)

#loop
lf <- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)){
  lf[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}



## -----------------------------------------------------------------------------

rsq <- function(mod) summary(mod)$r.squared

#in the exercise3
sapply(f1, rsq)
sapply(f2, rsq)
sapply(lf1, rsq)

#in the exercise4
sapply(f, rsq)
sapply(lf, rsq)



## -----------------------------------------------------------------------------

trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

# anonymous function
sapply(trials, function(x) x[["p.value"]])

# without anonymous function
sapply(trials, "[[", "p.value")



## -----------------------------------------------------------------------------

f <- function(x){
  return(x+1)
}

#parSapply
library(parallel)

cl <- makeCluster(getOption('cl.cores', 4));
parSapply(cl, 1:500, f)
stopCluster(cl)




## -----------------------------------------------------------------------------

library(Rcpp)

cppFunction('NumericVector rwMetropolisC(double sigma, double x0, double N) {
  NumericVector x(N+1);
  x[0] = x0;
  NumericVector u = runif(N);
  double k = 0;
  double y = 0;
  for (int i = 2; i < N+1; i++) {
    y = rnorm(1, x[i-2], sigma)[0];
    if (u[i-2] <= exp(-((abs(y))-(abs(x[i-2])))) ) {
      x[i-1] = y;
    } 
    else {
      x[i-1] = x[i-2];
      k++;
    }
  }
  x[N] = k;
  return(x);
}')

N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
set.seed(12345)
rwC1 <- rwMetropolisC(sigma[1], x0, N)
rwC2 <- rwMetropolisC(sigma[2], x0, N)
rwC3 <- rwMetropolisC(sigma[3], x0, N)
rwC4 <- rwMetropolisC(sigma[4], x0, N)
rwC <- cbind(rwC1[-(N+1)], rwC2[-(N+1)], rwC3[-(N+1)], rwC4[-(N+1)])

#the acceptance rates
print(c(1-rwC1[N+1]/N,1-rwC2[N+1]/N,1-rwC3[N+1]/N,1-rwC4[N+1]/N))

par(mar = c(4,4,4,2)+0.1,mfcol = c(2, 2)) 
refline <- c(log(2 * .025), -log(2 * (1 - .975)))
for (j in 1:4) {
  plot(rwC[, j], type = 'l', xlab = bquote(sigma == .(round(sigma[j], 3))), ylab = 'X', ylim = range(rwC[, j]))
  abline(h = refline)
}



## -----------------------------------------------------------------------------


dt <- function(t){
  exp(-abs(t))
}
 

rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (dt(y) / dt(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
        }
    }
  return(list(x=x, k=k))
  }

N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

rw <- c(rw1$k, rw2$k, rw3$k, rw4$k)
rej <- rw/N
rec <- 1-rej
m <- matrix(0, 2, 4)
rownames(m) <- c("rejection","reception") 
colnames(m) <- NULL
m[1,] <- rej
m[2,] <- rec
m

a <- 1:2000
par(mar = c(4,4,4,2)+0.1,mfcol=c(2,2))
plot(a, rw1$x, type = "s", main = "", xlab = "sigma = 0.05", ylab = "x")
plot(a, rw2$x, type = "s", main = "", xlab = "sigma = 0.5", ylab = "x")
plot(a, rw3$x, type = "s", main = "", xlab = "sigma = 2", ylab = "x")
plot(a, rw4$x, type = "s", main = "", xlab = "sigma = 16", ylab = "x")




## -----------------------------------------------------------------------------

par(mar = c(4,4,4,2)+0.1,mfcol = c(2, 2))
qqplot(rw1$x, rwC1[-(N+1)], xlab = 'rw1', ylab = 'rwC1')
qqplot(rw2$x, rwC2[-(N+1)], xlab = 'rw2', ylab = 'rwC2')
qqplot(rw3$x, rwC3[-(N+1)], xlab = 'rw3', ylab = 'rwC3')
qqplot(rw4$x, rwC4[-(N+1)], xlab = 'rw4', ylab = 'rwC4')



## -----------------------------------------------------------------------------

library(microbenchmark)

ts <- microbenchmark(rwC1 = rwMetropolisC(sigma[1], x0, N), rwC2 = rwMetropolisC(sigma[2], x0, N), rwC3 = rwMetropolisC(sigma[3], x0, N), rwC4 = rwMetropolisC(sigma[4], x0, N), rw1 = rw.Metropolis(sigma[1], x0, N), rw2 = rw.Metropolis(sigma[2], x0, N), rw3 = rw.Metropolis(sigma[3], x0, N), rw4 = rw.Metropolis(sigma[4], x0, N))
knitr::kable(summary(ts)[, c(1,3,5,6)])


