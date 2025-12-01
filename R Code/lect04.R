rm(list=ls())
setwd("~/Kurser/02418/2023Fall/lect04")
##################################################
## Example 4.16
y<-c(8,10,5,5,59,82,10,57,4,0.4,16,10,39,9,59,
     8,.7,10,1,10,4,2,0.4,38,22,5,2)

qqnorm(y)
qqline(y)

## box-cox transformation
bc.trans <- function(lambda,y){
    y.l <- (y^lambda-1)/lambda
    if(lambda==0){y.l <- log(y)}
    return(y.l)
}


## Try different lambda's

par(mfrow=c(2,2))
lambda <- 1
qqnorm(bc.trans(lambda ,y))
qqline(bc.trans(lambda ,y))

lambda <- 1/2
qqnorm(bc.trans(lambda ,y))
qqline(bc.trans(lambda ,y))

lambda <- 0
qqnorm(bc.trans(lambda ,y))
qqline(bc.trans(lambda ,y))

lambda <- 1/4
qqnorm(bc.trans(lambda ,y))
qqline(bc.trans(lambda ,y))

## 1/4 seems to be a good option


## profile likelihood for lambda
lp.lambda <- function(lambda,y){
    n <- length(y)
    y.l <- bc.trans(lambda ,y)
    sigmasq <- 1/n * sum((y.l-mean(y.l))^2)
    -n/2 * log(sigmasq) + (lambda-1)*sum(log(y))
}



par(mfrow=c(1,2))
## Plot the profile likelihood
lambda <- seq(-0.5,0.5,by=0.01)

lp <- sapply(lambda,lp.lambda, y=y)
plot(lambda,lp-max(lp),type="l")
lines(range(lambda),-qchisq(0.95,df=1)/2*c(1,1),lty=2,col=2)
## so lambda in the range -0.2 to 0.4 (with 95% confidence)

## Directly in R
library(MASS)
boxcox(lm(y~1),lambda=lambda)

## and finally we could optimize it
optimize(lp.lambda,c(-2,2),y=y,maximum=TRUE)
## Too find a lambda of 0.099 or 1/10

lambda <- 1/4
qqnorm(bc.trans(lambda ,y))
qqline(bc.trans(lambda ,y))

qqnorm(bc.trans(0.1 ,y))
qqline(bc.trans(0.1 ,y))

##################################################
##
##################################################


##################################################
## Example 4.17
##################################################

x<- c(-26.8, -3.5, -3.4, -1.2,  0.4, 1.3, 2.3, 2.7,
      3.0 , 3.2,  3.2,  3.5,  3.6, 3.9, 4.2, 4.4,
      5.0 , 6.5,  6.7,  7.1,  8.1, 10.5,10.7, 24.0, 32.8)

qqnorm(x)
qqline(x)
## Normal model poor fit with extremes

## likelihood for Gaussian model
n <- length(x)
s2 <- var(x) * (n - 1)/n ## MLE of sigma^2
normal.ll <- sum(dnorm(x, mean = mean(x),
                       sd = sqrt(s2), log = TRUE))

normal.ll

## likelihood for Cauchy model
## pars = c(mu,sigma)
## Cauchir dist qual student-t with df=1 
nll.cauchy <- function(pars,x){
    -sum(dcauchy(x,location = pars[1],scale = pars[2],log=TRUE))
}

opt <- nlminb(c(median(x),2), nll.cauchy,lower=c(-Inf,0), x = x)
opt
mean(x)

## Compare Cauchy and normal by AIC
-2 * normal.ll + 4
2 * opt$objective + 4 ## Negative log-likelihood
##################################################

## Profile likelihood for df i t-distribution
lp.ny <- function(ny,x){
    fun.tmp <- function(pars,ny,x){
    - sum(-log(pars[2]) + 
            dt((x-pars[1])/pars[2], df=ny,
               log = TRUE))
    }
    -nlminb(c(median(x),2),fun.tmp,
            lower=c(-Inf,0),ny=ny,x=x)$objective
}

par(mfrow=c(1,1))
ny <- seq(0.4,3,by=0.01)
llp <- sapply(ny,lp.ny,x=x)
plot(ny,llp-max(llp),type="l")
lines(c(0,4),-qchisq(0.95,df=1)/2*c(1,1),col=4)
## hence n=1 is a reasonable choice
## Reporting the final model include reporting mu and sigma
## and their uncertainties (profile likelihood in the
## book)
##################################################

##################################################
## Example 6.2
##################################################
dat <- read.table("hanford.dat",skip=6,
                  header=FALSE)

## Look at data
plot(dat[ ,2],dat[ ,3])

## design matrix, and obs
X <- cbind(1,dat[ ,2])
y <- dat[ ,3]
n <- length(y)
X
y

## correlation between estimates
cov2cor(t(X)%*%X)

## remove the correlation
X[ ,2] <- X[ ,2] - mean(X[ ,2])
cov2cor(t(X)%*%X)


## Parameter estimaes
beta <- solve(t(X) %*% X) %*% t(X) %*% y
beta

## variance parameter
yhat <- X %*% beta
sigmasq.hat <- sum((y-yhat)^2)/n
## Unbiased estimate
sigmasq.hat2 <- sum((y-yhat)^2)/(n - 2)

sigmasq.hat
sigmasq.hat2

## se for beta
se.beta <- sqrt(diag(sigmasq.hat2 * solve(t(X) %*% X)))
se.beta

## Directly in R
summary(lm(y ~ X[ ,2]))
cbind(beta,se.beta)
sqrt(sigmasq.hat2)

## qqplot of residuals
qqnorm(y-yhat)
qqline(y-yhat)
## Note the effect of the small dataset...


##################################################
## Example 6.3
##################################################
dat <- read.table("surgery.dat")

plot(dat) ## Bernulli trail

## negative log likelihood function
nll <- function(beta,y,X,n){
    p <- exp(X %*% beta) / (1 + exp(X %*% beta))
    -sum(dbinom(y,size=n,prob=p,log=TRUE))
}

## Observation and design matrix
y <- dat[ ,2]
X <- cbind(1,dat[ ,1]-mean(dat[ ,1]))

opt <- nlminb(c(-1,1),nll,y=y,X=X,n=1)

glm(y ~ -1 + X, family = binomial)
opt$par


## Parameter uncertainty
library(numDeriv)
H <- hessian(nll, opt$par, y=y, X=X, n=1)
se.beta <- sqrt(diag(solve(H)))

summary(glm(y~-1+X,family=binomial))
se.beta

lines(dat[ ,1], exp(X %*% opt$par)/(1+exp(X %*% opt$par)))
## Interpreted as the probability of death as a
## function of age

##################################################
## Example 6.4
##################################################
dat <- read.table("seed.dat")

y<-dat[,1]  # success
n<-dat[,2]  # number of planted seeds

## See Table 6.4
seed<- c(rep(0,11),rep(1,10))
extract<- c(rep(0,5),rep(1,6),rep(0,5),rep(1,5))

cbind(y,n,seed,extract)

(fit <- glm(cbind(y,n-y)~factor(seed)*factor(extract),
           family=binomial))

model.matrix(fit)
summary(fit)

##################################################
## Example 6.7
##################################################
dat <- read.table("compete.dat")
plot(dat[ ,1],dat[ ,2],ylim=c(0,125))

## Define x and y
y <- dat[ ,2] 
x <- log(dat[ ,1]) - mean(log(dat[ ,1]))
n <- length(y)
plot(x,y)

## Estimation with glm
summary(glm(y ~ x, family = Gamma))

## Our own implementation
X <- cbind(1,x) ## design matrix
nll <- function(beta,phi,X,y){
    h <- X %*% beta
    mu <- 1/h
    theta <- -1/mu
    A <- -log(-theta)
    c <- sum((log(y)-log(phi))/phi -log(gamma(1/phi)))
    -(sum( (y * theta - A)/phi) + c)
}

phi <- 1
opt <- nlminb(c(1, 0), nll, phi = 1, X = X, y=y)
H <- hessian(nll,opt$par, phi = 1, X = X, y=y)

## Estimamte phi by method of moments
h <- X%*%opt$par
mu <- 1/h
phi.mm <- sum((y-mu)^2/mu^2)/(n-2)

se.beta <- sqrt(phi.mm * diag(solve(H)))

summary(fit <- glm(y ~ x, family = Gamma))
cbind(opt$par,se.beta)


## profile for phi
phi <- seq(0.002,0.1, by =0.0001)
lp <- numeric(length(phi))

for(i in 1:length(lp)){
    lp[i] <- nll(opt$par,phi[i],X,y)
}

plot(phi,-(lp+max(-lp)),type="l",ylim=c(-3,0))
lines(range(phi),-qchisq(0.95,df=1)/2*c(1,1),
      lty=2,col=4)
phi.mm

## The fit 
plot(dat[ ,1],dat[ ,2],ylim=c(0,125))
lines(dat[ ,1],1/(model.matrix(fit) %*% coef(fit)))

################################################
## End
################################################
