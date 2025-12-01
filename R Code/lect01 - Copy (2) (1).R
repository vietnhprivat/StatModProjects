##################################
## remove all variables from work space
rm(list=ls())

##################################
## Likelihood function
L <- function(theta, y, n){
  theta^y * (1 - theta)^(n - y)
}

## Data fill in todays data
n <- 7
y <- 5

## Plot likelihood pf the parameter
theta<- seq(0,1,by=0.001)
plot(theta,L(theta,y,n)/max(L(theta,y,n)),type="l")
alpha <-0.05
lines(c(0,1),exp(-0.5*qchisq(1-alpha,df=1))*c(1,1),
      col=2,lty=2)

## Likelihood ratio
L(y/n,y,n)/L(1/6,y,n)

## What if we had 10 times the number of obs?
n2 <- n * 10
y2 <- y * 10
lines(theta,L(theta,y2,n2)/max(L(theta,y2,n2)),col=3)

## Log - likelihhod
ll <- function(theta,y,n){
  log(theta^y*(1-theta)^(n-y))
}

## Plot log-likelihood
theta<- seq(0.4, 1, by = 0.001)
plot(theta,ll(theta,y,n)-max(ll(theta,y,n)),type="l")
lines(theta,ll(theta,y2,n2)-max(ll(theta,y2,n2)),col=3)


############################################
## Optimize
?optimise

## Find optimal parameters by numerical optimization
(opt1 <- optimise(ll,c(0,1),y=y,n=n,maximum=TRUE))
y/n

## Another optimizer
?nlminb


## Invariance
## Log - likelihhod
ll2 <- function(theta,y,n){
    p <- exp(theta)/(1+exp(theta))
    log(p^y*(1-p)^(n-y))
}

## Optimize likelihood
(opt2 <- optimise(ll2,c(-10,10),y=y,n=n,maximum=TRUE))
## Notice same likelihood
(theta <- opt2$maximum)
## and same parameter
exp(theta)/(1+exp(theta))

ylim <- c(-10,0)
## Plot log-likelihood (original parameter domain)
par(mfrow=c(1,2))
theta<- seq(0, 1, by = 0.001) 
plot(theta,ll(theta,y,n)-max(ll(theta,y,n)),type="l",ylim=ylim)
lines(theta,ll(theta,y2,n2)-max(ll(theta,y2,n2)),col=3)


## Plot log-likelihood (transformed parameter domain)
theta<- seq(-2, 6, by = 0.001)
plot(theta,ll2(theta,y,n)-max(ll2(theta,y,n)),type="l",ylim=ylim)
lines(theta,ll2(theta,y2,n2)-max(ll2(theta,y2,n2)),col=3)
## Likelihood more regular in transformed domain.

##################################################
## Quadratic approximation

## Find hessian numarically
library(numDeriv)
H1 <- as.numeric(hessian(ll,opt1$maximum,y=y,n=n))
V1 <- -1/H1
V1

## What we would expect
opt1$maximum*(1-opt1$maximum)/n


## Plot log-likelihood (original parameter domain)
par(mfrow=c(1,2));ylim=c(-4,0)
theta<- seq(0.3, 0.99, by = 0.001) 
plot(theta,ll(theta,y,n)-max(ll(theta,y,n)),type="l",ylim=ylim)
## Quadratic approximation
lines(theta,0.5 * H1 * (theta -opt1$maximum)^2,lty=2)
## In this case quite bad (due to the small sample

## Likelihood based CI
alpha <-0.05
lines(c(0,1),-0.5*qchisq(1-alpha,df=1)*c(1,1),col=2,lty=2)



## Plot log-likelihood (transformed parameter domain)
theta<- log(theta/(1-theta))
plot(theta,ll2(theta,y,n)-max(ll2(theta,y,n)),type="l",ylim=c(-4,0))

## Hessian (and variance estimate)
H2 <- as.numeric(hessian(ll2,opt2$maximum,y=y,n=n))
V2 <- -1/H2
V2
lines(theta,0.5 * H2 * (theta -opt2$maximum)^2,lty=2)

## Likelihood based CI
alpha <-0.05
lines(range(theta),-0.5*qchisq(1-alpha,df=1)*c(1,1),col=2,lty=2)

## Wald CI
## Original domain:
opt1$maximum + c(-1,1) * qnorm(1-alpha/2) * sqrt(V1)

## Transformed domain:
(CI2 <- opt2$maximum + c(-1,1) * qnorm(1-alpha/2) * sqrt(V2)) ## CI transformed domain

exp(CI2)/(1+exp(CI2)) ## Back to original domain

##################################################
## End
##################################################

