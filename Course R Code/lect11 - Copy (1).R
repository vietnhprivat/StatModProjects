rm(list=ls())
setwd("~/undervisning/Kurser/02418/2025Fall/lect11/")
##################################################

## Data (including year 2007-2016)
dat <- read.table("earthquakes.txt",header=FALSE)
names(dat) <- c("year","eq")


## Plot data
par(mfrow=c(1,2))
tab <- table(dat$eq)
plot(tab,xlim=c(0,45),axes=FALSE,xlab="Count",ylab="Freq")
points(0:45,dpois(0:45,lambda=mean(dat$eq))*dim(dat)[1],pch=19)
mean(dat$eq)
var(dat$eq)
## hence overdispersion
acf(dat$eq)
## hence autocorreltaion

par(mfrow=c(1,1))
plot(dat$year,dat$eq,type="l")

##################################################
## Likelihood
##################################################

## book scripts
source("A1.R")
y <- dat$eq

## 1 - state 
## Initial values
m <- 1
lambda0 <- mean(y)
gamma0 <- 0

## optimize
fit1 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit1$AIC

## 2 - state 
## Initial values
m <- 2
lambda0 <- quantile(y,c(0.25,0.75))
gamma0 <- matrix(0.05,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## optimize
fit2 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit2$AIC

## 3 - state 
## Initial values
m <- 3
lambda0 <- quantile(y,c(0.25,0.5,0.75))
gamma0 <- matrix(0.05,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]
gamma0

## optimize
fit3 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit3$AIC


## 4 - state 
## Initial values
m <- 4
lambda0 <- quantile(y,c(0.2,0.4,0.6,0.8))
gamma0 <- matrix(0.025,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## optimize
fit4 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit4$AIC

AIC <- c(fit1$AIC,fit2$AIC,fit3$AIC,fit4$AIC)
ll <-  -c(fit1$mllk,fit2$mllk,fit3$mllk,fit4$mllk)
AIC
m <- c(1,2,3,4)
df <- m + (m^2-m) ## lambda + gamma
## What should we report
cbind(df, ll, AIC)
fit3

## AIC.pois.mix (mixture distribution)
##  841.8676 784.8275 782.1411 785.9570

##################################################
## working with the 3-state mmodel

## Finding the standard errors
m <- 3
## working parameters
parvect  <- pois.HMM.pn2pw(m,fit3$lambda,fit3$gamma)
## Optimize (hessian = TRUE return hessian)
mod <- nlm(pois.HMM.mllk,parvect,x=y,m=m,
            hessian=TRUE)  
mod

## Organize the result
parvect <- mod$estimate
names(parvect) <- c("lambda1","lambda2","lambda3","tau21",
                    "tau31","tau12","tau32","tau13","tau23")

se <- sqrt(diag(solve(mod$hessian)))

## Working pars + standard error
round(cbind(parvect,se),digits=2) ## note se of tau31
fit3$gamma

##################################################
## Parametric bootstrap from here

source("A2.R")
## Initailize
n <- length(y)
k <- 100
m <- 3
lambda0 <- quantile(y,c(0.25,0.5,0.75))
gamma0 <- matrix(0.025,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## Matrices to stores bootstap res.
GAMMA <- matrix(ncol=m*m,nrow=k)
Lambda <- matrix(ncol=m,nrow=k)
Delta <- matrix(ncol=m,nrow=k)
Code <- numeric(k)

## Parametric bootstrap using nlminb 
for(i in 1:k){
    set.seed(i)
    y.sim <- pois.HMM.generate_sample(n, m,
                                      fit3$lambda, fit3$gamma)
    lambda0 <- quantile(y.sim,c(0.25,0.5,0.75))
    fit3.tmp <- pois.HMM.mle.nlminb(y.sim, m,lambda0, gamma0)
    GAMMA[i, ] <- c(fit3.tmp$gamma[1, ],
                    fit3.tmp$gamma[2, ],
                    fit3.tmp$gamma[3, ])
    Lambda[i, ] <- fit3.tmp$lambda
    Delta[i, ] <- fit3.tmp$delta
    Code[i] <- fit3.tmp$code
    print(c(i,Code[i]))
 }
sum(Code!=0)

## Plot the results (i.e. statistical distribution of
## estimates)
par(mfrow=c(1,3))
hist(Lambda[ ,1],xlim=range(Lambda))
rug(fit3$lambda,lwd=2,col=2)
hist(Lambda[ ,2],xlim=range(Lambda))
rug(fit3$lambda,lwd=2,col=2)
hist(Lambda[ ,3],xlim=range(Lambda))
rug(fit3$lambda,lwd=2,col=2)
## There is a problem here (can you fix it?)

##################################################
## A bootstrap with 10.000 realisations
## Sim <- list(GAMMA=GAMMA,Lambda=Lambda,Delta=Delta,Dode=Code)
## save(Sim,file="sim_k_10000.Rdata")
load(file="sim_k_10000.Rdata")
par(mfrow=c(1,3))
lrange <- quantile(Sim$Lambda,prob=c(0,0.999))
hist(Sim$Lambda[ ,1],xlim=lrange,
     breaks=1:round(max(Sim$Lambda))+1)
rug(fit3$lambda,lwd=2,col=2)
hist(Sim$Lambda[ ,2],xlim=lrange,
     breaks=1:round(max(Sim$Lambda))+1)
rug(fit3$lambda,lwd=2,col=2)
hist(Sim$Lambda[ ,3],xlim=lrange,
     breaks=1:round(max(Sim$Lambda)+1))
rug(fit3$lambda,lwd=2,col=2)

## what is going on?
sum(Sim$Lambda[ ,1]>Sim$Lambda[ ,2])
sum(Sim$Lambda[ ,1]>Sim$Lambda[ ,3])
sum(Sim$Lambda[ ,2]>Sim$Lambda[ ,3])
## ie not sorted (you should fix that..)
##################################################

## Distribution of transition probability matrices
par(mfrow=c(3,3))
for(i in 1:9){
  hist(GAMMA[ ,i],xlim=c(0,1))
  rug(as.vector(t(fit3$gamma))[i],lwd=2,col=2)
}
## You should fix the sorting problem for this to make 
## absolute sense

##################################################
## Confidence intervals (90%)
## 95% CI for lambda 
apply(Lambda,2,quantile,prob=c(0.025,0.975))
## 95% CI for gamma
t(round(apply(GAMMA,2,quantile,prob=c(0.025,0.975)),
      digits=3))

## 95% CI for delta
round(t(apply(Delta,2,quantile,prob=c(0.025,0.975))),digits=3)

##################################################
## Profile likelihood for lambda 1
PL.lambda1 <- function(lambda1,m,y,lambda0,gamma0){
    ## Fun for inner optim
    fun.tmp <- function(pars,lambda1,y,m){
        parvect <- c(log(lambda1),pars)
        pois.HMM.mllk(parvect,y,m)
    }
    ## Initialize
    lambda0 <- c(lambda1,lambda0)
    parvect0 <- pois.HMM.pn2pw(m, lambda0, gamma0)
    parvect0 <- parvect0[-1]
    np <- length(parvect0)
    lower    <- rep(-10,np)
    upper    <- c(rep(max(y),m-1),rep(10,np+1-m))
    ## optimize to find profile likelihood
    nlminb(parvect0,fun.tmp, lambda1=lambda1,
           y=y, m=m, lower=lower,
           upper=upper)$objective    
}

## Initial values for estimation
lambda0 <- quantile(y,probs=c(1/3,2/3))
PL.lambda1(lambda1=29,m=m,y=y,lambda0,gamma0)

## Which lamdas should we look at
lambda1 <- seq(min(y),max(y),length=100)

## The profile liklielihood 
llp.lambda1 <- sapply(lambda1,PL.lambda1,m=m,y=y,
       lambda0=lambda0,gamma0=gamma0)

## Plot the profile likelihood
par(mfrow=c(1,1))
plot(lambda1,exp(-(llp.lambda1-fit3$mllk)),
     type="l")
lines(range(lambda1),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),col=2,lty=2,lwd=2)
rug(fit3$lambda,col=3,lwd=2)
## Looks bad!!!


## A fix?
PL.lambda1 <- function(lambda1,m,y,lambda0,gamma0){
  ## Fun for inner optim
  fun.tmp <- function(pars,lambda1,y,m){
    parvect <- c(log(lambda1),pars)
    pois.HMM.mllk(parvect,y,m)
  }
  ## Initialize
  lambda0 <- c(lambda1,lambda0)
  parvect0 <- pois.HMM.pn2pw(m, lambda0, gamma0)
  parvect0 <- parvect0[-1]
  np <- length(parvect0)
  lower    <- rep(-10,np)
  upper    <- c(rep(max(y),m-1),rep(10,np+1-m))
  ## optimize to find profile likelihood
  opt <- nlminb(parvect0,fun.tmp, lambda1=lambda1,
                y=y, m=m, lower=lower,
                upper=upper)
  opt ## Only fifference!! Return result of optimization, used for warm start.
}

## Using previous result for warm start...
lambda1 <- seq(min(y),max(y),length=100)
lambda0 <- quantile(y,probs=c(1/3,2/3))
gamma0 <- matrix(0.05,ncol=m,nrow=m)
llp.lambda1 <- numeric(length(lambda1))

for(i in 1:length(lambda1)){
  print(i)
  tmp <- PL.lambda1(lambda1[i],m,y,lambda0,gamma0)
  par0 <- pois.HMM.pw2pn(3,c(log(lambda1[i]),tmp$par*0.9))
  lambda0 <- par0$lambda[-1]
  gamma0 <- par0$gamma
  llp.lambda1[i]<-tmp$objective
}


## Plot the profile likelihood
par(mfrow=c(1,1))
plot(lambda1,exp(-(llp.lambda1-fit3$mllk)),
     type="l")
lines(range(lambda1),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),col=2,lty=2,lwd=2)
rug(fit3$lambda,col=3,lwd=2)
## Look better


## Wald statistic
cbind(mod$estimate,se)

## Quadratic (local) approximation (Wald CI)
cbind(exp(mod$estimate-1.96*se),
      exp(mod$estimate+1.96*se))[1:3, ]
lines(exp(mod$estimate[1]-1.96*se[1]*c(-1,1)),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2)
lines(exp(mod$estimate[2]-1.96*se[2]*c(-1,1)),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2)
lines(exp(mod$estimate[3]-1.96*se[3]*c(-1,1)),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2)
##################################################

## Why do we have 3 optimal values?

##################################################
## The end
##################################################