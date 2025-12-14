# setwd("~/undervisning/Kurser/02418/2025Fall/lect09")
rm(list=ls())
## Convergence of the average
N <- 10000
x1 <- rnorm(N)
m1 <- cumsum(x1)/(1:N)
plot(m1,type="l",ylim=c(-0.5,0.5))
lines(1:N,qnorm(0.975)/sqrt(1:N),col=2,lty=2)
lines(1:N,-qnorm(0.975)/sqrt(1:N),col=2,lty=2)

## Convergence with another distribution
a <- -1; b<- 1
mu <- (b+a)/2
sigma <- sqrt((b-a)^2/12)
x1 <- runif(N,a,b)
m1 <- cumsum(x1-mu)/(1:N)
plot(1:N,m1,type="l",ylim=c(-0.5,0.5))
lines(1:N,qnorm(0.975)/sqrt(1:N),col=2,lty=2)
lines(1:N,-qnorm(0.975)/sqrt(1:N),col=2,lty=2)

##################################################
## CLT
x <- rbinom(1000,size=5,prob=0.5)
qqnorm(x)
qqline(x)

x <- rbinom(1000,size=25,prob=0.5)
qqnorm(x)
qqline(x)

x <- rbinom(1000,size=250,prob=0.5)
qqnorm(x)
qqline(x)

x <- rbinom(1000,size=1000,prob=0.5)
qqnorm(x)
qqline(x)


##################################################
## Score test

## Simple Poisson example

## log likelihood
ll <- function(lambda,y){
    n <- length(y)
    n*(mean(y) * log(lambda) - lambda)
}

## Score function
score <- function(lambda,y){
    n <- length(y)
    n * (mean(y)/lambda-1)
}

## Observed information
ObsInfo <- function(lambda,y){
    n <- length(y)
    n * mean(y)/lambda^2
}

## Expected information
ExpInfo <- function(lambda,y){
    n <- length(y)
    n  / lambda
}

## Simulation example
set.seed(324)
(y <- rpois(10,lambda=2))
lambda0 <- 1 ## Null hypothesis
lambda <- seq(0.75,3.5,length=1e3)
plot(lambda,sapply(lambda,ll,y=y)-ll(mean(y),y),type="l",lwd=2) ## log likelihood
lines(lambda0*c(1,1),range(sapply(lambda,ll,y=y)-ll(mean(y),y)))## Null hypothesis

## Score
s <- score(lambda0,y) ## Score at h0
## Draw the tangent
a <- ll(lambda0,y) - s * lambda0 - ll(mean(y),y)
lines(lambda,a+s*lambda,col=2,lwd=2)


## Score test
z.exp <- score(lambda0,y)/sqrt(ExpInfo(lambda0,y))
z.exp
(p.score.exp <- 1-pchisq(z.exp^2,df=1)) ## p-value

## Score test (using observed information)
library(numDeriv)
(I <- -hessian(ll,lambda0,y=y))
ObsInfo(lambda0,y)
(z.obs <- score(lambda0,y)/sqrt(I))
(p.score.obs <- 1-pchisq(z.obs^2,df=1))## P-value
## As we see these are quite different in the present case

p.values1 <- c(score.obs=p.score.obs)


## Wald test:
## Quadratic approximation
Info <- ObsInfo(mean(y),y)

lines(lambda, - 0.5* Info * (lambda-mean(y))^2,lty=2)

## Wald test statistic
W <- (mean(y)-lambda0)^2/(1/Info) 
lines(range(lambda), -W/2*c(1,1),col=4,lty=2) ## log-likelihood at H_0
lines(mean(y)*c(1,1), -W/2*c(0,1)) ## Test stat.

## p-value
(p.wald <- 1-pchisq(W,df=1))
p.values1 <- c(p.values1,p.wald=p.wald)



## Confidence interval
((CI <- mean(y)+c(-1,1)*qnorm(0.975)/sqrt(Info)))

## LRT
lines(range(lambda),(ll(lambda0,y)-ll(mean(y),y))*c(1,1),lty=2,col="green",lwd=2)


## LRT
Q <- 2*(ll(mean(y),y)-ll(lambda0,y))
(p.LRT <- 1-pchisq(Q,df=1))
lines(mean(y)*c(1,1), c(0,-Q/2),lwd=2)
(p.values1 <- c(p.values1,p.LRT=p.LRT))


## likelihood CI (approximate)
Q <- 2*(ll(mean(y),y)-ll(lambda,y))
## abline(a=-qchisq(0.95,df=1)/2,b=0,col=4)
(CI <- range(lambda[Q<qchisq(0.95,df=1)]))

## CI illustration
lines(range(lambda),-qchisq(0.95,df=1)*c(1,1)/2)

CI.wald <- mean(y)+c(-1,1)*qnorm(0.975)/sqrt(Info) ## Wald CI

points(CI,-qchisq(0.95,df=1)*c(1,1)/2,col="green",pch=19) ## Profile likelihood CI
points(CI.wald,-qchisq(0.95,df=1)*c(1,1)/2,col="blue",pch=19) ## Wald CI


##################################################
## Transformed parameters, look at log(mu)
## x11()

## Likelihood
ll <- function(theta,y){
    n <- length(y)
    n*(mean(y) * theta - exp(theta))
}

## score function
score <- function(theta,y){
    n <- length(y)
    n * (mean(y)-exp(theta))
}

## Observed information
ObsInfo <- function(theta,y){
    n <- length(y)
    n * exp(theta)
}

## Expecteed information
ExpInfo <- function(theta,y){
    n <- length(y)
    n  * exp(theta)
}


theta0 <- 0 ## Null hypothesis (corresponding to lambda=1 above)

theta <- seq(log(0.75),log(3.5),by=0.01)
plot(theta,sapply(theta,ll,y=y)-ll(log(mean(y)),y=y),type="l",lwd=2)
lines(theta0*c(1,1),range(sapply(theta,ll,y=y)-ll(log(mean(y)),y)))

## Score
s <- score(theta0,y)
a <- ll(theta0,y) - s * theta0 -ll(log(mean(y)),y)
lines(theta,a+s*theta,col=2,lwd=2)


## Score test
z.exp2 <- score(theta0,y)/sqrt(ExpInfo(theta0,y))
z.exp2
(p.score.exp2 <- 1-pchisq(z.exp2^2,df=1)) ## p-value



## Score test (using observed information)
library(numDeriv)
(I <- -hessian(ll,theta0,y=y))
ObsInfo(theta0,y)
(z.obs2 <- score(theta0,y)/sqrt(I))
(p.score.obs2 <- 1-pchisq(z.obs2^2,df=1))## P-value


## Quadratic approximation
Info <- ObsInfo(log(mean(y)),y)
lines(theta,  - 0.5* Info * (theta-log(mean(y)))^2,lty=2)


## Wald
W <- (log(mean(y))-theta0)^2/(1/Info)
lines(range(theta), (-W/2)*c(1,1),col="blue",lwd=2,lty=2)
lines(log(mean(y))*c(1,1), c(0,-W/2))

(p.wald2 <- 1-pchisq(W,df=1))
p.values2 <- c(p.score.obs=p.score.obs2,p.wald2=p.wald2)


## LRT
lines(range(theta),ll(theta0,y)*c(1,1)-ll(log(mean(y)),y),lty=2,col="green",lwd=2)

Q2 <- 2*(ll(log(mean(y)),y)-ll(theta0,y))
(p.LRT2 <- 1-pchisq(Q2,df=1))
lines(log(mean(y))*c(1,1), c(0,-Q2/2),lwd=2)

p.values2 <- c(p.values2,p.LRT2=p.LRT2)

cbind(p.values1[-4],p.values2)
## We see better agreement between p.values 
## using the transformed parameters.

##################################################
## Example - random effect model
setwd("~/undervisning/Kurser/02418/2025Fall/lect09/")
## read data
y<- matrix(scan('random.dat'),ncol=5,byrow=T)
y <- as.vector(y)
x <- factor(rep(1:5,each=16))
y <-10*log10(y)
dat <- data.frame(y=y,x=x)
dat
## Plot data
plot(x,y)

## Usual analysis of variance
anova(lm(y~x))
## So difference between the women 


## Full likelihood
nll <- function(theta,y,group){
    S <- diag(rep(1,16))*theta[2]+
      matrix(1,ncol=16,nrow=16)*theta[3]
    nll <- 0
    for(i in 1:5){
        nll <- nll + 1/2 * log(det(S)) + 
          1/2 * t(y[group==i]-theta[1]) %*%
            solve(S) %*% (y[group==i]-theta[1])
    }
    nll
}

## check that it work
nll(c(mean(y),0.5,1.3),y,x)

## Find MLE
opt <- nlminb(c(mean(y),1,1),nll,y=y,group=x,lower=c(-Inf,0.1,0.1))

## Find parameter variance
V <- solve(hessian(nll,opt$par, y=y, group=x))
## Standard errors
se <- sqrt(diag(V))

## Some summary stat
sums <- cbind(opt$par,se,z.stat=opt$par/se)
rownames(sums) <- c("mu","sigma","sigma_a")
sums
## Hence sigma_a is not significant on a 5% level

## Profile for sigma_a
pnll.fun <- function(sigma.a, y, group){
    fun.tmp <- function(theta,sigma.a,y,group){
        nll(c(theta,sigma.a),y,group)
    }
    nlminb(c(mean(y),1),fun.tmp,sigma.a=sigma.a,y=y,group=x,
           lower=c(-Inf,0.1))$objective
}

## Choose vales of sigma_a for profiling
sigma.a <- seq(0,8,by=0.05)

## profile likeihood
pnll <- sapply(sigma.a,pnll.fun,y=y,group=x)

## Plot the profile likelihood
plot(sigma.a,-(pnll-min(pnll)),type="l",ylim=c(-10,0))
plot(sigma.a,exp(-(pnll-min(pnll))),type="l",
     ylim=c(0,1))
## 95% cut off
lines(range(sigma.a),c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2,lty=2)
## Quadractic approximation
lines(sigma.a,exp(-1/(2*V[3,3])*(sigma.a-opt$par[3])^2),
      lty=2)


## likelihood ratio test
(chisq <- 2*(-opt$objective-
               (-pnll.fun(0, y, group=x))))
pchisq(chisq,df=1,lower.tail=FALSE)
## So highly significant


## Profile likelihood for intraclass correlation
nll <- function(theta,y,group){
    S <- matrix(1,ncol=16,nrow=16) * theta[3] ##S_ij=rho
    diag(S) <- 1  ## S_ii=1
    S <- S * theta[2]
    nll <- 0
    for(i in 1:5){
        nll <- nll + 1/2 * log(det(S)) + 1/2 * t(y[group==i]-theta[1]) %*%
            solve(S) %*% (y[group==i]-theta[1])
    }
    nll
}


## check
nll(c(mean(y),0.5,0.5),y,x)

## Find MLE
opt <- nlminb(c(mean(y),1,1),nll,y=y,group=x,lower=c(-Inf,0.1,0),
              upper = c(Inf,Inf,0.99))

## Find se for parameters
V <- solve(hessian(nll,opt$par, y=y, group=x))
se <- sqrt(diag(V))

sums <- cbind(opt$par,se,z.stat=opt$par/se)
rownames(sums) <- c("mu","sigma+sigma_a","rho")
sums
## rho different from zero, but 95% CI for rho is
opt$par[3]+c(-1,1)*qnorm(0.975)*se[3]
## which include 1!


##################################################
## Profile for rho
pnll.fun <- function(sigma.a, y, group){
    fun.tmp <- function(theta,sigma.a,y,group){
        nll(c(theta,sigma.a),y,group)
    }
    nlminb(c(mean(y),1),fun.tmp,sigma.a=sigma.a,y=y,
           group=x,lower=c(-Inf,0.1))$objective
}


rho <- seq(0,0.99,by=0.01)

## profile likelihood
pnll <- sapply(rho,pnll.fun,y=y,group=x)

## plot the profile likeihood
plot(rho,-(pnll-min(pnll)),type="l",ylim=c(-10,0))
plot(rho,exp(-(pnll-min(pnll))),type="l",
     ylim=c(0,1))
## 95% cut off
lines(range(rho),c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2,lty=2)
## quadratic transformation
lines(rho,exp(-1/(2*V[3,3])*(rho-opt$par[3])^2),
      lty=2)

#####################################################
## End 
#####################################################
