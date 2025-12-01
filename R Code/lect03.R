n=7
y<- 5
nl.fun<- function(p,y,n){
  -dbinom(y,size=n,prob=p,log=TRUE)
}

p<-seq(0,1,by=0.01)
plot(p,-nl.fun(p,y,n),type="l",ylim=c(-5,0))


#######################################
rm(list=ls())
#######################################
## Example 4.2
## no of boys in families with 12 children

## Data
n <- 12
k<- 0:n
nk<- c(3, 24, 104, 286, 670, 1033, 1343,
      1112 , 829, 478, 181, 45, 7)

## Total no of families
(N<- sum(nk))

## Proportion in each cell
pk <- nk/N
round(pk,digits=4)

## Estimates (including se)
(theta <- sum(k*pk)/12)
(Ihat <-  N * n /(theta * (1 - theta)))
(se <- sqrt(1/Ihat))

## Observed variance
(m <- sum(k * nk) / N) ## obs mean
(sum((k - m)^2 * nk / (N - 1))) ## obs variance

## expected variance from model
n * theta * (1 - theta)
## Q: Is 3.48 significantly different from 3?

## A: look at the goodness of fit test
## residuals:
est <- N * dbinom(k,12,theta)     # expected frequencies
round(est,digits=2)

## Collapsing cells (to get e_(ij)>5)
est[2] <- sum(est[1:2])
est[12] <- sum(est[12:13])
(est <- est[-c(1,13)])

obs <-nk
obs[2] <- sum(obs[1:2])
obs[12] <- sum(obs[12:13])
(obs <- obs[-c(1,13)])
rbind(obs,round(est,digits=2))


## Goodness of fit test
(chi2 <- sum((obs-est)^2/est))

## "residuals"
(obs-est)/sqrt(est)

pchisq(chi2,df=length(obs)-1,lower.tail=FALSE)
## Hence assumption is not met (i.e overdispersion)


##################################################
## Comparing proportions
##################################################

## Die ex. (today)
n1 <- 7
x1 <- 5

## Large dice exp.
x2 <- 252 - 58
n2 <- 252

## chi² test directly in R
prop.test(c(x1,x2), c(n1, n2))


## Log odds ratio
(theta <- log(x1/(n1-x1)*(n2-x2)/x2))


## Plot
loglik.fun <- function(theta, eta, x, n){
    a <- exp(theta + eta)
    b <- exp(eta)
    x[1] * theta + (x[1] + x[2]) * eta - n[1] *
        log(1 + a) - n[2] * log(1 + b)
}


eta <- seq(0.75,1.75,length=100)
theta <- seq(-2,5,length=100)


par(mfrow=c(1,2))
ll<- outer(theta,eta,'loglik.fun', x=c(x1,x2), n=c(n1,n2))
like<- exp(ll-max(ll))
contour(theta,eta,like,level=c(0.05,.1,.3,.5,.7,.9),
        xlab=expression(theta),
        ylab=expression(eta))
  title(expression('Likelihood contour'))
lines(c(0,0),c(0,2),lty=2,col=2,lwd=2)

## profile likelihood
lp.theta <- function(theta,x,n){
    fun.tmp <- function(eta,theta,x,n){
        loglik.fun(theta, eta, x, n)
    }
    ## interval from plot
    optimise(fun.tmp,c(-20,20),theta=theta,
             x=x, n=n,maximum=TRUE)$objective
}


lp <- sapply(theta,lp.theta,x=c(x1,x2),n=c(n1,n2))

lp <- lp - max(lp)

plot(theta,exp(lp),type="l")
lines(c(-2,10),exp(-qchisq(0.95,df=1)/2)*c(1,1), col=2)
lines(c(0,0),c(0,1),col=4,lty=2,lwd=2)
title(expression('Profile Likelihood'))
 
## So we cannot reject that p1=p2

##################################################
## Poisson models
##################################################

## See EX4-7.R plus EX4-7C.R

##################################################

##################################################
## Revisit binomial example
## Data
n <- 12
k<- 0:n
nk<- c(3, 24, 104, 286, 670, 1033, 1343,
      1112 , 829, 478, 181, 45, 7)
## Total no of families
(N<- sum(nk))
## Proportion in each cell
pk <- nk/N
round(pk,digits=4)

## Estimates (including se)
(theta <- sum(k*pk)/12)
(Ihat <-  N * n /(theta * (1 - theta)))
(se <- sqrt(1/Ihat))


## Likelihood function, exponential dispersion fam.
ll <- function(theta,nk,n, phi){
    sum(n / phi * ((0:n)/n * log(theta/(1-theta)) * nk +
                   nk * log(1 - theta) ))
}

## Ml estimate
(opt <- optimize(ll,c(0,1),nk=nk,n=n,maximum=TRUE,phi=1))

##################################################
## Plot likelihood
d <- 0.005
theta <- seq(opt$maximum-d,opt$maximum+d,by=0.0001)
logL<-sapply(theta,ll,nk=nk,n=n,phi=1)
logL <- logL-max(logL)
plot(theta,exp(logL),type="l")


## Observed mean and variance
(m <- sum(k * nk) / N)
(s2 <- (sum((k-m)^2*nk/N)))

## expected variance from model
(v <- n*opt$maximum*(1-opt$maximum))

s2/v


logL<-sapply(theta,ll,nk=nk,n=n,phi=s2/v)## with overdispersion
logL <- logL-max(logL)
lines(theta,exp(logL),col=2,lwd=2)


## Using generalized linear model in R
z <- rep(0:12,nk)
resp <- cbind(z,12-z)

logis.glm<-glm(resp~1,family=binomial)
summary(logis.glm)
## Parameter in natural domain
1/(1+exp(-coef(logis.glm)))
opt$maximum



##################################################
## Hypothesis test of H_0:phi=1 (see ex. 4.14)
logis.glm2 <- glm(resp~1,family=quasibinomial)
summary(logis.glm2)
s2/v ## Is this different from 1?

phi.sim <- c()
for(i in 1:1000){
    r<-rbinom(N,size=12,prob=opt$maximum)
    p <- sum(r)/(N*12)
    phi.sim[i] <- var(r)/(12*p*(1-p))
}

##############################################
## Is phi != 0?
sum(phi.sim>s2/v)
quantile(phi.sim,probs=c(0.9,0.95,0.99,1))
## So yes phi !=1
##################################################
##
##################################################
