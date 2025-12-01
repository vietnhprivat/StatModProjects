setwd("~/Kurser/02418/2021Fall/lect08")
###################################################
## Score functions
##################################################

## Normal model

## likelihood
ll <- function(theta,sigma,x){
    sum(dnorm(x,mean=theta,sd=sigma))
}

## score function
s.norm <- function(theta,sigma,x){
    n <- length(x)
    n/sigma^2 * (mean(x) - theta)
}

par(mfrow=c(2,2))
## parameters
theta0 <- 0; sigma <- 1; n <- 160
## random sample
x <- rnorm(n,mean=theta0,sd=sigma)

## Plot the score function
theta <- seq(-1,1,by=0.1)
plot(theta,sapply(theta,s.norm,sigma=1,x=x),type="l",
     ylab="score",
     xlab=expression(theta),ylim=c(-25,25))
abline(a=0,b=0,col=gray(0.5))
abline(v=theta0,col=gray(0.5))

## Plot the score function for another realisation
x <- rnorm(n,mean=theta0,sd=sigma)
lines(theta,sapply(theta,s.norm,sigma=1,x=x))

## Plot for even more realisations
for(i in 1:10){
    x <- rnorm(n,mean=theta0,sd=sigma)
    lines(theta,sapply(theta,s.norm,sigma=1,x=x))
}

##################################################
## Poisson model

## Score function
s.pois <- function(theta,x){
    n <- length(x)
    n/theta * (mean(x) - theta)
}

par(mfrow=c(1,1))

## Parameters
theta0 <- 4; n <- 20

## random sample
x <- rpois(n,lambda=theta0)

## Plot the score function
theta <- seq(2,8,by=0.1)
plot(theta,sapply(theta,s.pois,x=x),type="l",
     ylab="score",xlab=expression(theta),ylim=c(-10,15))
abline(a=0,b=0,col=gray(0.5))
abline(v=theta0,col=gray(0.5))

## Plot the score function for another realisation
x <- rpois(n,lambda=theta0)
lines(theta,sapply(theta,s.pois,x=x))

## Plot for even more realisations
for(i in 1:10){
    x <- rpois(n,lambda=theta0)
    lines(theta,sapply(theta,s.pois,x=x))
}


##################################################
## Binomial model

## score function
s.binom <- function(theta,x,N){
    n <- length(x)
    n *(mean(x) - N * theta) /(theta*(1-theta))
}


par(mfrow=c(1,1))

## Parameters
theta0 <- 0.4; N <- 10; n <- 10

## Random sample
x <- rbinom(n,size = N,prob=theta0)

## Plot the score function
theta <- seq(0.1,0.7,by=0.01)
plot(theta,sapply(theta,s.binom,x=x,N=N),type="l",
     ylab="score",xlab=expression(theta),
     ylim=c(-300,500))
abline(a=0,b=0,col=gray(0.5))
abline(v=theta0,col=gray(0.5))

## Plot the score function for another realisation
x <- rbinom(n,size = N,prob=theta0)
lines(theta,sapply(theta,s.binom,x=x,N=N))

## Plot for even more realisations
for(i in 1:10){
    x <- rbinom(n,size = N,prob=theta0)
    lines(theta,sapply(theta,s.binom,x=x,N=N))
}

##################################################
## Cauchy model

## negative log likelihood
lL.cau <- function(theta,x){
    - sum(log(pi*(1+(x-theta)^2)))
}

## Score function
s.cau <- function(theta,x){
    sum(2*(x-theta)/(1+(x-theta)^2))
}

par(mfrow=c(1,1))

## Set parameters
theta0 <- 0; n <- 10

## random sample
x <- rcauchy(n,location = theta0)

## Plot the score function
theta <- seq(-10,10,by=0.01)
plot(theta,sapply(theta,s.cau,x=x),type="l",
     ylab="score",xlab=expression(theta),
     ylim=c(-15,15))
abline(a=0,b=0,col=gray(0.5))
abline(v=theta0,col=gray(0.5))

## Plot for more realisations
for(i in 1:19){
    x <- rcauchy(n,location = theta0)
    lines(theta,sapply(theta,s.cau,x=x))
}
## example of very complicated score function
##################################################



##################################################
## Survival regression model in R
xdat<- scan('rat.dat',skip=1,
            what=list(group=0,surv=0,status=0))
xdat


library(survival)
dat <- data.frame(t=xdat$surv,dead=xdat$status,
                  group=xdat$group)
(mod1 <- survreg(Surv(t, dead) ~ 1, data = dat,
                 dist = "exponential"))
summary(mod1)

sqrt(1/sum(dat$d)) ## SE of parameter

##################################################
##
#################################################

## Logistic regression model
N <- 10; n <- 20; theta0 <- 0.2
x <- rbinom(n,prob=theta0,size=N)
fit <- glm(cbind(x,N-x)~1,family=binomial)
summary(fit)

psi <- coef(fit)
I <- n*N*exp(psi)/(1+exp(psi))^2
sqrt(1/I)


##################################################
## A poisson simple example

## Simulate some data
n1 <- 10; n2 <- 20; N <- n1+n2
y1 <- rpois(n1,lambda=1)
y2 <- rpois(n2,lambda=2)
y <- c(y1,y2)
x <- c(rep(0,n1),rep(1,n2))

## fit the model using glm
fit1 <- glm(y~x,family=poisson)
fit1

## Information matrix
solve(summary(fit1)$cov.unscaled)
n1*mean(y1)+n2*mean(y2)
n2*mean(y2)

## The resulting covariance estimate
summary(fit1)$cov.unscaled

## diagonal elements
c(1/(n1*mean(y1)),(n1+mean(y2)/mean(y1)*n2)/(n1*n2*mean(y2)))



## The parameters
coef(fit1)
log(mean(y1))
log(mean(y2))-log(mean(y1))

## another link function
fit2 <- glm(y~x,family=poisson(link="identity"))
fit2
mean(y1)
mean(y2-y1)

solve(summary(fit2)$cov.unscaled)
n1/mean(y1)+n2/mean(y2)
n2/mean(y2)



##################################################
## A GLM example
smoke <- read.table("Data/smoke.csv",sep=",",
            header=TRUE)

age.cat <- 1:5
smoke <- cbind(age.cat,smoke)

## Data
smoke

par(mfrow=c(1,1))
## plot data
plot(smoke$age.cat,smoke$death/smoke$person.years,
     pch=rep(1:2,each=5),col=rep(1:2,each=5))

## an initial model
fit1 <- glm(deaths ~ smoking + age.cat +
                offset(log(person.years)),
            data = smoke, family = poisson)

summary(fit1)
par(mfrow=c(2,2))
plot(fit1)

fit2 <- glm(deaths ~ smoking + smoking:age.cat + age.cat +
                offset(log(person.years)), data=smoke,
            family=poisson)
summary(fit2)

par(mfrow=c(2,2))
plot(fit2)


anova(fit1,fit2,test="Chisq")

fit3 <- glm(deaths ~ smoking + smoking:age.cat + age.cat +
                I(age.cat^2) +
                offset(log(person.years)), data=smoke,
            family=poisson)
summary(fit3)

plot(fit3)

anova(fit1,fit2,fit3,test="Chisq")## Plot the result
par(mfrow=c(1,1))
plot(smoke$age.cat,smoke$death,
     pch=rep(3,each=5),col=rep(1:2,each=5))
points(smoke$age.cat,predict(fit3,type="response"),
       col=rep(1:2,each=5),pch=19)

par(mfrow=c(1,1))
plot(smoke$age.cat,smoke$death/smoke$person.years,
     pch=rep(3,each=5),col=rep(1:2,each=5))
points(smoke$age.cat,predict(fit3,type="response")/
                     smoke$person.years,
       col=rep(1:2,each=5),pch=19)
##################################################
## End
##################################################