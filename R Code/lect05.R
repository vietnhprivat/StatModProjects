setwd("~/undervisning/Kurser/02418/2025Fall/lect05")
rm(list=ls())
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

## Our own implementation (directly form slides)
X <- cbind(1,x) ## design matrix
nll <- function(beta,phi,X,y){
    h <- X %*% beta ## linear model
    mu <- 1/h  ## mean value model
    theta <- -1/mu ## theta
    A <- -log(-theta) ## Kumulant generating function
    ## Normalization constant
    c <- sum((log(y)-log(phi))/phi -log(gamma(1/phi)))    
    -(sum( (y * theta - A)/phi) + c) ## log-Likelihood
}

## Estimate with phi =1 (including parameter covariance)
library(numDeriv)

phi <- 1
(opt <- nlminb(c(1, 0.1), nll, phi = 1, X = X, y=y))
H <- hessian(nll,opt$par, phi = 1, X = X, y=y)
## Estimamte phi by method of moments
h <- X %*% opt$par
mu <- 1/h
## Compare model and obs
lines(x,mu,type="l")

phi.mm <- sum((y-mu)^2/mu^2)/(n-2)
se.beta <- sqrt(phi.mm * diag(solve(H))) ## Update parameter
                                         ## covariance
summary(glm(y ~ x, family = Gamma)) ## directly from R
cbind(opt$par,se.beta) ## Compare


## profile for phi (because we can calculate c)
phi <- seq(0.002,0.1, by =0.0001)
lp <- numeric(length(phi))
for(i in 1:length(lp)){
    lp[i] <- nll(opt$par,phi[i],X,y)
}
plot(phi,-(lp+max(-lp)),type="l",ylim=c(-3,0))
lines(range(phi),-qchisq(0.95,df=1)/2*c(1,1),lty=2,col=4)

## Hence phi around 0.025
phi.mm
##################################################


##################################################
## Example 6.3
##################################################
dat <- read.table("surgery.dat")
par(mfrow=c(1,1))
plot(dat) ## Bernulli trail

## negative log likelihood function
nll <- function(beta,y,X,n){
    p <- exp(X %*% beta) / (1 + exp(X %*% beta))
    -sum(dbinom(y, size = n, prob = p, log = TRUE))
}

## Observation and design matrix
y <- dat[ ,2]
X <- cbind(1,dat[ ,1]-mean(dat[ ,1]))

opt <- nlminb(c(-1,1),nll,y=y,X=X,n=1)

glm(y~-1+X,family=binomial)
opt$par


## Parameter uncertainty
library(numDeriv)
H <- hessian(nll, opt$par, y=y, X=X, n=1)
se.beta <- sqrt(diag(solve(H)))

summary(fit <- glm(y~-1+X,family=binomial))

## calculating deviances (null deviance p=0.5)
mu0 <- 0.5
Di0 <- y * log(y/mu0) 
Di0[y==0] <- (1-y[y==0]) * log((1-y[y==0])/(1-mu0))

mu1 <- predict(fit, type="response")
Di1 <- y * log(y/mu1) 
Di1[y==0] <- (1-y[y==0]) * log((1-y[y==0])/(1-mu1[y==0]))

2*sum(Di1)
2*sum(Di0)

summary(fit)

lines(dat[ ,1], exp(X %*% opt$par)/(1+exp(X %*% opt$par)))
## Interpreted as the probability of death as a
## function of age

## Or using glm
fit <- glm(y~-1+X,family=binomial)
plot(dat) ## Bernulli trail
lines(dat[ ,1], predict(fit, type="response"))

## More on the fit
?predict.glm
pred <- predict(fit, type="response", 
                inteval="cofidence",se.fit=TRUE)

## Producing confidence intervals
## Wald confidence intervals
lines(dat[ ,1], pred$fit + 2 * pred$se.fit,col=2,lty=2)
lines(dat[ ,1], pred$fit - 2 * pred$se.fit,col=2,lty=2)


## Intervals in linear domain (transformed back to orginal domain)
pred <- predict(fit, type="link", inteval="cofidence",se.fit=TRUE)
lines(dat[ ,1], exp(pred$fit + 2 * pred$se.fit) /
     (1 + exp(pred$fit + 2 * pred$se.fit)),col=3,lty=2)
lines(dat[ ,1], exp(pred$fit - 2 * pred$se.fit) /
     (1 + exp(pred$fit - 2 * pred$se.fit)),col=3,lty=2)




##################################################
# Example 6.6

x<- scan('accident.dat',
        what=list(acc0=0,acc1=0,year0=0,year1=0))
attach(x)

acc <- c(acc0, acc1)
year <- c(year0, year1)
site <- rep(1:8, 2)
treat <- factor(rep(c(0,1), c(8,8)))
cbind(acc, year, site, treat)

## Look at the mean
mean(acc[treat==0])
mean(acc[treat==1])
## Is this fair?



summary(xreg <- glm(acc ~ treat + offset(log(year)),family=poisson))




## Direct calculations of deviance
pred <- predict(xreg,type="response")
d <-  acc * log(acc / pred)
d[is.na(d)] <- 0
d <- d - (acc - pred)
2 * sum(d)

##################################################
## residuals (response)
par(mfrow=c(1,2))
plot(pred, acc - pred)
## Deviance
plot(pred, sign(acc - pred) * sqrt(d))
## More homogene variance
##################################################

##################################################
## Model selection
summary(xreg2 <- glm(acc~ treat + factor(site), offset = log(year), family = poisson))

anova(xreg, xreg2, test = "Chisq")
50.86 - 16.28
1 - pchisq(34.58, df = 7)
## So the bigger model should be used (including the effect of site)

summary(xreg2)
## Lack of fit/ goodness of fit
1 - pchisq(16.275,df=7) ## So indication of misfit


## Overdispersion
summary(xreg3 <- glm(acc~ treat + factor(site), offset = log(year),
             family=quasipoisson))

summary(xreg4 <- glm(acc~ treat,offset=log(year), family=quasipoisson))

anova(xreg4,xreg3,test="Chisq") ## The chisq test
anova(xreg4,xreg3,test="F") ## another test (to be prefered in this case)
## We go with the simpler model 

## Deviance plot
names(xreg4) ## What is in xreg4?
xreg4$deviance
summary(xreg4) ## Which beta1's should we look at 


beta1 <- seq(-2,0.3,length=100) ## beta1
D <- numeric(length(beta1)) ## Initialize deviance

for(i in 1:length(D)){
  D[i] <- glm(acc~ 1+ offset(log(year) + beta1[i]*(treat==1)), ## beta1 as
                                                               ## offset
                 family=quasipoisson)$deviance
}

## plot it 
plot(beta1, (D - xreg4$deviance) / summary(xreg4)$dispersion, type = "l")
lines(range(beta1),qchisq(0.95,df=1)*c(1,1),lty=2,col=2,lwd=2)
## confint use profile deviance
confint(xreg4)
lines(confint(xreg4)[2,1]*c(1,1),c(0,5))
lines(confint(xreg4)[2,2]*c(1,1),c(0,5))

## Wald confidence interval
summary(xreg4)$coefficients[2,1] +
                 summary(xreg4)$coefficients[2,2]*1.96*c(-1,1)


summary(xreg0 <- glm(acc~ 1,offset=log(year), family=quasipoisson))
anova(xreg0,xreg4,test="F")


##################################################
## Example 6.18 (Box-Cox)
rm(list=ls())

# Example 6.18a: to produce Figure 6.8 and the first table 
# in page 180

#Y	SO2 content of air in micrograms per cubic metre
#X1	Average annual temperature in oF
#X2	Number of manufacturing enterprises employing 20 or more workers
#X3	Population size (1970 census); in thousands
#X4	Average annual wind speed in miles per hour
#X5	Average annual precipitation in inches
#X6	Average annual of days with precipitation per year

x<- scan('air.dat',skip=9)
x<- matrix(x,byrow=T,ncol=7)
y<- x[,1]
x<- x[,3]
ord<- order(x)
x<- x[ord]
y<- y[ord]
n<- length(x)

par(mfrow=c(2,2))
plot(x,y,xlab='Industries',
         ylab='Sulphur dioxide',type='n')
  points(x,y,cex=.6)
  title(expression('(a) Pollution data'))

plot(x,y,xlab='Industries',log='x',
         ylab='Sulphur dioxide',type='n')
points(x,y,cex=.6)
title(expression('(b) Transform x-axis'))

lx <- log(x)
lx2 <- lx^2
xreg <- lm(y~lx+lx2)

lines(x,xreg$fit,lty=2)
asum <- summary(xreg)
print(asum)

qqnorm(xreg$res,plot=T)
qqline(xreg$res)
plot(xreg$fitted,xreg$res)
## The normal assumption is violated

library(MASS)
boxcox(xreg)
## Hence pointing to lambda=0

xreg<- lm(log(y) ~ lx + lx2)

qqnorm(xreg$res,plot=T)
qqline(xreg$res)
plot(xreg$fitted,xreg$res)
## Look ok 

## Graphical investigation
par(mfrow=c(2,2))
## lambda = 1
qqnorm(lm(y~lx+lx2)$residuals)
qqline(lm(y~lx+lx2)$residuals)
## lambda = 1/2
qqnorm(lm(sqrt(y)~lx+lx2)$residuals)
qqline(lm(sqrt(y)~lx+lx2)$residuals)
## lambda = 0
qqnorm(lm(log(y)~lx+lx2)$residuals)
qqline(lm(log(y)~lx+lx2)$residuals)



##################################################
## Another example
## Model volume of trees 
?trees

plot(trees)

fit0 <- lm(Volume~Girth,data=trees)
summary(fit0)

## ?????

fit0 <- lm(Volume~Girth*Height+I(Girth^2):Height,
           data=trees)
summary(fit0)
par(mfrow=c(2,2))
plot(fit0)

fit1 <- lm(Volume~-1+I(Girth^2):Height,
           data=trees)
summary(fit1)
anova(fit1,fit0)
plot(fit1)
library(MASS)
boxcox(fit0)

fit2 <- lm(log(Volume)~I(log(Girth))+I(log(Height)),
           data=trees)

par(mfrow=c(2,2))
plot(fit2)
summary(fit2)

##################################################
##
##################################################

