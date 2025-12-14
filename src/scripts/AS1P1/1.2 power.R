##################################################
## Assignment 1 – Normal (Eq 1 transformed), Beta, Gamma fits
## Based on Lecture 4 Example 4.17 (MLE style)
##################################################

rm(list = ls())

# -------------------------------------------------
# 1. Load and normalize data
# -------------------------------------------------
df_tuno <- read.table("tuno.txt", header = TRUE)
y <- df_tuno$pow.obs / 5000         # normalize to (0,1)
eps <- 1e-6
y[y <= 0] <- eps
y[y >= 1] <- 1 - eps
n <- length(y)

# -------------------------------------------------
# 2. Transformation Eq (1)
# -------------------------------------------------
# y(λ) = (1/λ) * log( y^λ / (1 - y^λ) )
trans1 <- function(lambda, y){
  (1 / lambda) * log((y^lambda) / (1 - y^lambda))
}

# Use λ found earlier (≈ 0.2 – 0.3 typically)
lambda.eq1 <- 0.326
y.trans <- trans1(lambda.eq1, y)

##################################################
## NORMAL DISTRIBUTION FIT (on transformed data)
##################################################

mu.norm <- mean(y.trans)
sd.norm <- sd(y.trans)
logLik.norm <- sum(dnorm(y.trans, mu.norm, sd.norm, log = TRUE))
AIC.norm <- -2 * logLik.norm + 2 * 2

hist(y.trans, breaks = 30, freq = FALSE, col = "gray90",
     main = sprintf("Normal Fit on Eq (1) Transformed Power (λ=%.2f)", lambda.eq1),
     xlab = "Transformed Power")
curve(dnorm(x, mu.norm, sd.norm), add = TRUE, col = "blue", lwd = 2, n = 2000)
legend("topleft",
       legend = sprintf("Normal (μ=%.1f, σ=%.1f)", mu.norm, sd.norm),
       col = "blue", lwd = 2, bty = "n")

##################################################
## BETA DISTRIBUTION FIT (MLE on original y)
##################################################

nll.beta <- function(pars, x){
  a <- pars[1]; b <- pars[2]
  if (a <= 0 || b <= 0) return(Inf)
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}
opt.beta <- nlminb(start = c(2, 5), objective = nll.beta,
                   x = y, lower = c(1e-6, 1e-6))

alpha.hat <- opt.beta$par[1]
beta.hat  <- opt.beta$par[2]
AIC.beta <- 2 * opt.beta$objective + 2 * 2

hist(y, breaks = 30, freq = FALSE, col = "gray90",
     main = "Beta Fit to Normalized Power",
     xlab = "Normalized Power (0–1)")
curve(dbeta(x, shape1 = alpha.hat, shape2 = beta.hat),
      add = TRUE, col = "darkgreen", lwd = 2, from = 1e-6, to = 1, n = 2000)
legend("topright",
       legend = sprintf("Beta fit (α=%.2f, β=%.2f)", alpha.hat, beta.hat),
       col = "darkgreen", lwd = 2, bty = "n")

##################################################
## GAMMA DISTRIBUTION FIT (MLE on original y)
##################################################

nll.gamma <- function(pars, x){
  s <- pars[1]; r <- pars[2]
  if (s <= 0 || r <= 0) return(Inf)
  -sum(dgamma(x, shape = s, rate = r, log = TRUE))
}
opt.gamma <- nlminb(start = c(1, 1), objective = nll.gamma,
                    x = y, lower = c(1e-6, 1e-6))

shape.hat <- opt.gamma$par[1]
rate.hat  <- opt.gamma$par[2]
AIC.gamma <- 2 * opt.gamma$objective + 2 * 2

hist(y, breaks = 30, freq = FALSE, col = "gray90",
     main = "Gamma Fit to Normalized Power",
     xlab = "Normalized Power (0–1)")
curve(dgamma(x, shape = shape.hat, rate = rate.hat),
      add = TRUE, col = "red", lwd = 2, from = 1e-6, to = 1, n = 2000)
legend("topright",
       legend = sprintf("Gamma fit (shape=%.2f, rate=%.2f)", shape.hat, rate.hat),
       col = "red", lwd = 2, bty = "n")

##################################################
## QQ-PLOTS
##################################################
##################################################
## QQ-PLOTS
##################################################

# --- All three QQ-plots: Normal (transformed), Beta, Gamma ---
par(mfrow = c(1,3))  # 3 plots side by side

# (a) Normal QQ (Eq. 1 transformed)
q.empir.trans <- sort(y.trans)
q.norm <- qnorm(ppoints(n), mean = mu.norm, sd = sd.norm)
plot(q.norm, q.empir.trans,
     main = "Normal Fit (Eq. 1 Transformed)",
     xlab = "Theoretical quantiles (Normal)",
     ylab = "Empirical quantiles (Transformed)")
abline(0, 1, col = "blue", lwd = 2)

# (b) Beta QQ (untransformed)
q.empir <- sort(y)
q.beta <- qbeta(ppoints(n), shape1 = alpha.hat, shape2 = beta.hat)
plot(q.beta, q.empir,
     main = "Beta Fit (Normalized y)",
     xlab = "Theoretical quantiles (Beta)",
     ylab = "Empirical quantiles")
abline(0, 1, col = "darkgreen", lwd = 2)

# (c) Gamma QQ (untransformed)
q.gamma <- qgamma(ppoints(n), shape = shape.hat, rate = rate.hat)
plot(q.gamma, q.empir,
     main = "Gamma Fit (Normalized y)",
     xlab = "Theoretical quantiles (Gamma)",
     ylab = "Empirical quantiles")
abline(0, 1, col = "red", lwd = 2)

par(mfrow = c(1,1))  # reset layout


##################################################
## AIC COMPARISON
##################################################
cat("\n--------------------------------------------------\n")
cat("Model comparison (lower AIC = better fit):\n")
cat(sprintf("  Normal (Eq 1 transformed) AIC = %.2f\n", AIC.norm))
cat(sprintf("  Beta   (normalized)       AIC = %.2f\n", AIC.beta))
cat(sprintf("  Gamma  (normalized)       AIC = %.2f\n", AIC.gamma))
