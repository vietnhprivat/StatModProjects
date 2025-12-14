## Assignment 1, Part 1 – Transformations for Normalized Power
## Based on Lecture 4 (Example 4.16) and assignment equations (1) & (2)

# ---------------------------------------------------------------
# 1. Load and normalize data
# ---------------------------------------------------------------
df_tuno <- read.table("tuno.txt", header = TRUE)
df_tuno$pow.obs <- df_tuno$pow.obs / 5000        # normalize to (0,1)
y <- df_tuno$pow.obs
df_tuno
summary(df_tuno)
# ---------------------------------------------------------------
# 2. Define transformation functions
# ---------------------------------------------------------------

# --- Standard Box–Cox (Lecture 4)
bc.trans <- function(lambda, y) {
  if (lambda == 0) log(y) else (y^lambda - 1) / lambda
}

# --- Assignment transformation (1)
# y(λ) = (1/λ) * log( y^λ / (1 - y^λ) )
trans1 <- function(lambda, y) {
  (1 / lambda) * log((y^lambda) / (1 - y^lambda))
}

# --- Assignment transformation (2)
# y(λ) = 2 * log( y^λ / (1 - y)^(1 - λ) )
trans2 <- function(lambda, y) {
  2 * log((y^lambda) / ((1 - y)^(1 - lambda)))
}

# ---------------------------------------------------------------
# 3. Profile log-likelihood function (from Lecture 4)
# ---------------------------------------------------------------
lp.lambda <- function(lambda, y, trans.fun) {
  y.l <- trans.fun(lambda, y)
  if (any(!is.finite(y.l))) return(NA)     # skip invalid λ
  n <- length(y)
  sigmasq <- mean((y.l - mean(y.l))^2)
  -n/2 * log(sigmasq) + (lambda - 1) * sum(log(y))
}

# ---------------------------------------------------------------
# 4. Compute profile likelihoods
# ---------------------------------------------------------------
lambda.seq.bc <- seq(-0.5, 0.5, by = 0.01)
lambda.seq.t1 <- seq(0.001, 2, by = 0.01)
lambda.seq.t2 <- seq(0.001, 0.999, by = 0.01)

lp.bc <- sapply(lambda.seq.bc, lp.lambda, y = y, trans.fun = bc.trans)
lp.t1 <- sapply(lambda.seq.t1, lp.lambda, y = y, trans.fun = trans1)
lp.t2 <- sapply(lambda.seq.t2, lp.lambda, y = y, trans.fun = trans2)

# ---------------------------------------------------------------
# 5. Plot profile likelihoods
# ---------------------------------------------------------------
par(mfrow = c(1,3))
plot(lambda.seq.bc, lp.bc - max(lp.bc, na.rm=TRUE), type = "l",
     main = "Profile Likelihood (Box–Cox)",
     xlab = expression(lambda), ylab = "Relative log-likelihood")
abline(h = -qchisq(0.95, df = 1)/2, col = 2, lty = 2)

plot(lambda.seq.t1, lp.t1 - max(lp.t1, na.rm=TRUE), type = "l",
     main = "Profile Likelihood (Eq. 1)",
     xlab = expression(lambda), ylab = "Relative log-likelihood")
abline(h = -qchisq(0.95, df = 1)/2, col = 2, lty = 2)

plot(lambda.seq.t2, lp.t2 - max(lp.t2, na.rm=TRUE), type = "l",
     main = "Profile Likelihood (Eq. 2)",
     xlab = expression(lambda), ylab = "Relative log-likelihood")
abline(h = -qchisq(0.95, df = 1)/2, col = 2, lty = 2)
par(mfrow = c(1,1))

# ---------------------------------------------------------------
# 6. Find optimal λ for each transformation (Lecture 4 syntax)
# ---------------------------------------------------------------
opt.bc <- optimize(lp.lambda, c(-2, 2), y = y, trans.fun = bc.trans, maximum = TRUE)$maximum
opt.t1 <- optimize(lp.lambda, c(0.001, 2), y = y, trans.fun = trans1, maximum = TRUE)$maximum
opt.t2 <- optimize(lp.lambda, c(0.001, 0.999), y = y, trans.fun = trans2, maximum = TRUE)$maximum

opt.bc; opt.t1; opt.t2

# ---------------------------------------------------------------
# 7. Apply transformations using best λ
# ---------------------------------------------------------------
y.bc <- bc.trans(opt.bc, y)
y.t1 <- trans1(opt.t1, y)
y.t2 <- trans2(opt.t2, y)

# ---------------------------------------------------------------
# 8. Plot histograms before and after transformation
# ---------------------------------------------------------------
par(mfrow = c(2,2))
hist(y, breaks=30, col="lightblue", main="Original normalized power",
     xlab="Normalized power (0–1)")
hist(y.bc, breaks=30, col="lightgreen",
     main=paste0("Box–Cox (λ=", round(opt.bc,3),")"),
     xlab="Transformed power")
hist(y.t1, breaks=30, col="orange",
     main=paste0("Eq. (1) (λ=", round(opt.t1,3),")"),
     xlab="Transformed power")
hist(y.t2, breaks=30, col="pink",
     main=paste0("Eq. (2) (λ=", round(opt.t2,3),")"),
     xlab="Transformed power")
par(mfrow = c(1,1))

# ---------------------------------------------------------------
# 9. QQ-plots before and after transformation
# ---------------------------------------------------------------
par(mfrow = c(2,2))
qqnorm(y, main="Original normalized power"); qqline(y)
qqnorm(y.bc, main=paste0("Box–Cox (λ=", round(opt.bc,3),")")); qqline(y.bc)
qqnorm(y.t1, main=paste0("Eq. (1) (λ=", round(opt.t1,3),")")); qqline(y.t1)
qqnorm(y.t2, main=paste0("Eq. (2) (λ=", round(opt.t2,3),")")); qqline(y.t2)
par(mfrow = c(1,1))

# ===============================================================
# 10. Model fitting for normalized power (continuation)
# ===============================================================

# --- Choose the best transformation (Eq. 1)
y.trans <- y.t1   # Eq. (1)-transformed power
y.raw   <- y      # Untransformed normalized power (0–1)

# ===============================================================
# 10.1 Fit distributions using Lecture 4 methods (MASS::fitdistr)
# ===============================================================
library(MASS)

# --- (a) Normal fit (for transformed data)
fit.normal <- fitdistr(y.trans, densfun = "normal")
normal.mean <- fit.normal$estimate["mean"]
normal.sd   <- fit.normal$estimate["sd"]
AIC.normal  <- -2 * fit.normal$loglik + 2 * 2

# --- (b) Gamma fit (for untransformed normalized data)
fit.gamma <- fitdistr(y.raw, densfun = "gamma")
gamma.shape <- fit.gamma$estimate["shape"]
gamma.rate  <- fit.gamma$estimate["rate"]
AIC.gamma   <- -2 * fit.gamma$loglik + 2 * 2

# --- (c) Beta fit (for untransformed normalized data)
fit.beta <- fitdistr(y.raw, densfun = "beta",
                     start = list(shape1 = 2, shape2 = 5))
beta.a  <- fit.beta$estimate["shape1"]
beta.b  <- fit.beta$estimate["shape2"]
AIC.beta <- -2 * fit.beta$loglik + 2 * 2

# ===============================================================
# 11. Visual evaluation (histograms + fitted PDFs)
# ===============================================================

# --- (a) Normal fit on transformed data
hist(y.trans, breaks = 30, freq = FALSE, col = "gray90",
     main = "Eq. (1)-Transformed Power (Normal Fit)",
     xlab = "Transformed Power")
curve(dnorm(x, mean = normal.mean, sd = normal.sd),
      add = TRUE, col = "blue", lwd = 2)
legend("topright", legend = "Normal fit", col = "blue", lwd = 2, bty = "n")

# ===============================================================
# 12. Quantitative comparison (log-likelihoods and AIC)
# ===============================================================
results <- data.frame(
  Model  = c("Normal (transformed)", "Gamma (raw)", "Beta (raw)"),
  LogLik = c(fit.normal$loglik, fit.gamma$loglik, fit.beta$loglik),
  AIC    = c(AIC.normal, AIC.gamma, AIC.beta)
)
print(results)


# ===============================================================
# Assignment 1, Part 2 – Transformations for ws30
# ===============================================================

# ---------------------------------------------------------------
# 1. Load and normalize data
# ---------------------------------------------------------------
df_tuno <- read.table("tuno.txt", header = TRUE)

# Normalize ws30 to (0,1)
df_tuno$ws30 <- df_tuno$ws30 / max(df_tuno$ws30, na.rm = TRUE)
y <- df_tuno$ws30

# Avoid 0 or 1 values
eps <- 1e-6
y[y <= 0] <- eps
y[y >= 1] <- 1 - eps

# ---------------------------------------------------------------
# 2. Compute profile likelihoods for ws30
# ---------------------------------------------------------------
lambda.seq.bc <- seq(-0.5, 0.5, by = 0.01)
lambda.seq.t1 <- seq(0.001, 2, by = 0.01)
lambda.seq.t2 <- seq(0.001, 0.999, by = 0.01)

lp.bc <- sapply(lambda.seq.bc, lp.lambda, y = y, trans.fun = bc.trans)
lp.t1 <- sapply(lambda.seq.t1, lp.lambda, y = y, trans.fun = trans1)
lp.t2 <- sapply(lambda.seq.t2, lp.lambda, y = y, trans.fun = trans2)

# ---------------------------------------------------------------
# 3. Plot profile likelihoods
# ---------------------------------------------------------------
par(mfrow = c(1,3))
plot(lambda.seq.bc, lp.bc - max(lp.bc, na.rm=TRUE), type="l",
     main = "Profile Likelihood (Box–Cox, ws30)",
     xlab = expression(lambda), ylab = "Relative log-likelihood")
abline(h = -qchisq(0.95, df = 1)/2, col=2, lty=2)

plot(lambda.seq.t1, lp.t1 - max(lp.t1, na.rm=TRUE), type="l",
     main = "Profile Likelihood (Eq. 1, ws30)",
     xlab = expression(lambda), ylab = "Relative log-likelihood")
abline(h = -qchisq(0.95, df = 1)/2, col=2, lty=2)

plot(lambda.seq.t2, lp.t2 - max(lp.t2, na.rm=TRUE), type="l",
     main = "Profile Likelihood (Eq. 2, ws30)",
     xlab = expression(lambda), ylab = "Relative log-likelihood")
abline(h = -qchisq(0.95, df = 1)/2, col=2, lty=2)
par(mfrow = c(1,1))

# ---------------------------------------------------------------
# 4. Find optimal λ
# ---------------------------------------------------------------
opt.bc <- optimize(lp.lambda, c(-2, 2), y = y, trans.fun = bc.trans, maximum = TRUE)$maximum
opt.t1 <- optimize(lp.lambda, c(0.001, 2), y = y, trans.fun = trans1, maximum = TRUE)$maximum
opt.t2 <- optimize(lp.lambda, c(0.001, 0.999), y = y, trans.fun = trans2, maximum = TRUE)$maximum

opt.bc; opt.t1; opt.t2

# ---------------------------------------------------------------
# 5. Apply transformations using best λ
# ---------------------------------------------------------------
y.bc <- bc.trans(opt.bc, y)
y.t1 <- trans1(opt.t1, y)
y.t2 <- trans2(opt.t2, y)

# ---------------------------------------------------------------
# 6. Histograms
# ---------------------------------------------------------------
par(mfrow = c(2,2))
hist(y, breaks=30, col="lightblue", main="Original ws30 (normalized)",
     xlab="Normalized ws30 (0–1)")
hist(y.bc, breaks=30, col="lightgreen",
     main=paste0("Box–Cox (λ=", round(opt.bc,3),")"),
     xlab="Transformed ws30")
hist(y.t1, breaks=30, col="orange",
     main=paste0("Eq. (1) (λ=", round(opt.t1,3),")"),
     xlab="Transformed ws30")
hist(y.t2, breaks=30, col="pink",
     main=paste0("Eq. (2) (λ=", round(opt.t2,3),")"),
     xlab="Transformed ws30")
par(mfrow = c(1,1))

# ---------------------------------------------------------------
# 7. QQ-plots
# ---------------------------------------------------------------
par(mfrow = c(2,2))
qqnorm(y, main="Original ws30 (normalized)"); qqline(y)
qqnorm(y.bc, main=paste0("Box–Cox (λ=", round(opt.bc,3),")")); qqline(y.bc)
qqnorm(y.t1, main=paste0("Eq. (1) (λ=", round(opt.t1,3),")")); qqline(y.t1)
qqnorm(y.t2, main=paste0("Eq. (2) (λ=", round(opt.t2,3),")")); qqline(y.t2)
par(mfrow = c(1,1))

# ---------------------------------------------------------------
# 8. Fit distributions and compare (choose best transformation)
# ---------------------------------------------------------------
y.trans <- y.t1   # Example: Eq. (1) works well again
y.raw   <- y

# --- Normal fit
fit.normal <- fitdistr(y.trans, densfun = "normal")
normal.mean <- fit.normal$estimate["mean"]
normal.sd   <- fit.normal$estimate["sd"]
AIC.normal  <- -2 * fit.normal$loglik + 2 * 2

# --- Gamma fit
fit.gamma <- fitdistr(y.raw, densfun = "gamma")
gamma.shape <- fit.gamma$estimate["shape"]
gamma.rate  <- fit.gamma$estimate["rate"]
AIC.gamma   <- -2 * fit.gamma$loglik + 2 * 2

# --- Beta fit
fit.beta <- fitdistr(y.raw, densfun = "beta", start = list(shape1 = 2, shape2 = 5))
beta.a  <- fit.beta$estimate["shape1"]
beta.b  <- fit.beta$estimate["shape2"]
AIC.beta <- -2 * fit.beta$loglik + 2 * 2

# ---------------------------------------------------------------
# 9. Plot Normal fit for transformed data
# ---------------------------------------------------------------
hist(y.trans, breaks=30, freq=FALSE, col="gray90",
     main="Eq. (1)-Transformed ws30 (Normal Fit)",
     xlab="Transformed ws30")
curve(dnorm(x, mean=normal.mean, sd=normal.sd), add=TRUE, col="blue", lwd=2)
legend("topright", legend="Normal fit", col="blue", lwd=2, bty="n")

# ---------------------------------------------------------------
# 10. Compare model AIC values
# ---------------------------------------------------------------
results_ws30 <- data.frame(
  Model  = c("Normal (transformed)", "Gamma (raw)", "Beta (raw)"),
  LogLik = c(fit.normal$loglik, fit.gamma$loglik, fit.beta$loglik),
  AIC    = c(AIC.normal, AIC.gamma, AIC.beta)
)
print(results_ws30)


