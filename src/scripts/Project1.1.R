################################################################################
# INTEGRATED WIND DATA ANALYSIS: POWER, SPEED, AND DIRECTION
# Merged from: 1.1 power.R, 1.2 power.R, 1.1 ws.R, 1.2 ws.R, 1.2.2 ws.R, 1.1 wd.R
################################################################################

rm(list = ls())
library(MASS)

# ==============================================================================
# 0. GLOBAL SETUP & HELPER FUNCTIONS
# ==============================================================================

# 0.1 Load Data (Once for the whole script)
# ------------------------------------------------------------------------------
df_tuno <- read.table("data/tuno.txt", header = TRUE)
cat("Data Loaded. Dimensions:", dim(df_tuno), "\n")
head(df_tuno)

# 0.2 Define Transformation Functions
# ------------------------------------------------------------------------------
# Standard Box-Cox
bc.trans <- function(lambda, y) {
  if (abs(lambda) < 1e-6) log(y) else (y^lambda - 1) / lambda
}

# Transformation Eq (1): y(λ) = (1/λ) * log( y^λ / (1 - y^λ) )
trans1 <- function(lambda, y) {
  (1 / lambda) * log((y^lambda) / (1 - y^lambda))
}

# Transformation Eq (2): y(λ) = 2 * log( y^λ / (1 - y)^(1 - λ) )
trans2 <- function(lambda, y) {
  2 * log((y^lambda) / ((1 - y)^(1 - lambda)))
}

# 0.3 Generic Profile Log-Likelihood Function
# ------------------------------------------------------------------------------
lp.lambda <- function(lambda, y, trans.fun) {
  y.l <- trans.fun(lambda, y)
  # Handle potential numerical issues with extreme lambdas
  if (any(!is.finite(y.l))) return(NA)
  
  n <- length(y)
  # MLE for sigma^2
  sigma2 <- mean((y.l - mean(y.l))^2)
  
  # Profile log-likelihood
  # Jacobian adjustment depends on specific transformation, 
  # but for Box-Cox and derived assignments, the standard form is often used:
  # -n/2 * log(sigma2) + (lambda - 1) * sum(log(y)) 
  # Note: The specific Jacobian term (lambda-1)*sum(log(y)) is specific to Box-Cox.
  # For exact correctness with Eq1/Eq2, one should derive the specific Jacobian.
  # However, following the provided source code structure, we apply the logic used in 1.1 power.R:
  -n/2 * log(sigma2) + (lambda - 1) * sum(log(y))
}


# ==============================================================================
# PART 1: WIND POWER ANALYSIS (pow.obs)
# ==============================================================================
cat("\n================ PART 1: WIND POWER ================\n")

# 1.1 Data Preparation
# ------------------------------------------------------------------------------
# Normalize to (0,1)
y.pow <- df_tuno$pow.obs / 5000 
# Handle boundary values for Log/Beta transformations
eps <- 1e-6
y.pow[y.pow <= 0] <- eps
y.pow[y.pow >= 1] <- 1 - eps
n.pow <- length(y.pow)

# 1.2 Transformations & Profile Likelihood
# ------------------------------------------------------------------------------
# We test three transformations: Box-Cox, Eq1, Eq2
lambda.seq.bc <- seq(-0.5, 0.5, by = 0.01)
lambda.seq.t1 <- seq(0.001, 2, by = 0.01)
lambda.seq.t2 <- seq(0.001, 0.999, by = 0.01)

lp.bc <- sapply(lambda.seq.bc, lp.lambda, y = y.pow, trans.fun = bc.trans)
lp.t1 <- sapply(lambda.seq.t1, lp.lambda, y = y.pow, trans.fun = trans1)
lp.t2 <- sapply(lambda.seq.t2, lp.lambda, y = y.pow, trans.fun = trans2)

# Find optimal Lambdas
opt.bc <- optimize(lp.lambda, c(-2, 2), y = y.pow, trans.fun = bc.trans, maximum = TRUE)$maximum
opt.t1 <- optimize(lp.lambda, c(0.001, 2), y = y.pow, trans.fun = trans1, maximum = TRUE)$maximum
opt.t2 <- optimize(lp.lambda, c(0.001, 0.999), y = y.pow, trans.fun = trans2, maximum = TRUE)$maximum

cat(sprintf("Optimal Lambdas -> BoxCox: %.3f, Eq1: %.3f, Eq2: %.3f\n", opt.bc, opt.t1, opt.t2))

# 1.3 Visualizing Transformations
# ------------------------------------------------------------------------------
par(mfrow = c(1, 3))
hist(bc.trans(opt.bc, y.pow), breaks=30, col="lightgreen", main=sprintf("Box-Cox (L=%.2f)", opt.bc), xlab="Transformed")
hist(trans1(opt.t1, y.pow), breaks=30, col="orange", main=sprintf("Eq(1) (L=%.2f)", opt.t1), xlab="Transformed")
hist(trans2(opt.t2, y.pow), breaks=30, col="pink", main=sprintf("Eq(2) (L=%.2f)", opt.t2), xlab="Transformed")
par(mfrow = c(1,1))

# 1.4 Model Fitting & AIC Comparison
# ------------------------------------------------------------------------------
# Strategy: Compare Normal (on Eq1 transformed data) vs Gamma (raw) vs Beta (raw)

# A. Normal Fit on Transformed Data (Using Eq 1 as it typically performs best for Power)
y.trans.pow <- trans1(opt.t1, y.pow)
fit.norm.pow <- fitdistr(y.trans.pow, densfun = "normal")
AIC.norm.pow <- -2 * fit.norm.pow$loglik + 2 * 2

# B. Gamma Fit on Raw Data
fit.gamma.pow <- fitdistr(y.pow, densfun = "gamma")
AIC.gamma.pow <- -2 * fit.gamma.pow$loglik + 2 * 2

# C. Beta Fit on Raw Data
fit.beta.pow  <- fitdistr(y.pow, densfun = "beta", start = list(shape1 = 2, shape2 = 5))
AIC.beta.pow  <- -2 * fit.beta.pow$loglik + 2 * 2

# 1.5 Results Output
# ------------------------------------------------------------------------------
cat("\n--- Power Model Comparison (AIC) ---\n")
cat(sprintf("Normal (Eq1 Transformed): %.2f\n", AIC.norm.pow))
cat(sprintf("Gamma  (Raw Data):        %.2f\n", AIC.gamma.pow))
cat(sprintf("Beta   (Raw Data):        %.2f\n", AIC.beta.pow))

# 1.6 QQ Plots Comparison
# ------------------------------------------------------------------------------
par(mfrow = c(1,3))
# Normal (Transformed)
qqnorm(y.trans.pow, main="QQ: Normal (Eq1 Transformed)"); qqline(y.trans.pow, col="blue")
# Gamma (Raw)
q.gamma <- qgamma(ppoints(n.pow), shape = fit.gamma.pow$estimate["shape"], rate = fit.gamma.pow$estimate["rate"])
plot(q.gamma, sort(y.pow), main="QQ: Gamma (Raw)", xlab="Theoretical", ylab="Empirical"); abline(0,1, col="red")
# Beta (Raw)
q.beta <- qbeta(ppoints(n.pow), shape1 = fit.beta.pow$estimate["shape1"], shape2 = fit.beta.pow$estimate["shape2"])
plot(q.beta, sort(y.pow), main="QQ: Beta (Raw)", xlab="Theoretical", ylab="Empirical"); abline(0,1, col="green")
par(mfrow = c(1,1))


# ==============================================================================
# PART 2: WIND SPEED ANALYSIS (ws30)
# ==============================================================================
cat("\n================ PART 2: WIND SPEED ================\n")

# 2.1 Data Preparation
# ------------------------------------------------------------------------------
y.ws <- df_tuno$ws30
y.ws <- y.ws[y.ws > 0] # strictly positive
n.ws <- length(y.ws)

# 2.2 Transformations (Box-Cox)
# ------------------------------------------------------------------------------
# Find optimal lambda for Box-Cox
lambda.grid <- seq(-1, 1, by = 0.01)
lp.ws <- sapply(lambda.grid, lp.lambda, y = y.ws, trans.fun = bc.trans)
lambda.hat.ws <- lambda.grid[which.max(lp.ws)]

cat(sprintf("Optimal Box-Cox Lambda for Wind Speed: %.2f\n", lambda.hat.ws))

# 2.3 Comparison A: Standard Distributions
# ------------------------------------------------------------------------------
# Compares: Normal (Log-transformed), Gamma (Raw), Beta (Rescaled Raw)

y.ws.log <- log(y.ws)

# Rescale for Beta (0,1)
y.ws.min <- min(y.ws); y.ws.max <- max(y.ws)
y.ws.beta <- (y.ws - y.ws.min) / (y.ws.max - y.ws.min)
y.ws.beta[y.ws.beta <= 0] <- eps; y.ws.beta[y.ws.beta >= 1] <- 1 - eps

# Fits
fit.norm.ws  <- fitdistr(y.ws.log, "normal")
fit.gamma.ws <- fitdistr(y.ws, "gamma")
fit.beta.ws  <- fitdistr(y.ws.beta, "beta", start = list(shape1 = 2, shape2 = 5))

# AICs
AIC.norm.ws  <- -2 * fit.norm.ws$loglik + 2 * 2
AIC.gamma.ws <- -2 * fit.gamma.ws$loglik + 2 * 2
AIC.beta.ws  <- -2 * fit.beta.ws$loglik + 2 * 2

cat("\n--- Speed Comparison A: Standard Fits ---\n")
cat(sprintf("Normal (Log-trans): %.2f\n", AIC.norm.ws))
cat(sprintf("Gamma  (Raw):       %.2f\n", AIC.gamma.ws))
cat(sprintf("Beta   (Rescaled):  %.2f\n", AIC.beta.ws))

# 2.4 Comparison B: Heavy Tail Models (Log Domain)
# ------------------------------------------------------------------------------
# Compares Normal vs Cauchy vs Student-t on the Log-Transformed data

# Cauchy Fit
fit.cauchy.ws <- fitdistr(y.ws.log, "cauchy")
AIC.cauchy.ws <- -2 * fit.cauchy.ws$loglik + 2 * 2

# Student-t Fit (requires custom optimization as fitdistr doesn't fully support 't')
nll.t <- function(pars, x) {
  mu <- pars[1]; sigma <- pars[2]; nu <- pars[3]
  if (sigma <= 0 || nu <= 0) return(Inf)
  -sum(dt((x - mu)/sigma, df = nu, log = TRUE) - log(sigma))
}

opt.t <- nlminb(start = c(mean(y.ws.log), sd(y.ws.log), 5),
                objective = nll.t,
                x = y.ws.log,
                lower = c(-Inf, 1e-6, 0.5))

AIC.t.ws <- 2 * opt.t$objective + 2 * 3

cat("\n--- Speed Comparison B: Log-Domain Heavy Tails ---\n")
cat(sprintf("Normal:    %.2f\n", AIC.norm.ws))
cat(sprintf("Cauchy:    %.2f\n", AIC.cauchy.ws))
cat(sprintf("Student-t: %.2f (nu = %.2f)\n", AIC.t.ws, opt.t$par[3]))

# 2.5 Visualizing the Best Heavy Tail Fit
# ------------------------------------------------------------------------------
par(mfrow = c(1,3))
# Normal QQ
q.norm <- qnorm(ppoints(n.ws), mean=fit.norm.ws$estimate[1], sd=fit.norm.ws$estimate[2])
plot(q.norm, sort(y.ws.log), main="Normal Fit (Log)", ylab="Empirical"); abline(0,1, col="blue")

# Cauchy QQ
q.cauchy <- qcauchy(ppoints(n.ws), location=fit.cauchy.ws$estimate[1], scale=fit.cauchy.ws$estimate[2])
plot(q.cauchy, sort(y.ws.log), main="Cauchy Fit (Log)", ylab="Empirical"); abline(0,1, col="red")

# Student-t QQ
mu.t <- opt.t$par[1]; sigma.t <- opt.t$par[2]; nu.t <- opt.t$par[3]
q.t <- mu.t + sigma.t * qt(ppoints(n.ws), df = nu.t)
plot(q.t, sort(y.ws.log), main=sprintf("Student-t (nu=%.1f)", nu.t), ylab="Empirical"); abline(0,1, col="darkgreen")
par(mfrow = c(1,1))


# ==============================================================================
# PART 3: WIND DIRECTION ANALYSIS (wd30)
# ==============================================================================
cat("\n================ PART 3: WIND DIRECTION ================\n")

# 3.1 Data Preparation
# ------------------------------------------------------------------------------
theta <- df_tuno$wd30 # in radians

# 3.2 Von Mises Functions
# ------------------------------------------------------------------------------
# Negative log-likelihood (optimizing kappa for fixed mu)
nll_kappa <- function(kappa, mu, x) {
  if (kappa <= 0) return(Inf)
  -sum( kappa*cos(x - mu) - log(2*pi) - log(besselI(kappa, 0)) )
}

# Profile log-likelihood for mu (optimizing kappa inside)
lp_mu <- function(mu, x) {
  opt <- optimize(f = nll_kappa, mu = mu, x = x, interval = c(0.0001, 50))
  return(opt$objective) # Note: this returns negative log-likelihood (minimized)
}

# 3.3 Fitting Parameters
# ------------------------------------------------------------------------------
# Grid search for Mu
mu_grid <- seq(-pi, pi, length.out = 400)
# We use 'lp_mu' which returns negative log-lik, so we minimize it
nll_vals <- sapply(mu_grid, lp_mu, x = theta) 
mu_hat <- mu_grid[which.min(nll_vals)]

# Optimize Kappa at best Mu
opt_kappa <- optimize(f = nll_kappa, mu = mu_hat, x = theta, interval = c(0.0001, 50))
kappa_hat <- opt_kappa$minimum

cat(sprintf("Von Mises Estimates -> Mu: %.3f, Kappa: %.3f\n", mu_hat, kappa_hat))

# 3.4 Visualization
# ------------------------------------------------------------------------------
vonmises_pdf <- function(x, mu, kappa) {
  exp(kappa*cos(x - mu)) / (2*pi*besselI(kappa, 0))
}

hist(theta, breaks = 40, freq = FALSE, col = "lightgray",
     main = "Wind Direction: Von Mises Fit",
     xlab = "Direction (radians)")
curve(vonmises_pdf(x, mu_hat, kappa_hat), add = TRUE, lwd = 3, col = "red")