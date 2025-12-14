##################################################
## Continuation – Model fitting (Normal, Gamma, Beta)
## Using log-transformed and untransformed wind speeds
##################################################

library(MASS)

# -------------------------------------------------
# 1. Select transformation and data
# -------------------------------------------------
y.trans <- y.log    # log-transformed
y.raw   <- y        # untransformed (positive)
n <- length(y)

# -------------------------------------------------
# 2. Fit Normal distribution (to log-transformed data)
# -------------------------------------------------
fit.norm <- fitdistr(y.trans, densfun = "normal")
mu.norm  <- fit.norm$estimate["mean"]
sd.norm  <- fit.norm$estimate["sd"]
AIC.norm <- -2 * fit.norm$loglik + 2 * 2

# -------------------------------------------------
# 3. Fit Gamma distribution (to untransformed data)
# -------------------------------------------------
fit.gamma <- fitdistr(y.raw, densfun = "gamma")
shape.hat <- fit.gamma$estimate["shape"]
rate.hat  <- fit.gamma$estimate["rate"]
AIC.gamma <- -2 * fit.gamma$loglik + 2 * 2

# -------------------------------------------------
# 4. Fit Beta distribution (to raw data scaled to (0,1))
# -------------------------------------------------
# Beta requires (0,1) support; rescale to (0,1)
y.min <- min(y.raw)
y.max <- max(y.raw)
y.beta <- (y.raw - y.min) / (y.max - y.min)
eps <- 1e-6
y.beta[y.beta <= 0] <- eps
y.beta[y.beta >= 1] <- 1 - eps

fit.beta <- fitdistr(y.beta, densfun = "beta",
                     start = list(shape1 = 2, shape2 = 5))
alpha.hat <- fit.beta$estimate["shape1"]
beta.hat  <- fit.beta$estimate["shape2"]
AIC.beta  <- -2 * fit.beta$loglik + 2 * 2

# -------------------------------------------------
# 5. Histograms with fitted PDFs
# -------------------------------------------------

# (a) Normal on log-transformed
hist(y.trans, breaks = 30, freq = FALSE, col = "gray90",
     main = "Normal Fit (log-transformed)",
     xlab = "log(Wind Speed)")
curve(dnorm(x, mean = mu.norm, sd = sd.norm),
      add = TRUE, col = "blue", lwd = 2)
legend("topright",
       legend = sprintf("Normal (μ=%.2f, σ=%.2f)", mu.norm, sd.norm),
       col = "blue", lwd = 2, bty = "n")

# (b) Gamma on untransformed
hist(y.raw, breaks = 30, freq = FALSE, col = "gray90",
     main = "Gamma Fit (ws30)",
     xlab = "Wind Speed (m/s)")
curve(dgamma(x, shape = shape.hat, rate = rate.hat),
      add = TRUE, col = "red", lwd = 2)
legend("topright",
       legend = expression(Gamma~(alpha~","~theta)),
       col = "red", lwd = 2, bty = "n")

# (c) Beta on raw (scaled internally)
hist(y.raw, breaks = 30, freq = FALSE, col = "gray90",
     main = "Beta Fit (ws30)",
     xlab = "Wind Speed (m/s)")
curve(dbeta((x - y.min)/(y.max - y.min),
            shape1 = alpha.hat, shape2 = beta.hat) / (y.max - y.min),
      add = TRUE, col = "darkgreen", lwd = 2)
legend("topright",
       legend = expression(Beta~(alpha~","~beta)),
       col = "darkgreen", lwd = 2, bty = "n")

# -------------------------------------------------
# 6. QQ-Plots
# -------------------------------------------------
par(mfrow = c(1,3))

# (a) Normal QQ
q.norm <- qnorm(ppoints(n), mean = mu.norm, sd = sd.norm)
plot(q.norm, sort(y.trans),
     main = "Normal Fit (log-transformed)",
     xlab = "Theoretical quantiles (Normal)",
     ylab = "Empirical quantiles (log(ws30))")
abline(0, 1, col = "blue", lwd = 2)

# (b) Gamma QQ
q.gamma <- qgamma(ppoints(n), shape = shape.hat, rate = rate.hat)
plot(q.gamma, sort(y.raw),
     main = "Gamma Fit (ws30)",
     xlab = "Theoretical quantiles (Gamma)",
     ylab = "Empirical quantiles (ws30)")
abline(0, 1, col = "red", lwd = 2)

# (c) Beta QQ
q.beta <- qbeta(ppoints(n), shape1 = alpha.hat, shape2 = beta.hat)
plot(y.min + q.beta * (y.max - y.min), sort(y.raw),
     main = "Beta Fit (ws30)",
     xlab = "Theoretical quantiles (Beta)",
     ylab = "Empirical quantiles (ws30)")
abline(0, 1, col = "darkgreen", lwd = 2)

par(mfrow = c(1,1))
