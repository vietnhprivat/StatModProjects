##################################################
## Assignment 1 – Comparing Normal, Cauchy & Student-t fits
## Using log-transformed wind speeds (ws30)
## Based on Lecture 4 Example 4.17 (MLE style)
##################################################

rm(list = ls())
library(MASS)

# -------------------------------------------------
# 1. Load and transform data
# -------------------------------------------------
df_tuno <- read.table("tuno.txt", header = TRUE)
y <- df_tuno$ws30
y <- y[y > 0]                 # remove zeros for log-transform
y.log <- log(y)
n <- length(y.log)

# -------------------------------------------------
# 2. Fit Normal distribution (to log(ws30))
# -------------------------------------------------
fit.norm <- fitdistr(y.log, densfun = "normal")
mu.norm  <- fit.norm$estimate["mean"]
sd.norm  <- fit.norm$estimate["sd"]
logLik.norm <- fit.norm$loglik
AIC.norm <- -2 * logLik.norm + 2 * 2

# -------------------------------------------------
# 3. Fit Cauchy distribution (to log(ws30))
# -------------------------------------------------
fit.cauchy <- fitdistr(y.log, densfun = "cauchy")
loc.cauchy <- fit.cauchy$estimate["location"]
scale.cauchy <- fit.cauchy$estimate["scale"]
logLik.cauchy <- fit.cauchy$loglik
AIC.cauchy <- -2 * logLik.cauchy + 2 * 2

# -------------------------------------------------
# 4. Fit Student-t distribution (to log(ws30))
# -------------------------------------------------
nll.t <- function(pars, x) {
  mu <- pars[1]; sigma <- pars[2]; nu <- pars[3]
  if (sigma <= 0 || nu <= 0) return(Inf)
  -sum(dt((x - mu)/sigma, df = nu, log = TRUE) - log(sigma))
}

opt.t <- nlminb(start = c(mean(y.log), sd(y.log), 5),
                objective = nll.t,
                x = y.log,
                lower = c(-Inf, 1e-6, 0.5))

mu.t     <- opt.t$par[1]
sigma.t  <- opt.t$par[2]
nu.t     <- opt.t$par[3]
logLik.t <- -opt.t$objective
AIC.t    <- -2 * logLik.t + 2 * 3   # 3 parameters: μ, σ, ν

# -------------------------------------------------
# 5. QQ-Plots
# -------------------------------------------------
par(mfrow = c(1,3))

# (a) Normal QQ
q.norm <- qnorm(ppoints(n), mean = mu.norm, sd = sd.norm)
plot(q.norm, sort(y.log),
     main = "Normal Fit (log(ws30))",
     xlab = "Theoretical quantiles (Normal)",
     ylab = "Empirical quantiles (log(ws30))")
abline(0, 1, col = "blue", lwd = 2)

# (b) Cauchy QQ
q.cauchy <- qcauchy(ppoints(n), location = loc.cauchy, scale = scale.cauchy)
plot(q.cauchy, sort(y.log),
     main = "Cauchy Fit (log(ws30))",
     xlab = "Theoretical quantiles (Cauchy)",
     ylab = "Empirical quantiles (log(ws30))")
abline(0, 1, col = "red", lwd = 2)

# (c) Student-t QQ
q.t <- mu.t + sigma.t * qt(ppoints(n), df = nu.t)
plot(q.t, sort(y.log),
     main = sprintf("Student-t Fit (log(ws30), ν=%.2f)", nu.t),
     xlab = "Theoretical quantiles (t)",
     ylab = "Empirical quantiles (log(ws30))")
abline(0, 1, col = "darkgreen", lwd = 2)

par(mfrow = c(1,1))

# -------------------------------------------------
# 6. AIC Comparison
# -------------------------------------------------
cat("\n--------------------------------------------------\n")
cat("Model comparison (lower AIC = better fit):\n")
cat(sprintf("  Normal     AIC = %.2f\n",  AIC.norm))
cat(sprintf("  Cauchy     AIC = %.2f\n",  AIC.cauchy))
cat(sprintf("  Student-t  AIC = %.2f\n",  AIC.t))
cat("--------------------------------------------------\n")

# Optional: print parameter estimates
cat(sprintf("\nNormal:   μ = %.3f, σ = %.3f\n", mu.norm, sd.norm))
cat(sprintf("Cauchy:   location = %.3f, scale = %.3f\n", loc.cauchy, scale.cauchy))
cat(sprintf("Student-t: μ = %.3f, σ = %.3f, ν = %.3f\n", mu.t, sigma.t, nu.t))
