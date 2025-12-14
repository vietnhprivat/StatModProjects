###############################################################
# Fit von Mises using PROFILE LIKELIHOOD
# Produce: histogram + fitted density curve
###############################################################

df_tuno <- read.table("tuno.txt", header = TRUE)
theta <- df_tuno$wd30   # wind direction (radians)

###############################################################
# von Mises negative log-likelihood (optimize only over kappa)
###############################################################
nll_kappa <- function(kappa, mu, x) {
  if (kappa <= 0) return(Inf)
  -sum( kappa*cos(x - mu) - log(2*pi) - log(besselI(kappa, 0)) )
}

###############################################################
# Profile log likelihood for μ
###############################################################
lp_mu <- function(mu, x) {
  opt <- optimize(
    f = nll_kappa,
    mu = mu,
    x = x,
    interval = c(0.0001, 50)
  )
  return(opt$objective)
}

###############################################################
# Compute μ̂ via profiling
###############################################################
mu_grid <- seq(-pi, pi, length.out = 400)
logLp <- sapply(mu_grid, lp_mu, x = theta)
mu_hat <- mu_grid[which.min(logLp)]   # best μ

###############################################################
# Compute κ̂ at μ̂
###############################################################
opt_kappa <- optimize(
  f = nll_kappa,
  mu = mu_hat,
  x = theta,
  interval = c(0.0001, 50)
)
kappa_hat <- opt_kappa$minimum

cat("Estimated mu:", mu_hat, "\n")
cat("Estimated kappa:", kappa_hat, "\n")

###############################################################
# von Mises PDF for plotting
###############################################################
vonmises_pdf <- function(x, mu, kappa) {
  exp(kappa*cos(x - mu)) / (2*pi*besselI(kappa, 0))
}

###############################################################
# Plot histogram + fitted von Mises density
###############################################################
hist(theta,
     breaks = 40,
     freq = FALSE,
     col = "lightgray",
     main = "Wind Direction with Fitted von Mises Model",
     xlab = "Wind direction (radians)")

curve(vonmises_pdf(x, mu_hat, kappa_hat),
      add = TRUE,
      lwd = 3,
      col = "red")
