##################################################
## Assignment 1 – Transformations of Wind Speed (ws30)
## Following Lecture 4 (Example 4.16)
##################################################

rm(list = ls())

# -------------------------------------------------
# 1. Load data
# -------------------------------------------------
df_tuno <- read.table("tuno.txt", header = TRUE)
y <- df_tuno$ws30
y <- y[y > 0]        # Box-Cox requires strictly positive data
n <- length(y)

# -------------------------------------------------
# 2. Define Box-Cox transformation
# -------------------------------------------------
bc.trans <- function(lambda, y){
  if (lambda == 0) return(log(y))
  (y^lambda - 1) / lambda
}

# -------------------------------------------------
# 3. Find optimal lambda (profile likelihood)
# -------------------------------------------------
lp.lambda <- function(lambda, y){
  n <- length(y)
  y.l <- bc.trans(lambda, y)
  sigma2 <- mean((y.l - mean(y.l))^2)
  -n/2 * log(sigma2) + (lambda - 1) * sum(log(y))
}

lambda.grid <- seq(-1, 1, by = 0.01)
lp <- sapply(lambda.grid, lp.lambda, y = y)
lambda.hat <- lambda.grid[which.max(lp)]
lambda.hat  # print estimated lambda

# -------------------------------------------------
# 4. Transformations
# -------------------------------------------------
y.log  <- log(y)
y.box  <- bc.trans(lambda.hat, y)

# -------------------------------------------------
# 5. Histograms
# -------------------------------------------------
par(mfrow = c(1, 3))

hist(y, breaks = 30, col = "gray90", main = "No transformation",
     xlab = "Wind Speed (m/s)")
hist(y.log, breaks = 30, col = "lightblue", main = "Log-transformed",
     xlab = "log(Wind Speed)")
hist(y.box, breaks = 30, col = "lightgreen",
     main = sprintf("Box–Cox (λ = %.2f)", lambda.hat),
     xlab = "Box–Cox transformed speed")

par(mfrow = c(1,1))

# -------------------------------------------------
# 6. QQ-plots
# -------------------------------------------------
par(mfrow = c(1,3))

qqnorm(y, main = "QQ-Plot: No transformation")
qqline(y, col = "blue", lwd = 2)

qqnorm(y.log, main = "QQ-Plot: Log-transformed")
qqline(y.log, col = "blue", lwd = 2)

qqnorm(y.box, main = sprintf("QQ-Plot: Box–Cox (λ = %.2f)", lambda.hat))
qqline(y.box, col = "blue", lwd = 2)

par(mfrow = c(1,1))

# -------------------------------------------------
# 7. Profile likelihood plot for λ (optional)
# -------------------------------------------------
plot(lambda.grid, lp - max(lp), type = "l",
     main = "Profile log-likelihood for λ (Box–Cox)",
     xlab = "λ", ylab = "log-likelihood (centered)")
abline(v = lambda.hat, col = "red", lwd = 2, lty = 2)
lines(range(lambda.grid), -0.5 * qchisq(0.95, df = 1) * c(1,1),
      col = "darkgreen", lty = 2)

