

data_finance <- read.table("data/finance_data.csv", header = TRUE, sep = ";")

# Udtræk kolonnen USO
x <- as.numeric(data_finance$USO)

summary(x)
hist(x, breaks = 30, main = "Histogram of data", xlab = "x", probability = TRUE)
lines(density(x), col = "blue", lwd = 2)

######### præsenter data
hist(x, breaks = 30, probability = TRUE, 
     main = "Histogram of weekly returns (USO)", 
     xlab = "Weekly return (USO)")
lines(density(x), col = "blue", lwd = 2)
legend("topright", legend = c("Normal"), col = c("blue"), lwd = 2)

# 3. Estimer parametre for en normalmodel
mu_hat <- mean(x)
sigma_hat <- sd(x)
mu_hat; sigma_hat

# 4. Tjek normalitet via et QQ-plot
qqnorm(x, main = "QQ-plot (normal model)",)
qqline(x, col = "red")

#Vi kan se tunge haler ved normalfordelingen

########## Derfor kunne vi bruge en t-fordeling (god hvis der er tunge haler)
#### t-fordelingen håndterer ektreme afkast bedre end normalfordelingen,
#### som gør at fordelingen håndterer tunge haler bedre
library(MASS)

# Fit en t-fordeling (maksimum likelihood-estimering (MLE) via fitdistr() fra pakken MASS):
fit_t <- fitdistr(x, "t") 
fit_t

m_hat  <- fit_t$estimate["m"]
s_hat  <- fit_t$estimate["s"]
df_hat <- fit_t$estimate["df"]


### Tjekker fit via et QQ-plot
# Gem estimerede parametre
m_hat  <- fit_t$estimate["m"] # middelværdi
s_hat  <- fit_t$estimate["s"] # skala (standardafvigelse)
df_hat <- fit_t$estimate["df"] # frihedsgrader

# QQ plot mod t-fordeling
qqplot(
  qt(ppoints(length(x)), df = df_hat),      # teoretiske kvantiler
  (x - m_hat) / s_hat,                      # standardiserede data
  main = "QQ-plot (t-distribution)",
  xlab = "Theoretical quantiles",
  ylab = "Sample quantiles"
)
abline(0, 1, col = "red", lwd = 2)

### Dette QQ-plot passer bedre

# Sammenlign fit (normalfordeling mod t-fordeling)
hist(x, breaks = 30, probability = TRUE, main = "Normal vs t-distribution", xlab = "x")
curve(dnorm(x, mean = mu_hat, sd = sigma_hat), col = "red", add = TRUE, lwd = 2)
curve(dt((x - fit_t$estimate["m"])/fit_t$estimate["s"], df = fit_t$estimate["df"]) / fit_t$estimate["s"],
      col = "blue", add = TRUE, lwd = 2)
legend("topright", legend = c("Normal", "t-distribution"), col = c("red", "blue"), lwd = 2)


############# Præsenterer endelig model


cat("Normalmodel:\n")
cat("  μ =", mu_hat, "\n")
cat("  σ =", sigma_hat, "\n\n")

cat("t-fordeling:\n")
print(fit_t)

### Her defineres log likelihoods hhv. for normalfordelingen og t-fordelingen
logLik_norm <- sum(dnorm(x, mu_hat, sigma_hat, log = TRUE))
logLik_t <- fit_t$loglik

### Sammenligner med AIC:
AIC_norm <- -2 * logLik_norm + 2 * 2  # 2 parametre
AIC_t <- -2 * logLik_t + 2 * 3         # 3 parametre (mu, sigma, df)

cat("AIC (Normal):", AIC_norm, "\n")
cat("AIC (t):", AIC_t, "\n")


### Vi får meget negative AIC-værdier og det fordi log-likelihood værdierne er meget positive 
### (dvs modellerne passer dataene rigtigt godt). (se formel for AIC)

### T-fordeling er 40 mindre så den er signifikant bedre fit til vores finance data. 

### Generelt er fit af fordeling til model bedre jo lavere AIC er. 









