# ==============================================================================
# Assignment 2: Statistical Modeling
# ==============================================================================
# Project 1: Wind Power Forecast
# 
# Introduction
# In this project you will analyze a dataset from Tunø Knob wind power plant. 
# Wind power is the response variable, with wind speed and wind direction as 
# explanatory variables.
#
# Regression Models
# ==============================================================================

# ==============================================================================
# Load Data and transformation done in Assigment 1 for wind speed and wind Power.
# 
# (1) Wind-Power <- Box-cox.transformation()
# (2) Wind-Speed <- log.transformation()
# ==============================================================================

# Task 1.1: Formulate Initial Model
df_tuno <- read.table("data/tuno.txt", header = TRUE)
df_tuno$pow.obs <- df_tuno$pow.obs / 5000

# Transformation for wind power
box_cox_trans <- function(lambda, y){
  (1 / lambda) * log((y^lambda) / (1 - y^lambda))
}

# Transformation of Wind direction
df_tuno$wd_sin <- sin(df_tuno$wd30)
df_tuno$wd_cos <- cos(df_tuno$wd30)

# Defining variables for simplification
pow.obs <- df_tuno$pow.obs
ws30 <- df_tuno$ws30
wd_sin <- df_tuno$wd_sin
wd_cos <- df_tuno$wd_cos 
wd_rad <- df_tuno$wd_rad

# Power transformation from assigment 1
lambda.eq1 <- 0.326
boxpow <- box_cox_trans(lambda.eq1, pow.obs)

# Log transformation for wind speed
logws <- log(ws30)

# Plotting
par(bg = "white")
plot(pow.obs, ws30, main = "Power vs Wind Speed", xlab = "Wind Speed (m/s)", ylab = "Power Output")
head(df_tuno)

# Plot transformed data
par(bg = "white")
plot(boxpow, logws, main = "Transformed Power vs Log Wind Speed", xlab = "Log Wind Speed", ylab = "Box-Cox Transformed Power")

# ==============================================================================
# Base model
# ==============================================================================

# Used later. For plotting
plot_model <- function(model, title = "Model Diagnostics") {
  width <- 6
  height <- 8
  par(bg = "white")
  options(repr.plot.width = width, repr.plot.height = height)
  print(summary(model))
  plot(model, which = 1, main = title)
}

# Base model: pow.obs ~ ws30 + ws30^2
basemodel <- lm(pow.obs ~ ws30 + I(ws30^2))
plot_model(basemodel, "Base Model: Power ~ Wind Speed")

# ==============================================================================
# Question 1.2
# You might consider non-normal models and/or normal model with data 
# transformation. Further you might consider including wind direction. 
# You should develop a suited model for prediction of daily power production.
# ==============================================================================

# Task 1.2: Model Development with Transformations
# Transformed model: boxpow ~ logws + logws^2
Box_cox_model <- lm(boxpow ~ logws + I(logws^2))
plot_model(Box_cox_model, "Box-Cox Model: Transformed Power ~ Log Wind Speed")

# Box-Cox model with wind direction
Box_cox_model_wd <- lm(boxpow ~ logws + I(logws^2) + wd_sin + wd_cos)
plot_model(Box_cox_model_wd, "Transformed Power ~ Log Wind Speed + Wind Direction")

# Gamma Generalized Linear Model (GLM)
gamma_model <- glm(pow.obs ~ ws30 + I(ws30^2), family = Gamma)
plot_model(gamma_model, "Gamma GLM Model: Power ~ Wind Speed")

# ==============================================================================
# Question 1.3
# Present the parameters of the final model, this include quantification of 
# the uncertainty of the parameters.
# ==============================================================================

# Extract coefficients and standard errors
beta <- coef(Box_cox_model_wd)
se <- summary(Box_cox_model_wd)$coefficients[, "Std. Error"]

# Calculate Lower and Upper bounds manually
lower <- beta - 1.96 * se
upper <- beta + 1.96 * se

# Display
cbind(Estimate = beta, SE = se, Lower = lower, Upper = upper)

# ==============================================================================
# Question 1.4
# Give an interpretation of the parameters in particular this should include 
# presentation of any nonlinear functions (series expansions) of the 
# explanatory variables.
# ==============================================================================

# ==============================================================================
# Question 1.5
# Present the final model, e.g. some graphical presentation of predictions 
# under different scenarios of wind speed and wind direction.
# ==============================================================================

# 1. Setup
lambda <- 0.326 

# 2. Plot the raw data
par(bg = "white") 
width <- 16
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
plot(df_tuno$ws30, df_tuno$pow.obs, 
     pch = 16, cex = 0.9, col = "darkgrey",
     xlab = "Wind Speed (m/s)", 
     ylab = "Power Output", 
     main = "Final Model Fit")

# 3. Create prediction sequence
ws_seq <- seq(0, 25, length.out = 200)

# Keep cos and sin konstant at their mean value 
new_data <- data.frame(
  logws = log(ifelse(ws_seq == 0, 1e-4, ws_seq)), 
  wd_sin = mean(df_tuno$wd_sin),
  wd_cos = mean(df_tuno$wd_cos)
)

# 4. Predict and Back-Transform (Using the LOGIT-POWER Inverse)
# Step A: Get the model prediction (z)
z <- predict(Box_cox_model_wd, newdata = new_data)

# Step B: Apply the specific inverse formula for your transformation
# inverse formula USING SYMPY CHECK INVERSE.IPYNB
y_pred <- (exp(lambda * z) / (1 + exp(lambda * z)))^(1/lambda)

# Plot the line
lines(ws_seq, y_pred, col = "red", lwd = 2)
