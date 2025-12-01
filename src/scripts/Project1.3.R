# ==============================================================================
# Assignment 3: Statistical Modeling
# ==============================================================================
# Project 1: Wind Power Forecast
# ==============================================================================

# ==============================================================================
# Load Data and transformation done in Assigment 1 for wind speed and wind Power.
# 
# (1) Wind-Power <- Box-cox.transformation()
# (2) Wind-Speed <- log.transformation()
# ==============================================================================

# Task 1.1: Load Data and Define Transformations
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

# Display data summary
head(df_tuno)

# ==============================================================================
# Final model from assignment 2
# ==============================================================================

# Box-Cox model with wind direction
Box_cox_model_wd <- lm(boxpow ~ logws + I(logws^2) + wd_sin + wd_cos)
summary(Box_cox_model_wd)
