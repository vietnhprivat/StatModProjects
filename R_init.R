# R Initialization Script for Jupyter
# Install and load required packages for Jupyter R kernel

# Install IRkernel if not already installed
if (!require("IRkernel", quietly = TRUE)) {
  install.packages("IRkernel", repos = "https://cloud.r-project.org/")
  library(IRkernel)
}

# Register R kernel with Jupyter
IRkernel::installspec(user = TRUE)
