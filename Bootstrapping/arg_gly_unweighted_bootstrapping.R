# Load necessary libraries
{
  library(readxl)
  library(openxlsx)
  library(tidyr)
  library(dplyr)
}

# Load processed data
data <- readRDS("Bootstrapping/data/donorF_arg_gly.rds")

# Define dataframes, matrices and variables 
# be sure to only include CD4+ or CD8+ data columns 2:4 are CD4+, 5:7 are CD8+

{met <- data[,1]
  ln_conc <- as.matrix(cbind(data[,5:7]))
  nmet <- nrow(met)
}


# create a data frame to store the ci and sem
{ci_results <- data.frame(Metabolite = met, Lower.boot = NA, Upper.boot = NA, Mean.boot = NA, SEM = NA, Geommean = NA, LB = NA, UB = NA)
  
  # Define empty vectors
  #ln_met <- vector(length = nmet)
  confidence_interval <- matrix(nrow = nmet, ncol = 2)
  stdev <- matrix(nrow = nmet, ncol = 1)
}

# Bootstrap each metabolite and calculate means and ci
for (x in 1:nmet){
  # Count sample size for that metabolite
  sample_size = 3
  nrep = 10000

  means <- vector(length=nrep)
  
  for (i in 1:nrep){
    # sample with weights with replacing 1000 times with probabilities
    bootstrap_sample <- sample(ln_conc[x,], size = sample_size, replace = TRUE)
    means[i] <- mean(bootstrap_sample)
  }
  
  # Bootstrapped 95% ci
  confidence_interval[x,] <- quantile(means,probs=c(0.025,0.975)) # add standard deviation calculation
  cat("The 95% confidence interval is", confidence_interval[x,], "\n")
  stdev[x,] <- sd(means)
  cat("The sd is ", stdev[x,], "\n") # Verify that the calculated standard deviation is approximately equal to the measured sd from the 95% ci
  
  #Store results
  ci_results$Lower.boot[x] <- confidence_interval[x,1]
  ci_results$Upper.boot[x] <- confidence_interval[x,2]
  ci_results$Mean.boot[x] <- mean(means)
  ci_results$SEM[x] <- (ci_results$Upper.boot[x] - ci_results$Lower.boot[x])/4  
  cat("The SEM is", ci_results$SEM[x], "\n")
  
  #Calculate adjusted upper and lower bounds and translate to normal space
  ci_results$Geommean[x] <- exp(ci_results$Mean.boot[x])
  ci_results$LB[x] <- exp(ci_results$Lower.boot[x])
  ci_results$UB[x] <- exp(ci_results$Upper.boot[x])
}
