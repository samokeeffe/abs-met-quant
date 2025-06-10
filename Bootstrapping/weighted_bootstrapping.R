# Load necessary libraries
{
  library(readxl)
  library(openxlsx)
  library(tidyr)
  library(dplyr)
}

# Import data
data <- readRDS("Bootstrapping/data/Jurkat_dttp.rds")

# Define dataframes, matrices and variables - ** be sure to only include Jurkat, CD4+ or CD8+ data **
{met <- data[,1]
  ln_conc <- as.matrix(cbind(data[,2:9])) # set columns 
  ratios <- as.matrix(cbind(data[,11:18])) # set columns ratios are labeled "RXX" where XX is the sample label
  nmet <- nrow(met)
  ecoli_error <- data[["ecoli sem ln"]] # this column should be present in every data file
}

# FUNCTIONS
#Gaussian fit based on non linear least squares fitting
gaussian_so <- function(x) {
  # parameters from fitting
  a <- 0.831451 
  c <- 4.575862
  
  #function
  a*exp(-((x)^2)/(2*c^2))
}

# Define parabolic likeliness function which converts a given ratio to a R^2 value
likelihood <- function (x)
{
  ifelse(!(is.na(x)), gaussian_so(x), 0) #not sure this is the best way to do this # took out ! is.infinite
}

  


# create a data frame to store the ci and sem
{ci_results <- data.frame(Metabolite = met, Lower.boot = NA, Upper.boot = NA, Mean.boot = NA, SEM = NA, error_prop = NA, LB.ln = NA, UB.ln = NA, Geommean = NA, LB = NA, UB = NA)
  
# Define empty vectors
#ln_met <- vector(length = nmet)
confidence_interval <- matrix(nrow = nmet, ncol = 2)
stdev <- matrix(nrow = nmet, ncol = 1)
}

# Convert ratios to R^2 values
R2_new <- likelihood(log(ratios)) 
  
# Define probability matrix
probability <- R2_new/rowSums(R2_new)
  
# Verify the probabilities for each metabolite is 1 and everything >0
for (x in 1:nmet) {
  
  # Skip rows that contain NA or NaN values
  if (is.na(rowSums(probability)[x]) || any(is.nan(rowSums(probability)[x]))) {
    next
  }
  
  # Check if row sum is approximately equal to 1
  if (!all.equal(rowSums(probability)[x], 1)) { 
    print(x)
    print(rowSums(probability)[x])
    cat("Row", x, "probability not equal to 1.\n")
  }
  
  # Check for negative probabilities
  for (y in 1:ncol(probability)) {
    if (!is.na(probability[x, y]) && !is.nan(probability[x, y])) { 
      if (probability[x, y] < 0) {
        cat("Probability in position", x, y, "is less than zero\n")
      }
    }
  }
}

  
# Bootstrap each metabolite and calculate means and ci
for (x in 1:nmet){
    # Count sample size for that metabolite
    sample_size = sum(probability[x,] != 0 & !is.na(probability[x,]))
    nrep = 10000
    
    if (is.na(rowSums(probability)[x]) || any(is.nan(rowSums(probability)[x]))) {
      next
    }
    means <- vector(length=nrep)
    
    for (i in 1:nrep){
    # sample with weights with replacing 1000 times with probabilities
    bootstrap_sample <- sample(ln_conc[x,], size = sample_size, replace = TRUE, prob = probability[x,])
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
    ci_results$SEM[x] <- (ci_results$Upper.boot[x] - ci_results$Lower.boot[x])/4  #compare to sd 
    cat("The SEM is", ci_results$SEM[x], "\n")
    
    #Error propagation
    # Validate SEM and ecoli_error before computing error_prop 
    if (!is.na(ecoli_error[x]) && ci_results$SEM[x]^2 > ecoli_error[x]^2) {
      ci_results$error_prop[x] <- sqrt((ci_results$SEM[x]^2) - ecoli_error[x]^2)
       } else {
      ci_results$error_prop[x] <- ecoli_error[x]  
       }
    
    #Calculate adjusted upper and lower bounds and translate to normal space
    ci_results$LB.ln[x] <- ci_results$Mean.boot[x] - 2*ci_results$error_prop[x]
    ci_results$UB.ln[x] <- ci_results$Mean.boot[x] + 2*ci_results$error_prop[x]
    ci_results$Geommean[x] <- exp(ci_results$Mean.boot[x])
    ci_results$LB[x] <- exp(ci_results$LB.ln[x])
    ci_results$UB[x] <- exp(ci_results$UB.ln[x])
}
  