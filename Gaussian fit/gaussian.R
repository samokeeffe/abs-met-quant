# Define data
df <- data.frame("ratio" = c(1, 4, 20, 0.25, 0.05), 
                 "R^2" = c(0.8453953, 0.7848674, 0.6734205, 0.7848674, 0.6734205))

# Log-transform the ratio
df$ln_ratio <- log(df$ratio)

# Extract x and y
x <- df$ln_ratio
y <- df$R.2

#### FITTING
# Non-linear least squares regression to fit to a Gaussian function 
starting_width <- sd(x)
nlsout_test <- nls(y ~ a*exp(-((x)^2)/(2*c^2)), start = c(a = 0.83, c = starting_width))
summary(nlsout_test)


# Verifying output from fitting and calculating R^2
a <- 0.831451 
c <- 4.575862

gaussian <- function(x) {
  a*exp(-(x^2)/(2*c^2))
}

y_gaussian <- gaussian(x)

#R^2 calculation
residuals <- y - y_gaussian
SSR <- sum(residuals^2)
SST <- sum((df$R.2 - mean(df$R.2))^2)
print(SST)
print(SSR)
R2 <- 1-SSR/SST
print(R2)


