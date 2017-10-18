# Adam Tonks, Yale-NUS College

### PART A ###
# load package for Inverse-Gamma
library(LearnBayes)

# setup parameters of distributions to simulate from
x <- c(10, 11, 12, 11, 9)
y_bar <- mean(x)
n <- length(x)
S <- sum((x-y_bar)^2)

# generate the samples
set.seed(412)
sigma_sq <- rigamma(10000, n/2, S/2)
set.seed(412)
mu_samp1 <- rnorm(10000, y_bar, sqrt(sigma_sq/n))

# find mean and variance
mean(mu_samp1)
var(mu_samp1)

### PART D ###
# function for computing posterior joint PDF
posterior_p <- function(mu, lambda) {
  F_X1 <- pnorm(x + 0.5, mu, exp(lambda/2))
  F_X2 <- pnorm(x - 0.5, mu, exp(lambda/2))
  return(prod(F_X1-F_X2))
}

# setup the grid values
mu_grid <- seq(5, 16, 0.01)
lambda_grid <- seq(-3, 6, 0.01)
probs <- rep(NA, length(mu_grid) * length(lambda_grid))

# populate the points for discrete approximation
for(i in seq_along(mu_grid)) {
  for(j in seq_along(lambda_grid)) {
    probs[(i-1)*length(lambda_grid) + j] <- posterior_p(mu_grid[i], lambda_grid[j]) * 0.01^2 / 0.001682444
  }
}

# put points into a grid
probs_mat <- matrix(probs, length(lambda_grid), length(mu_grid))

# examine the shape of the PDF, to check that we sample from correct region
op <- par(no.readonly=TRUE)
par(mar=c(5, 5, 5, 5))
filled.contour(lambda_grid, mu_grid, probs_mat, main="Posterior joint density",
               xlab="lambda", ylab="mu")
par(op)

# setup the lambda and mu values to be sampled
lambda_mat <- rep(lambda_grid, length(mu_grid))
mu_mat <- rep(mu_grid, each = length(lambda_grid))

# sample from points on discrete approximation according to the probabilites
set.seed(432)
theta_index <- sample(seq_len(length(lambda_grid)*length(mu_grid)), 10000, replace=TRUE, prob=probs)

# plot the samples onto the contour plot to check that we did not exclude regions where draws are plausible
par(mar=c(5, 5, 5, 5))
filled.contour(lambda_grid, mu_grid, probs_mat, main="Posterior joint density",
               xlab="lambda", ylab="mu", plot.axes={points(lambda_mat[theta_index], mu_mat[theta_index])})
par(op)

# find mean and variance
mean(mu_samp2)
var(mu_samp2)

### PART E ###
# compare the means and variances
(mean(mu_samp2) - mean(mu_samp1)) / mean(mu_samp1)
(var(mu_samp2) - var(mu_samp1)) / var(mu_samp1)
