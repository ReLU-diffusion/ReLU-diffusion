library(nleqslv)

# Choose c value here
c_val <- 1  # or 0.5 for the other case

# Simulation function: theta0 ~ Unif(0, c * sqrt(theta2))
sim_relu <- function(mu, theta1, theta2, n_trials = 10, dt = 0.001, c = c_val) {
  T_vals <- numeric(n_trials)
  for (i in 1:n_trials) {
    theta0 <- runif(1, 0, c * sqrt(theta2))
    rho <- mu - theta1 / 2
    t <- 0
    xi <- theta0
    m <- min(0, xi)
    while ((xi - m) < theta2) {
      t <- t + dt
      xi <- xi + rho * dt + sqrt(dt) * rnorm(1)
      m <- min(m, xi)
      if (t > 1000) break
    }
    T_vals[i] <- t
  }
  return(T_vals)
}

# Updated expected_T using the formula under theta0 ~ Unif(0, c * sqrt(theta2))
expected_T <- function(rho, theta2, c = c_val) {
  term1 <- exp(-2 * rho * theta2)
  term2 <- (1 - exp(-2 * rho * c * sqrt(theta2))) / (2 * rho * c * sqrt(theta2))
  term3 <- 2 * rho * (theta2 - (c / 2) * sqrt(theta2))
  (1 / (2 * rho^2)) * (term1 - term2 + term3)
}

# Estimation using numerical solver
estimate_theta <- function(T0_bar, T1_bar, c = c_val) {
  theta1_approx <- sqrt(log(T0_bar) / T1_bar)
  theta2_approx <- sqrt(T1_bar * log(T0_bar))
  initial_guess <- c(theta1_approx, theta2_approx)
  
  obj_fn <- function(par) {
    theta1 <- par[1]
    theta2 <- par[2]
    rho0 <- -theta1 / 2
    rho1 <- theta1 / 2
    val0 <- expected_T(rho0, theta2, c)
    val1 <- expected_T(rho1, theta2, c)
    return(c(val0 - T0_bar, val1 - T1_bar))
  }
  
  result <- nleqslv(x = initial_guess, fn = obj_fn)
  return(result$x)
}

# Main simulation loop
set.seed(123)
M <- 100
n_per_group <- 10
true_theta1 <- 1
true_theta2 <- 5

results <- data.frame(
  theta1_hat = numeric(M),
  theta2_hat = numeric(M),
  T0_mean = numeric(M),
  T1_mean = numeric(M),
  T0_hat = numeric(M),
  T1_hat = numeric(M)
)

for (m in 1:M) {
  T0_vals <- sim_relu(mu = 0, theta1 = true_theta1, theta2 = true_theta2, n_trials = n_per_group, c = c_val)
  T1_vals <- sim_relu(mu = true_theta1, theta1 = true_theta1, theta2 = true_theta2, n_trials = n_per_group, c = c_val)
  
  T0_bar <- mean(T0_vals)
  T1_bar <- mean(T1_vals)
  
  theta_hat <- estimate_theta(T0_bar, T1_bar, c = c_val)
  theta1_hat <- theta_hat[1]
  theta2_hat <- theta_hat[2]
  
  rho0_hat <- -theta1_hat / 2
  rho1_hat <- theta1_hat / 2
  
  T0_hat <- expected_T(rho0_hat, theta2_hat, c = c_val)
  T1_hat <- expected_T(rho1_hat, theta2_hat, c = c_val)
  
  results[m, ] <- c(theta1_hat, theta2_hat, T0_bar, T1_bar, T0_hat, T1_hat)
}

# Compute summary statistics
theta1_mean <- mean(results$theta1_hat)
theta2_mean <- mean(results$theta2_hat)
theta1_sd <- sd(results$theta1_hat)
theta2_sd <- sd(results$theta2_hat)
theta1_rmse <- sqrt(mean((results$theta1_hat - true_theta1)^2))
theta2_rmse <- sqrt(mean((results$theta2_hat - true_theta2)^2))

# Present all summary results
summary_df <- data.frame(
  theta1_mean = theta1_mean,
  theta2_mean = theta2_mean,
  theta1_sd = theta1_sd,
  theta2_sd = theta2_sd,
  theta1_rmse = theta1_rmse,
  theta2_rmse = theta2_rmse
)

print(summary_df)

#n=10
#theta1_mean theta2_mean theta1_sd theta2_sd theta1_rmse theta2_rmse
#1    1.009291    5.014361 0.1199388 0.3923618   0.1196988   0.3906591

#n=100
#theta1_mean theta2_mean  theta1_sd theta2_sd theta1_rmse theta2_rmse
#1   0.9929312    5.025567 0.04041266  0.124954   0.0408267   0.1269294

