rm(list = ls())
set.seed(1)

#--------#
# prepare
#--------#
library(foreach)
library(magrittr)
library(codetools)
library(doParallel)
library(HighBLP)

#--------------#
# set constants
#--------------#
# number of markets
T <- 100
# total number of products
J <- 10
# number of consumers for simulation
N <- 1000
# number of consumers for estimation
M <- 1000
# number of product characteristics excluding p
K <- 5

#---------------------#
# set hyper parameters
#---------------------#
# probability that a product is available
prob_availability <- 0.7
# partial correlation between p and xi
cor_xi_p <- 0.6
# partial correlation between p and marginal cost shocks
cor_mc_p <- 0.5

#-----------------------------#
# simulate exogenosu variables
#-----------------------------#
# product characteristics
X <- matrix(rnorm(J * K), nrow = J)
X[, 1] <- 1
# availability of the products
availability <-
  foreach (t = 1:T) %do% {
    availability_t <- purrr::rbernoulli(J, p = prob_availability)
    availability_t <- which(availability_t)
    availability_t <- sort(availability_t)
    return(availability_t)
  }
# product characteristics across markets
X <-
  foreach (t = 1:T) %do% {
    X_t <- X[availability[[t]], , drop = FALSE]
    return(X_t)
  }
# unobserved fixed effects
xi <-
  foreach (t = 1:T) %do% {
    xi_t <- matrix(rnorm(J))
    xi_t <- xi_t[availability[[t]], , drop = FALSE]
    return(xi_t)
  }
# idiosyncratic shock
epsilon <-
  foreach (t = 1:T) %do% {
    epsilon_t <- matrix(evd::rgev(J * N), nrow = J)
    epsilon_t <- epsilon_t[availability[[t]], , drop = FALSE]
    return(epsilon_t)
  }
# taste shocks on X
nu <-
  foreach (t = 1:T) %do% {
    nu_t <- matrix(rnorm(K * N), nrow = K)
    return(nu_t)
  }
# taste shocks on p
upsilon <-
  foreach (t = 1:T) %do% {
    upsilon_t <- matrix(rnorm(N), nrow = 1)
    return(upsilon_t)
  }
# marginal cost shocks
mc <-
  foreach (t = 1:T) %do% {
    mc_t <- matrix(rnorm(J), nrow = J)
    mc_t <- mc_t[availability[[t]], , drop = FALSE]
    return(mc_t)
  }
# p
p <- 
  foreach (t = 1:T) %do% {
    p_t <- matrix(exp(rnorm(length(availability[[t]])) + cor_xi_p * xi[[t]] + cor_mc_p * mc[[t]]))

    return(p_t)
  }
# use mc as instrumental variables
Z <- mc

#---------------#
# set parameters
#---------------#
# mean taste for money
alpha <- - 0.5
# mean taste for X
beta <- matrix(rnorm(K))
# standard deviation of taste for money
sigma_upsilon <- 0.1
# standard deviation of taste for X
sigma_nu <- exp(rnorm(K))
# bundle parameters
theta_linear <- matrix(c(beta, alpha))
theta_nonlinear <- matrix(c(sigma_nu, sigma_upsilon))
# random l
rl <- as.integer(theta_nonlinear > 0)

#-------------------------------#
# simulate endogenous vavriables
#-------------------------------#
# compute mean utility
mean_utility <- compute_mean_utility(beta, alpha, X, p, xi)
mean_utility_rcpp <- compute_mean_utility_rcpp(beta, alpha, X, p, xi)
max(abs(unlist(mean_utility) - unlist(mean_utility_rcpp)))

# compute individual utility
individual_utility <-
  compute_individual_utility(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
individual_utility_rcpp <-
  compute_individual_utility_rcpp(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
max(abs(unlist(individual_utility) - unlist(individual_utility_rcpp)))

# compute individual share
individual_share <-
  compute_individual_share(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
individual_share_rcpp <-
  compute_individual_share_rcpp(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
max(abs(unlist(individual_share) - unlist(individual_share_rcpp)))

# compute share
share <- compute_share(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
share_rcpp <- compute_share_rcpp(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
max(abs(unlist(share) - unlist(share_rcpp)))

# compute invidual utility from delta
individual_utility_delta <- 
  compute_individual_utility_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
individual_utility_delta_rcpp <- 
  compute_individual_utility_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
max(abs(unlist(individual_utility) - unlist(individual_utility_delta)))
max(abs(unlist(individual_utility_delta) - unlist(individual_utility_delta_rcpp)))

# compute individual share from delta
individual_share_delta <- 
  compute_individual_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
individual_share_delta_rcpp <- 
  compute_individual_share_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
max(abs(unlist(individual_share) - unlist(individual_share_delta)))
max(abs(unlist(individual_share_delta) - unlist(individual_share_delta_rcpp)))

# compute share
share_delta <- 
  compute_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
share_delta_rcpp <- 
  compute_share_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
max(abs(unlist(share_delta) - unlist(share)))
max(abs(unlist(share_delta) - unlist(share_delta_rcpp)))

# invert share
mean_utility <- invert_share(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
mean_utility_rcpp <- invert_share_rcpp(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
max(abs(unlist(mean_utility) - unlist(mean_utility_rcpp)))

# make initial weight
X_vec <- 
  X %>%
  purrr::reduce(rbind)
X_vec_rcpp <- vstack_rcpp(X)
max(abs(X_vec - X_vec_rcpp))
Z_vec <-
  Z %>%
  purrr::reduce(rbind)
XZ <- cbind(X_vec, Z_vec)
XZ <- cbind(XZ, XZ[, 2:ncol(XZ)]^2, XZ[, 2:ncol(XZ)]^3)
W <- crossprod(XZ, XZ)
W <- diag(ncol(XZ))

# estimate the linear parameters
theta_linear_hat <- estimate_linear_parameters(mean_utility, X, p, Z, W)
theta_linear_hat_rcpp <- estimate_linear_parameters_rcpp(mean_utility, X, p, Z, W)
cbind(theta_linear_hat, theta_linear, theta_linear_hat - theta_linear)
max(abs(theta_linear_hat - theta_linear))
max(abs(theta_linear_hat - theta_linear_hat_rcpp))

# elicit xi
xi_hat <- elicit_xi(theta_linear_hat, mean_utility, X, p)
xi_hat_rcpp <- elicit_xi_rcpp(theta_linear_hat, mean_utility, X, p)
max(abs(unlist(xi_hat) - unlist(xi_hat_rcpp)))

# compute moments
moments <- compute_moments(theta_nonlinear, share, mean_utility, X, p, Z, nu, upsilon, W)

# compute objective function
objective <- compute_objective(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W)

# compute derivatives of share with respect to delta
share_derivatives_wrt_mean_utility <- compute_share_derivatives_wrt_mean_utility(individual_share)
share_derivatives_wrt_mean_utility_numDeriv <- compute_share_derivatives_wrt_mean_utility_numDeriv(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
share_derivatives_wrt_mean_utility_rcpp <- compute_share_derivatives_wrt_mean_utility_rcpp(individual_share)
max(abs(unlist(share_derivatives_wrt_mean_utility) - unlist(share_derivatives_wrt_mean_utility_numDeriv)))
max(abs(unlist(share_derivatives_wrt_mean_utility) - unlist(share_derivatives_wrt_mean_utility_rcpp)))

# compute derivatives of share with respect to non-linear parameters
share_derivatives_wrt_theta_nonlinear <- compute_share_derivatives_wrt_theta_nonlinear(individual_share, X, p, nu, upsilon)
share_derivatives_wrt_theta_nonlinear_numDeriv <- compute_share_derivatives_wrt_theta_nonlinear_numDeriv(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
share_derivatives_wrt_theta_nonlinear_rcpp <- compute_share_derivatives_wrt_theta_nonlinear_rcpp(individual_share, X, p, nu, upsilon)
max(abs(unlist(share_derivatives_wrt_theta_nonlinear) - unlist(share_derivatives_wrt_theta_nonlinear_numDeriv)))
max(abs(unlist(share_derivatives_wrt_theta_nonlinear) - unlist(share_derivatives_wrt_theta_nonlinear_rcpp)))

# compute derivatives of mean utility with respect to non-linear parameters
mean_utility_derivatives_wrt_theta_nonlinear <- compute_mean_utility_derivatives_wrt_theta_nonlinear(individual_share, X, p, nu, upsilon)
mean_utility_derivatives_wrt_theta_nonlinear_rcpp <- compute_mean_utility_derivatives_wrt_theta_nonlinear_rcpp(individual_share, X, p, nu, upsilon)
max(abs(unlist(mean_utility_derivatives_wrt_theta_nonlinear) - unlist(mean_utility_derivatives_wrt_theta_nonlinear_rcpp)))

# compute derivatives of the objective function with respect to non-linear parameters
objective_derivatives_wrt_theta_nonlinear <-
  compute_objective_derivatives_wrt_theta_nonlinear(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W)
objective_derivatives_wrt_theta_nonlinear_numDeriv <-
  compute_objective_derivatives_wrt_theta_nonlinear_numDeriv(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W)
max(abs(unlist(objective_derivatives_wrt_theta_nonlinear) - unlist(objective_derivatives_wrt_theta_nonlinear_numDeriv)))

# estiamte parameters
solution <- estimate_parameters(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W) 
theta_nonlinear_hat <- solution$par
max(abs(theta_nonlinear_hat - theta_nonlinear))

# compute efficient weighting matrix
W_efficient <- compute_efficient_weighting_matrix(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W)

# compute standard errors
covariance_theta <- compute_covariance_theta(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W)

se_theta <- sqrt(diag(covariance_theta))
theta <- c(beta, alpha, theta_nonlinear)

















