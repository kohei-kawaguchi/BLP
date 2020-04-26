rm(list = ls())
set.seed(1)

#--------#
# prepare
#--------#
library(foreach)
library(magrittr)
library(codetools)
library(doParallel)

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
# number of product characteristics excluding price
K <- 5

#---------------------#
# set hyper parameters
#---------------------#
# probability that a product is available
prob_availability <- 0.7
# partial correlation between price and xi
cor_xi_price <- 0.6
# partial correlation between price and marginal cost shocks
cor_mc_price <- 0.5

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
# taste shocks on price
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
# price
price <- 
  foreach (t = 1:T) %do% {
    price_t <- matrix(exp(rnorm(length(availability[[t]])) + cor_xi_price * xi[[t]] + cor_mc_price * mc[[t]]))

    return(price_t)
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

#-------------------------------#
# simulate endogenous vavriables
#-------------------------------#
# compute mean utility
mean_utility <- compute_mean_utility(beta, alpha, X, price, xi)

# compute individual utility
individual_utility <-
  compute_individual_utility(beta, alpha, sigma_nu, sigma_upsilon, X, price, xi, nu, upsilon)

# compute individual share
individual_share <-
  compute_individual_share(beta, alpha, sigma_nu, sigma_upsilon, X, price, xi, nu, upsilon)

# compute share
share <- compute_share(beta, alpha, sigma_nu, sigma_upsilon, X, price, xi, nu, upsilon)

# compute invidual utility from delta
individual_utility_delta <- 
  compute_individual_utility_delta(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)
max(abs(unlist(individual_utility) - unlist(individual_utility_delta)))

# compute individual share from delta
individual_share_delta <- 
  compute_individual_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)
max(abs(unlist(individual_share) - unlist(individual_share_delta)))

# compute share
share_delta <- 
  compute_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)
max(abs(unlist(share_delta) - unlist(share)))

# invert share
mean_utility <- invert_share(share, mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)

# make initial weight
X_vec <- 
  X %>%
  purrr::reduce(rbind)
Z_vec <-
  Z %>%
  purrr::reduce(rbind)
XZ <- cbind(X_vec, Z_vec)
W <- crossprod(XZ, XZ)

# estimate the linear parameters
theta_linear_hat <- estimate_linear_parameters(mean_utility, X, price, Z, W)
max(abs(theta_linear_hat - theta_linear))

# elicit xi
xi_hat <- elicit_xi(theta_linear_hat, mean_utility, X, price)

# compute moments
moments <- compute_moments(theta_nonlinear, share, mean_utility, X, price, Z, nu, upsilon, W)

# compute objective function
objective <- compute_objective(theta_nonlinear, share, mean_utility, X, price, Z, nu, upsilon, W)

# estiamte parameters
solution <- estimate_parameters(theta_nonlinear, share, mean_utility, X, price, Z, nu, upsilon, W) 
theta_nonlinear_hat <- solution$par
max(abs(theta_nonlinear_hat - theta_nonlinear))

# initial mean_utility
initial_mean_utility <-
  share %>%
  purrr::map(., ~ log(. / (1 - sum(.))))
solution <- estimate_parameters(theta_nonlinear, share, initial_mean_utility, X, price, Z, nu, upsilon, W) 
theta_nonlinear_hat <- solution$par
max(abs(theta_nonlinear_hat - theta_nonlinear))






