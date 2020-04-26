# compute mean utility
compute_mean_utility <-
  function(beta, alpha, X, price, xi) {
    mean_indirect_utility <-
      foreach (t = 1:length(X)) %do% {
        X_t <- X[[t]]
        price_t <- price[[t]]
        mean_indirect_utility_t <-
          X_t %*% beta + price_t * alpha
        return(mean_indirect_utility_t)
      }
    return(mean_indirect_utility)
  }

# compute invidual utility
compute_individual_utility <-
  function(beta, alpha, sigma_nu, sigma_upsilon,
           X, price, xi, nu, upsilon) {
    # copute mean utility
    mean_utility <- compute_mean_utility(beta, alpha, X, price, xi)
    # add individual tastes
    utility <-
      foreach (t = 1:length(X)) %do% {
        X_t <- X[[t]]
        price_t <- price[[t]]
        nu_t <- nu[[t]]
        xi_t <- xi[[t]]
        upsilon_t <- upsilon[[t]]
        mean_utility_t <- mean_utility[[t]]
        shock_utility_t <- 
          X_t %*% diag(sigma_nu) %*% nu_t +
          price_t %*% sigma_upsilon %*% upsilon_t
        utility_t <-
          mean_utility_t %*% t(rep(1, ncol(nu_t))) +
          shock_utility_t
        return(utility_t)
      }
    return(utility)
  }

# compute individual share
compute_individual_share <-
  function(beta, alpha, sigma_nu, sigma_upsilon, X, price, xi, nu, upsilon) {
    # compute individual utility
    individual_utility <-
      compute_individual_utility(beta, alpha, sigma_nu, sigma_upsilon, X, price, xi, nu, upsilon)
    # compute individual shareindividual_utility
    individual_share <-
      foreach (t = 1:length(individual_utility)) %do% {
        individual_utility_t <- individual_utility[[t]]
        individual_share_t <- exp(individual_utility_t) 
        denominator_t <- matrix(apply(1 + individual_share_t, 2, sum), nrow = 1)
        denominator_t <- matrix(rep(1, nrow(individual_share_t))) %*% denominator_t
        individual_share_t <- individual_share_t / denominator_t
        return(individual_share_t)
      }
    return(individual_share)
  }

# compute share
compute_share <-
  function(beta, alpha, sigma_nu, sigma_upsilon, X, price, xi, nu, upsilon) {
    # compute individual share
    individual_share <-
      compute_individual_share(beta, alpha, sigma_nu, sigma_upsilon, X, price, xi, nu, upsilon)
    # compute share
    share <-
      individual_share %>%
      purrr::map(., ~ apply(., 1, mean)) %>%
      purrr::map(., matrix)
    return(share)
  }

# compute invidual utility from delta
compute_individual_utility_delta <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon) {
    indirect_utility <-
      foreach (t = 1:length(mean_utility)) %do% {
        X_t <- X[[t]]
        price_t <- price[[t]]
        nu_t <- nu[[t]]
        upsilon_t <- upsilon[[t]]
        mean_indirect_utility_t <- mean_utility[[t]]
        shock_indirect_utility_t <- 
          X_t %*% diag(sigma_nu) %*% nu_t +
          price_t %*% sigma_upsilon %*% upsilon_t
        indirect_utility_t <-
          mean_indirect_utility_t %*% t(rep(1, ncol(nu_t))) +
          shock_indirect_utility_t
        return(indirect_utility_t)
      }
    return(indirect_utility)
  }

# compute individual share from delta
compute_individual_share_delta <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon) {
    # compute individual utility
    individual_utility <-
      compute_individual_utility_delta(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)
    # compute individual shareindividual_utility
    individual_share <-
      foreach (t = 1:length(individual_utility)) %do% {
        individual_utility_t <- individual_utility[[t]]
        individual_share_t <- exp(individual_utility_t) 
        denominator_t <- matrix(apply(1 + individual_share_t, 2, sum), nrow = 1)
        denominator_t <- matrix(rep(1, nrow(individual_share_t))) %*% denominator_t
        individual_share_t <- individual_share_t / denominator_t
        return(individual_share_t)
      }
    return(individual_share)
  }

# compute share
compute_share_delta <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon) {
    # compute individual share
    individual_share <-
      compute_individual_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)
    # compute share
    share <-
      individual_share %>%
      purrr::map(., ~ apply(., 1, mean)) %>%
      purrr::map(., matrix)
    return(share)
  }

# invert share
invert_share <-
  function(share, mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon) {
    # distance 
    distance <- 100
    while (distance > 1e-8) {
      share_delta <- 
        compute_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)
      update <-
        purrr::map2(share, share_delta, ~ log(.x) - log(.y))
      mean_utility <- 
        purrr::map2(mean_utility, update, ~ .x + .y)
      distance <- max(abs(unlist(update)))
    }
    return(mean_utility)
  }

# estimate the linear parameters
estimate_linear_parameters <-
  function(mean_utility, X, price, Z, W) {
    mean_utility_vec <- 
      mean_utility %>%
      purrr::reduce(rbind)
    X_vec <- 
      X %>%
      purrr::reduce(rbind)
    price_vec <-
      price %>%
      purrr::reduce(rbind)
    Z_vec <-
      Z %>%
      purrr::reduce(rbind)
    XP <- cbind(X_vec, price_vec)
    XZ <- cbind(X_vec, Z_vec)
    W <- crossprod(XZ, XZ)
    A <- crossprod(XP, XZ) %*% solve(W, crossprod(XZ, XP))
    b <- crossprod(XP, XZ) %*% solve(W, crossprod(XZ, mean_utility_vec))
    theta_linear_hat <- solve(A, b)
    return(theta_linear_hat)
  }

# elicit xi
elicit_xi <-
  function(linear_parameters, mean_utility, X, price) {
    beta_hat <- linear_parameters[1:(length(linear_parameters) - 1)]
    alpha_hat <- linear_parameters[length(linear_parameters)]
    xi_hat <-
      foreach (t = 1:length(mean_utility)) %do% {
        xi_hat_t <-
          mean_utility[[t]] - X[[t]] %*% beta_hat - price[[t]] * alpha_hat
        return(xi_hat_t)
      }
    return(xi_hat) 
  }

# compute moments
compute_moments <-
  function(theta_nonlinear, share, mean_utility, X, price, Z, nu, upsilon, W) {
    # extract
    sigma_nu <- theta_nonlinear[1:(length(theta_nonlinear) - 1)]
    sigma_upsilon <- theta_nonlinear[length(theta_nonlinear)]
    # invert share
    mean_utility <- invert_share(share, mean_utility, sigma_nu, sigma_upsilon, X, price, nu, upsilon)
    # estimate the linear parameters
    theta_linear_hat <- estimate_linear_parameters(mean_utility, X, price, Z, W)
    # elicit xi
    xi_hat <- elicit_xi(theta_linear_hat, mean_utility, X, price)
    # XZ
    XZ <-
      purrr::map2(X, Z, cbind)
    # compute moments
    moments <-
      purrr::map2(xi_hat, XZ, crossprod) %>%
      purrr::reduce(rbind)
    moments <-
      moments %>%
      apply(., 2, mean) %>%
      matrix()
    # return
    return(moments)
  }

# compute objective function
compute_objective <-
  function(theta_nonlinear, share, mean_utility, X, price, Z, nu, upsilon, W) {
    # compute moments
    moments <- compute_moments(theta_nonlinear, share, mean_utility, X, price, Z, nu, upsilon, W)
    # compute objective 
    objective <- crossprod(moments, solve(W, moments)) / 2
    return(objective)
  }

# estiamte parameters
estimate_parameters <-
  function(theta_nonlinear, share, mean_utility, X, price, Z, nu, upsilon, W) {
    solution <-
      optim(
        par = theta_nonlinear,
        fn = compute_objective,
        method = "Nelder-Mead",
        share = share,
        mean_utility = mean_utility,
        X = X,
        price = price,
        Z = Z,
        nu = nu,
        upsilon = upsilon,
        W = W
      )
    return(solution)  
  }
