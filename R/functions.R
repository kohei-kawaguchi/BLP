# compute invidual utility
compute_individual_utility <-
  function(beta, alpha, sigma_nu, sigma_upsilon,
           X, price, nu, upsilon) {
    indirect_utility <-
      foreach (t = 1:T) %do% {
        X_t <- X[[t]]
        price_t <- price[[t]]
        nu_t <- nu[[t]]
        upsilon_t <- upsilon[[t]]
        mean_indirect_utility_t <-
          X_t %*% beta + price_t * alpha
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