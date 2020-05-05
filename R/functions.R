# compute mean utility
compute_mean_utility <-
  function(beta, alpha, X, p, xi) {
    mean_indirect_utility <-
      foreach (t = 1:length(X)) %do% {
        X_t <- X[[t]]
        p_t <- p[[t]]
        xi_t <- xi[[t]]
        mean_indirect_utility_t <-
          X_t %*% beta + p_t * alpha + xi_t
        return(mean_indirect_utility_t)
      }
    return(mean_indirect_utility)
  }

# compute invidual utility
compute_individual_utility <-
  function(beta, alpha, sigma_nu, sigma_upsilon,
           X, p, xi, nu, upsilon) {
    # copute mean utility
    mean_utility <- compute_mean_utility(beta, alpha, X, p, xi)
    # add individual tastes
    utility <-
      foreach (t = 1:length(X)) %do% {
        X_t <- X[[t]]
        p_t <- p[[t]]
        nu_t <- nu[[t]]
        upsilon_t <- upsilon[[t]]
        mean_utility_t <- mean_utility[[t]]
        shock_utility_t <- 
          X_t %*% diag(sigma_nu) %*% nu_t +
          p_t %*% sigma_upsilon %*% upsilon_t
        utility_t <-
          mean_utility_t %*% t(rep(1, ncol(nu_t))) +
          shock_utility_t
        return(utility_t)
      }
    return(utility)
  }

# compute individual share
compute_individual_share <-
  function(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon) {
    # compute individual utility
    individual_utility <-
      compute_individual_utility(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
    # compute individual shareindividual_utility
    individual_share <-
      foreach (t = 1:length(individual_utility)) %do% {
        individual_utility_t <- individual_utility[[t]]
        individual_share_t <- exp(individual_utility_t) 
        denominator_t <- 1 + matrix(apply(individual_share_t, 2, sum), nrow = 1)
        denominator_t <- matrix(rep(1, nrow(individual_share_t))) %*% denominator_t
        individual_share_t <- individual_share_t / denominator_t
        return(individual_share_t)
      }
    return(individual_share)
  }

# compute share
compute_share <-
  function(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon) {
    # compute individual share
    individual_share <-
      compute_individual_share(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
    # compute share
    share <-
      individual_share %>%
      purrr::map(., ~ apply(., 1, mean)) %>%
      purrr::map(., matrix)
    return(share)
  }

# compute invidual utility from delta
compute_individual_utility_delta <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
    indirect_utility <-
      foreach (t = 1:length(mean_utility)) %do% {
        X_t <- X[[t]]
        p_t <- p[[t]]
        nu_t <- nu[[t]]
        upsilon_t <- upsilon[[t]]
        mean_indirect_utility_t <- mean_utility[[t]]
        shock_indirect_utility_t <- 
          X_t %*% diag(sigma_nu) %*% nu_t +
          p_t %*% sigma_upsilon %*% upsilon_t
        indirect_utility_t <-
          mean_indirect_utility_t %*% t(rep(1, ncol(nu_t))) +
          shock_indirect_utility_t
        return(indirect_utility_t)
      }
    return(indirect_utility)
  }

# compute individual share from delta
compute_individual_share_delta <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
    # compute individual utility
    individual_utility <-
      compute_individual_utility_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute individual share
    individual_share <-
      foreach (t = 1:length(individual_utility)) %do% {
        individual_utility_t <- individual_utility[[t]]
        individual_share_t <- exp(individual_utility_t) 
        denominator_t <- 1 + matrix(apply(individual_share_t, 2, sum), nrow = 1)
        denominator_t <- matrix(rep(1, nrow(individual_share_t))) %*% denominator_t
        individual_share_t <- individual_share_t / denominator_t
        return(individual_share_t)
      }
    return(individual_share)
  }

# compute share
compute_share_delta <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
    # compute individual share
    individual_share <-
      compute_individual_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute share
    share <-
      individual_share %>%
      purrr::map(., ~ apply(., 1, mean)) %>%
      purrr::map(., matrix)
    return(share)
  }

# invert share
invert_share <-
  function(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
    # distance 
    distance <- 100
    while (distance > 1e-12) {
      share_delta <- 
        compute_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
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
  function(mean_utility, X, p, Z, W) {
    mean_utility_vec <- 
      mean_utility %>%
      purrr::reduce(rbind)
    X_vec <- 
      X %>%
      purrr::reduce(rbind)
    p_vec <-
      p %>%
      purrr::reduce(rbind)
    Z_vec <-
      Z %>%
      purrr::reduce(rbind)
    XP <- cbind(X_vec, p_vec)
    XZ <- cbind(X_vec, Z_vec)
    XZ <- cbind(XZ, XZ[, 2:ncol(XZ)]^2, XZ[, 2:ncol(XZ)]^3)
    A <- crossprod(XP, XZ) %*% qr.solve(W, crossprod(XZ, XP))
    b <- crossprod(XP, XZ) %*% qr.solve(W, crossprod(XZ, mean_utility_vec))
    theta_linear_hat <- qr.solve(A, b)
    return(theta_linear_hat)
  }

# elicit xi
elicit_xi <-
  function(linear_parameters, mean_utility, X, p) {
    beta_hat <- linear_parameters[1:(length(linear_parameters) - 1)]
    alpha_hat <- linear_parameters[length(linear_parameters)]
    xi_hat <-
      foreach (t = 1:length(mean_utility)) %do% {
        xi_hat_t <-
          mean_utility[[t]] - X[[t]] %*% beta_hat - p[[t]] * alpha_hat
        return(xi_hat_t)
      }
    return(xi_hat) 
  }

# compute moments
compute_moments <-
  function(theta_nonlinear, share, mean_utility, X, p, Z, nu, upsilon, W) {
    # extract
    sigma_nu <- theta_nonlinear[1:(length(theta_nonlinear) - 1)]
    sigma_upsilon <- theta_nonlinear[length(theta_nonlinear)]
    # invert share
    mean_utility <- invert_share_rcpp(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # estimate the linear parameters
    theta_linear_hat <- estimate_linear_parameters_rcpp(mean_utility, X, p, Z, W)
    # elicit xi
    xi_hat <- elicit_xi_rcpp(theta_linear_hat, mean_utility, X, p)
    # XZ
    XZ <-
      purrr::map2(X, Z, cbind)
    XZ <- 
      XZ %>%
      purrr::map(., ~ cbind(., .[, 2:ncol(.)]^2, .[, 2:ncol(.)]^3))
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
  function(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W) {
    # adjust parameters
    theta_nonlinear_full <- rep(0, length(rl))
    theta_nonlinear_full[which(rl == 1)] <- theta_nonlinear
    theta_nonlinear_full <- matrix(theta_nonlinear_full)
    # compute moments
    moments <- compute_moments(theta_nonlinear_full, share, mean_utility, X, p, Z, nu, upsilon, W)
    # compute objective 
    objective <- crossprod(moments, solve(W, moments)) / 2
    return(objective)
  }

# estiamte parameters
estimate_parameters <-
  function(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W) {
    solution <-
      optim(
        par = theta_nonlinear,
        fn = compute_objective,
        method = "Nelder-Mead",
        rl = rl,
        share = share,
        mean_utility = mean_utility,
        X = X,
        p = p,
        Z = Z,
        nu = nu,
        upsilon = upsilon,
        W = W
      )
    return(solution)  
  }

# compute derivatives of share with respect to delta
compute_share_derivatives_wrt_mean_utility <-
  function(individual_share) {
    share_derivatives_wrt_mean_utility <-
      foreach (t = 1:length(individual_share)) %do% {
        s_t <- individual_share[[t]]
        share_derivatives_wrt_mean_utility_t <-
          foreach (i = 1:ncol(s_t)) %do% {
            s_ti <- s_t[, i, drop = FALSE]
            ss_ti <- diag(as.numeric(s_ti)) - tcrossprod(s_ti, s_ti)
            return(ss_ti)
          }
        share_derivatives_wrt_mean_utility_t <-
          share_derivatives_wrt_mean_utility_t %>%
          purrr::reduce(., `+`)
        share_derivatives_wrt_mean_utility_t <- 
          share_derivatives_wrt_mean_utility_t / ncol(s_t)
        return(share_derivatives_wrt_mean_utility_t)
      }
    return(share_derivatives_wrt_mean_utility)
  }

# compute derivatives of share with respect to delta by numDeriv
compute_share_derivatives_wrt_mean_utility_numDeriv <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
    
    share_derivatives_wrt_mean_utility <-
      foreach (t = 1:length(mean_utility)) %do% {
        
        mean_utility_t <- mean_utility[[t]]
        
        fn <- 
          function(mean_utility_t) {
            # compute individual share
            individual_share <-
              compute_individual_share_delta(list(mean_utility_t), sigma_nu, sigma_upsilon, X[t], p[t], nu[t], upsilon[t])
            # compute share
            share <-
              individual_share %>%
              purrr::map(., ~ apply(., 1, mean)) %>%
              purrr::map(., matrix)
            share <- share[[1]]
            return(share)
          }
        
        share_derivatives_wrt_mean_utility_t <-
          numDeriv::jacobian(func = fn, x = mean_utility_t)
        
        return(share_derivatives_wrt_mean_utility_t)
      }
    return(share_derivatives_wrt_mean_utility)
  }


# compute derivatives of share with respect to non-linear parameters
compute_share_derivatives_wrt_theta_nonlinear <-
  function(individual_share, X, p, nu, upsilon) {
    share_derivatives_wrt_theta_nonlinear <-
      foreach (t = 1:length(individual_share)) %do% {
        s_t <- individual_share[[t]]
        X_t <- X[[t]]
        p_t <- p[[t]]
        nu_t <- nu[[t]]
        upsilon_t <- upsilon[[t]]
        Xp_t <- cbind(X_t, p_t)
        taste_t <- rbind(nu_t, upsilon_t)
        share_derivatives_wrt_theta_nonlinear_t <-
          foreach (i = 1:ncol(s_t)) %do% {
            s_ti <- s_t[, i, drop = FALSE]
            taste_ti <- taste_t[, i, drop = FALSE]
            xs <- crossprod(s_ti, Xp_t)
            xs <- Xp_t - matrix(rep(1, nrow(Xp_t))) %*% xs
            xs <- matrix(rep(1, nrow(Xp_t))) %*% t(taste_ti) * xs
            xs <- s_ti %*% matrix(rep(1, ncol(Xp_t)), nrow = 1) * xs
            return(xs)
          }
        share_derivatives_wrt_theta_nonlinear_t <-
          share_derivatives_wrt_theta_nonlinear_t %>%
          purrr::reduce(., `+`)
        share_derivatives_wrt_theta_nonlinear_t <-
          share_derivatives_wrt_theta_nonlinear_t / ncol(s_t)
        return(share_derivatives_wrt_theta_nonlinear_t)
      }
    return(share_derivatives_wrt_theta_nonlinear)
  }



# compute derivatives of share with respect to non-linear parameters by numDeriv
compute_share_derivatives_wrt_theta_nonlinear_numDeriv <-
  function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
    
    theta_nonlinear <- c(sigma_nu, sigma_upsilon)
    
    share_derivatives_wrt_mean_utility <-
      foreach (t = 1:length(mean_utility)) %do% {
        
        fn <- 
          function(theta_nonlinear) {
            sigma_nu <- theta_nonlinear[1:(length(theta_nonlinear) - 1)]
            sigma_upsilon <- theta_nonlinear[length(theta_nonlinear)]
            # compute individual share
            individual_share <-
              compute_individual_share_delta(mean_utility[t], sigma_nu, sigma_upsilon, X[t], p[t], nu[t], upsilon[t])
            # compute share
            share <-
              individual_share %>%
              purrr::map(., ~ apply(., 1, mean)) %>%
              purrr::map(., matrix)
            share <- share[[1]]
            return(share)
          }
        
        share_derivatives_wrt_mean_utility_t <-
          numDeriv::jacobian(func = fn, x = theta_nonlinear)
        
        return(share_derivatives_wrt_mean_utility_t)
      }
    return(share_derivatives_wrt_mean_utility)
  }

# compute derivatives of mean utility with respect to non-inear parameters
compute_mean_utility_derivatives_wrt_theta_nonlinear <-
  function(individual_share, X, p, nu, upsilon) {
    # compute derivatives of share with respect to delta
    share_derivatives_wrt_mean_utility <- compute_share_derivatives_wrt_mean_utility(individual_share)
    # compute derivatives of share with respect to non-linear parameters
    share_derivatives_wrt_theta_nonlinear <- compute_share_derivatives_wrt_theta_nonlinear(individual_share, X, p, nu, upsilon)
    # derivatives
    mean_utility_derivatives_wrt_theta_nonlinear <-
      purrr::map2(share_derivatives_wrt_mean_utility, share_derivatives_wrt_theta_nonlinear, ~ - qr.solve(.x, .y))
    return(mean_utility_derivatives_wrt_theta_nonlinear)
  }

# compute derivatives of the objective function with respect to non-linear parameters
compute_objective_derivatives_wrt_theta_nonlinear <-
  function(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W) {
    # adjust parameters
    theta_nonlinear_full <- rep(0, length(rl))
    theta_nonlinear_full[which(rl == 1)] <- theta_nonlinear
    theta_nonlinear_full <- matrix(theta_nonlinear_full)
    # compute moments
    # extract
    sigma_nu <- theta_nonlinear_full[1:(length(theta_nonlinear_full) - 1)]
    sigma_upsilon <- theta_nonlinear_full[length(theta_nonlinear_full)]
    # invert share
    mean_utility <- invert_share_rcpp(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute individual share from delta
    individual_share_delta <- 
      compute_individual_share_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute derivatives of mean utility with respect to non-linear parameters
    mean_utility_derivatives_wrt_theta_nonlinear <- compute_mean_utility_derivatives_wrt_theta_nonlinear_rcpp(individual_share, X, p, nu, upsilon)
    # estimate the linear parameters
    theta_linear_hat <- estimate_linear_parameters_rcpp(mean_utility, X, p, Z, W)
    # elicit xi
    xi_hat <- elicit_xi_rcpp(theta_linear_hat, mean_utility, X, p)
    # XZ
    XZ <-
      purrr::map2(X, Z, cbind)
    XZ <- 
      XZ %>%
      purrr::map(., ~ cbind(., .[, 2:ncol(.)]^2, .[, 2:ncol(.)]^3))
    # compute moments
    moments <-
      purrr::map2(xi_hat, XZ, crossprod) %>%
      purrr::reduce(rbind)
    moments <-
      moments %>%
      apply(., 2, mean) %>%
      matrix()
    # compute derivatives of moments
    moments_derivatives <-
      purrr::map2(mean_utility_derivatives_wrt_theta_nonlinear, XZ, crossprod) %>%
      purrr::reduce(., `+`)
    moments_derivatives <-
      moments_derivatives / length(share)
    # compute objective 
    objective_derivatives <- moments_derivatives %*% solve(W, moments)
    # return
    return(objective_derivatives)
  }

# compute derivatives of the objective function with respect to non-linear parameters using numDeriv
compute_objective_derivatives_wrt_theta_nonlinear_numDeriv <-
  function(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W) {
    # adjust parameters
    theta_nonlinear_full <- rep(0, length(rl))
    theta_nonlinear_full[which(rl == 1)] <- theta_nonlinear
    theta_nonlinear_full <- matrix(theta_nonlinear_full)
    
    fn <- function(theta_nonlinear_full) {
      # compute moments
      moments <- compute_moments(theta_nonlinear_full, share, mean_utility, X, p, Z, nu, upsilon, W)
      # compute objective 
      objective <- crossprod(moments, solve(W, moments)) / 2
      return(objective)
    }
    
    objective_derivatives <-
      numDeriv::grad(func = fn, x = theta_nonlinear_full)
    
    return(objective_derivatives)
  }

# compute efficient weighting matrix
compute_efficient_weighting_matrix <-
  function(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W) {
    # adjust parameters
    theta_nonlinear_full <- rep(0, length(rl))
    theta_nonlinear_full[which(rl == 1)] <- theta_nonlinear
    theta_nonlinear_full <- matrix(theta_nonlinear_full)
    # compute moments
    # extract
    sigma_nu <- theta_nonlinear_full[1:(length(theta_nonlinear_full) - 1)]
    sigma_upsilon <- theta_nonlinear_full[length(theta_nonlinear_full)]
    # invert share
    mean_utility <- invert_share_rcpp(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute individual share from delta
    individual_share_delta <- 
      compute_individual_share_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute derivatives of mean utility with respect to non-linear parameters
    mean_utility_derivatives_wrt_theta_nonlinear <- compute_mean_utility_derivatives_wrt_theta_nonlinear_rcpp(individual_share, X, p, nu, upsilon)
    # estimate the linear parameters
    theta_linear_hat <- estimate_linear_parameters_rcpp(mean_utility, X, p, Z, W)
    # elicit xi
    xi_hat <- elicit_xi_rcpp(theta_linear_hat, mean_utility, X, p)
    
    # XZ
    XZ <-
      purrr::map2(X, Z, cbind) %>%
      purrr::reduce(., rbind)
    XZ <- cbind(XZ, XZ[, 2:ncol(XZ)]^2, XZ[, 2:ncol(XZ)]^3)
    # compute variance of moments
    moments <-
      xi_hat %>%
      purrr::reduce(rbind)
    moments <- moments %*% matrix(rep(1, ncol(XZ)), nrow = 1)
    moments <- moments *XZ
    Omega <- crossprod(moments, moments) / nrow(XZ)
    
    # return
    return(Omega)
  }

# compute standard errors
compute_covariance_theta <-
  function(theta_nonlinear, rl, share, mean_utility, X, p, Z, nu, upsilon, W) {
    # adjust parameters
    theta_nonlinear_full <- rep(0, length(rl))
    theta_nonlinear_full[which(rl == 1)] <- theta_nonlinear
    theta_nonlinear_full <- matrix(theta_nonlinear_full)
    # compute moments
    # extract
    sigma_nu <- theta_nonlinear_full[1:(length(theta_nonlinear_full) - 1)]
    sigma_upsilon <- theta_nonlinear_full[length(theta_nonlinear_full)]
    # invert share
    mean_utility <- invert_share_rcpp(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute individual share from delta
    individual_share_delta <- 
      compute_individual_share_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
    # compute derivatives of mean utility with respect to non-linear parameters
    mean_utility_derivatives_wrt_theta_nonlinear <- compute_mean_utility_derivatives_wrt_theta_nonlinear_rcpp(individual_share, X, p, nu, upsilon)
    # estimate the linear parameters
    theta_linear_hat <- estimate_linear_parameters_rcpp(mean_utility, X, p, Z, W)
    # elicit xi
    xi_hat <- elicit_xi_rcpp(theta_linear_hat, mean_utility, X, p)
    
    # XZ
    XZ <-
      purrr::map2(X, Z, cbind) %>%
      purrr::reduce(., rbind)
    XZ <-
      cbind(XZ, XZ[, 2:ncol(XZ)]^2, XZ[, 2:ncol(XZ)]^3)
    Xp <-
      purrr::map2(X, p, cbind) %>%
      purrr::reduce(., rbind)
    # compute variance of moments
    moments <-
      xi_hat %>%
      purrr::reduce(rbind)
    moments <- moments %*% matrix(rep(1, ncol(XZ)), nrow = 1)
    moments <- moments * XZ
    Omega <- crossprod(moments, moments) / nrow(XZ)
    
    # compute derivatives of moments wrt non-linear parameters
    G <-purrr::reduce(mean_utility_derivatives_wrt_theta_nonlinear, rbind)
    # compute derivatives of moments wrt linear parameters
    G <- cbind(G, Xp)
    
    G <- crossprod(XZ, G) / nrow(XZ)
    # compute asymptotic covariance of the estimator
    L <- crossprod(G, qr.solve(W, G))
    C <- crossprod(G, qr.solve(W, Omega)) %*% qr.solve(W, G)
    covariance <- qr.solve(L, C) %*% qr.solve(L)
    covariance <- covariance / nrow(XZ)
    
    # return
    return(covariance)
  }
