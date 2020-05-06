#include <Rcpp.h>
#include <RcppEigen.h>

// # compute mean utility
// compute_mean_utility <-
// function(beta, alpha, X, p, xi) {
//   mean_indirect_utility <-
//     foreach (t = 1:length(X)) %do% {
//       X_t <- X[[t]]
//       p_t <- p[[t]]
//       mean_indirect_utility_t <-
//         X_t %*% beta + p_t * alpha
//       return(mean_indirect_utility_t)
//     }
//     return(mean_indirect_utility)
// }
// [[Rcpp::export]]
Rcpp::List compute_mean_utility_rcpp(
    Eigen::MatrixXd beta, 
    double alpha, 
    Rcpp::List X, 
    Rcpp::List p,
    Rcpp::List xi
) {
  Rcpp::List mean_indirect_utility;
  for (int t = 0; t < X.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> X_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (X.at(t)));
    Eigen::Map<Eigen::MatrixXd> p_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (p.at(t)));
    Eigen::Map<Eigen::MatrixXd> xi_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (xi.at(t)));
    Eigen::MatrixXd mean_indirect_utility_t = X_t * beta + p_t * alpha + xi_t;
    mean_indirect_utility.push_back(mean_indirect_utility_t);
  }
  return(mean_indirect_utility);
}

// # compute invidual utility
// compute_individual_utility <-
//   function(beta, alpha, sigma_nu, sigma_upsilon,
//            X, p, xi, nu, upsilon) {
// # copute mean utility
//     mean_utility <- compute_mean_utility(beta, alpha, X, p, xi)
// # add individual tastes
//     utility <-
//       foreach (t = 1:length(X)) %do% {
//         X_t <- X[[t]]
//         p_t <- p[[t]]
//         nu_t <- nu[[t]]
//         upsilon_t <- upsilon[[t]]
//         mean_utility_t <- mean_utility[[t]]
//         shock_utility_t <- 
//           X_t %*% diag(sigma_nu) %*% nu_t +
//           p_t %*% sigma_upsilon %*% upsilon_t
//         utility_t <-
//           mean_utility_t %*% t(rep(1, ncol(nu_t))) +
//           shock_utility_t
//         return(utility_t)
//       }
//       return(utility)
//   }
// [[Rcpp::export]]
Rcpp::List compute_individual_utility_rcpp(
    Eigen::MatrixXd beta, 
    double alpha, 
    Eigen::VectorXd sigma_nu, 
    double sigma_upsilon, 
    Rcpp::List X, 
    Rcpp::List p, 
    Rcpp::List xi, 
    Rcpp::List nu, 
    Rcpp::List upsilon
) {
  // compute mean utility
  Rcpp::List mean_utility = compute_mean_utility_rcpp(beta, alpha, X, p, xi);
  // add individual tastes
  Rcpp::List utility;
  for (int t = 0; t < X.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> X_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (X.at(t)));
    Eigen::Map<Eigen::MatrixXd> p_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (p.at(t)));
    Eigen::Map<Eigen::MatrixXd> nu_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (nu.at(t)));
    Eigen::Map<Eigen::MatrixXd> upsilon_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (upsilon.at(t)));
    Eigen::Map<Eigen::MatrixXd> mean_utility_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (mean_utility.at(t)));
    Eigen::MatrixXd shock_utility_t = X_t * sigma_nu.asDiagonal() * nu_t +
      p_t * sigma_upsilon * upsilon_t;
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(1, nu_t.cols());
    Eigen::MatrixXd utility_t = mean_utility_t * ones + shock_utility_t;
    utility.push_back(utility_t);
  }
  return(utility);
}

// # compute individual share
// compute_individual_share <-
//   function(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon) {
// # compute individual utility
//     individual_utility <-
//       compute_individual_utility(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
// # compute individual shareindividual_utility
//     individual_share <-
//       foreach (t = 1:length(individual_utility)) %do% {
//         individual_utility_t <- individual_utility[[t]]
//         individual_share_t <- exp(individual_utility_t) 
//         denominator_t <- 1 + matrix(apply(individual_share_t, 2, sum), nrow = 1)
//         denominator_t <- matrix(rep(1, nrow(individual_share_t))) %*% denominator_t
//         individual_share_t <- individual_share_t / denominator_t
//         return(individual_share_t)
//       }
//       return(individual_share)
//   }
// [[Rcpp::export]]
Rcpp::List compute_individual_share_rcpp(
    Eigen::MatrixXd beta,
    double alpha,
    Eigen::VectorXd sigma_nu,
    double sigma_upsilon,
    Rcpp::List X,
    Rcpp::List p,
    Rcpp::List xi,
    Rcpp::List nu,
    Rcpp::List upsilon
) {
  // compute individual utility
  Rcpp::List individual_utility =
    compute_individual_utility_rcpp(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon);
  // compute individual shareindividual_utility
  Rcpp::List individual_share;
  for (int t = 0; t < X.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> individual_utility_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (individual_utility.at(t)));
    Eigen::MatrixXd individual_share_t = individual_utility_t.array().exp();
    Eigen::MatrixXd denominator_t = individual_share_t.colwise().sum();
    denominator_t = Eigen::MatrixXd::Ones(1, denominator_t.cols()) + denominator_t;
    denominator_t = Eigen::MatrixXd::Ones(individual_share_t.rows(), 1) * denominator_t;
    individual_share_t = (individual_share_t.array() / denominator_t.array()).matrix();
    individual_share.push_back(individual_share_t);
  }
  return(individual_share);
}

// [[Rcpp::export]]
Eigen::MatrixXd test(Eigen::MatrixXd X) {
  Eigen::MatrixXd Y = X.colwise().sum();
  return(Y);
}

// # compute share
// compute_share <-
//   function(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon) {
// # compute individual share
//     individual_share <-
//       compute_individual_share(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon)
// # compute share
//     share <-
//       individual_share %>%
//       purrr::map(., ~ apply(., 1, mean)) %>%
//       purrr::map(., matrix)
//     return(share)
//   }
// [[Rcpp::export]]
Rcpp::List compute_share_rcpp(
    Eigen::MatrixXd beta,
    double alpha,
    Eigen::VectorXd sigma_nu,
    double sigma_upsilon,
    Rcpp::List X,
    Rcpp::List p,
    Rcpp::List xi,
    Rcpp::List nu,
    Rcpp::List upsilon
) {
  // compute individual share
  Rcpp::List individual_share =
    compute_individual_share_rcpp(beta, alpha, sigma_nu, sigma_upsilon, X, p, xi, nu, upsilon);
  // compute share
  Rcpp::List share;
  for (int t = 0; t < individual_share.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> individual_share_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (individual_share.at(t)));
    Eigen::MatrixXd share_t = individual_share_t.rowwise().mean();
    share.push_back(share_t);
  }
  return(share);
}

// # compute invidual utility from delta
// compute_individual_utility_delta <-
//   function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
//     indirect_utility <-
//       foreach (t = 1:length(mean_utility)) %do% {
//         X_t <- X[[t]]
//         p_t <- p[[t]]
//         nu_t <- nu[[t]]
//         upsilon_t <- upsilon[[t]]
//         mean_indirect_utility_t <- mean_utility[[t]]
//         shock_indirect_utility_t <- 
//           X_t %*% diag(sigma_nu) %*% nu_t +
//           p_t %*% sigma_upsilon %*% upsilon_t
//         indirect_utility_t <-
//           mean_indirect_utility_t %*% t(rep(1, ncol(nu_t))) +
//           shock_indirect_utility_t
//         return(indirect_utility_t)
//       }
//       return(indirect_utility)
//   }
// [[Rcpp::export]]
Rcpp::List compute_individual_utility_delta_rcpp(
    Rcpp::List mean_utility, 
    Eigen::VectorXd sigma_nu, 
    double sigma_upsilon, 
    Rcpp::List X, 
    Rcpp::List p, 
    Rcpp::List nu, 
    Rcpp::List upsilon
) {
  Rcpp::List indirect_utility;
  for (int t = 0; t < mean_utility.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> X_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (X.at(t)));
    Eigen::Map<Eigen::MatrixXd> p_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (p.at(t)));
    Eigen::Map<Eigen::MatrixXd> nu_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (nu.at(t)));
    Eigen::Map<Eigen::MatrixXd> upsilon_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (upsilon.at(t)));
    Eigen::Map<Eigen::MatrixXd> mean_indirect_utility_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (mean_utility.at(t)));
    Eigen::MatrixXd shock_indirect_utility_t =
      X_t * sigma_nu.asDiagonal() * nu_t +
      p_t * sigma_upsilon * upsilon_t;
    Eigen::MatrixXd indirect_utility_t = mean_indirect_utility_t * Eigen::MatrixXd::Ones(1, nu_t.cols()) +
      shock_indirect_utility_t;
    indirect_utility.push_back(indirect_utility_t);
  }
  return(indirect_utility);
}

// # compute individual share from delta
// compute_individual_share_delta <-
//   function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
// # compute individual utility
//     individual_utility <-
//       compute_individual_utility_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
// # compute individual shareindividual_utility
//     individual_share <-
//       foreach (t = 1:length(individual_utility)) %do% {
//         individual_utility_t <- individual_utility[[t]]
//         individual_share_t <- exp(individual_utility_t) 
//         denominator_t <- 1 + matrix(apply(individual_share_t, 2, sum), nrow = 1)
//         denominator_t <- matrix(rep(1, nrow(individual_share_t))) %*% denominator_t
//         individual_share_t <- individual_share_t / denominator_t
//         return(individual_share_t)
//       }
//       return(individual_share)
//   }
// [[Rcpp::export]]
Rcpp::List compute_individual_share_delta_rcpp(
    Rcpp::List mean_utility, 
    Eigen::VectorXd sigma_nu, 
    double sigma_upsilon, 
    Rcpp::List X, 
    Rcpp::List p, 
    Rcpp::List nu, 
    Rcpp::List upsilon  
) {
  // compute individual utility
  Rcpp::List individual_utility = 
    compute_individual_utility_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon);
  // compute individual share
  Rcpp::List individual_share;
  for (int t = 0; t < X.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> individual_utility_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (individual_utility.at(t)));
    Eigen::MatrixXd individual_share_t = individual_utility_t.array().exp();
    Eigen::MatrixXd denominator_t = individual_share_t.colwise().sum();
    denominator_t = Eigen::MatrixXd::Ones(1, denominator_t.cols()) + denominator_t;
    denominator_t = Eigen::MatrixXd::Ones(individual_share_t.rows(), 1) * denominator_t;
    individual_share_t = (individual_share_t.array() / denominator_t.array()).matrix();
    individual_share.push_back(individual_share_t);
  }
  return(individual_share);
}

// # compute share
// compute_share_delta <-
//   function(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
// # compute individual share
//     individual_share <-
//       compute_individual_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
// # compute share
//     share <-
//       individual_share %>%
//       purrr::map(., ~ apply(., 1, mean)) %>%
//       purrr::map(., matrix)
//     return(share)
//   }
// [[Rcpp::export]]
Rcpp::List compute_share_delta_rcpp(
    Rcpp::List mean_utility, 
    Eigen::VectorXd sigma_nu, 
    double sigma_upsilon, 
    Rcpp::List X, 
    Rcpp::List p, 
    Rcpp::List nu, 
    Rcpp::List upsilon  
) {
  // compute individual share
  Rcpp::List individual_share =
    compute_individual_share_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon);
  // compute share
  Rcpp::List share;
  for (int t = 0; t < individual_share.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> individual_share_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (individual_share.at(t)));
    Eigen::MatrixXd share_t = individual_share_t.rowwise().mean();
    share.push_back(share_t);
  }
  return(share);
}

// # invert share
// invert_share <-
//   function(share, mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon) {
// # distance 
//     distance <- 100
//     while (distance > 1e-12) {
//       share_delta <- 
//         compute_share_delta(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon)
//       update <-
//         purrr::map2(share, share_delta, ~ log(.x) - log(.y))
//       mean_utility <- 
//         purrr::map2(mean_utility, update, ~ .x + .y)
//       distance <- max(abs(unlist(update)))
//     }
//     return(mean_utility)
//   }
// [[Rcpp::export]]
Rcpp::List invert_share_rcpp(
    Rcpp::List share,
    Rcpp::List mean_utility, 
    Eigen::VectorXd sigma_nu, 
    double sigma_upsilon, 
    Rcpp::List X, 
    Rcpp::List p, 
    Rcpp::List nu, 
    Rcpp::List upsilon  
) {
  double distance = 100;
  while (distance > 1e-12) {
    Rcpp::List share_delta = 
      compute_share_delta_rcpp(mean_utility, sigma_nu, sigma_upsilon, X, p, nu, upsilon);
    distance = 0;
    Rcpp::List mean_utility_new;
    for (int t = 0; t < mean_utility.size(); t++) {
      Eigen::Map<Eigen::MatrixXd> mean_utility_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (mean_utility.at(t)));
      Eigen::Map<Eigen::MatrixXd> share_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (share.at(t)));
      Eigen::Map<Eigen::MatrixXd> share_delta_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (share_delta.at(t)));
      Eigen::MatrixXd update_t = share_t.array().log() - share_delta_t.array().log();
      distance = std::max(distance, update_t.array().abs().maxCoeff());
      mean_utility_t = mean_utility_t + update_t;
      mean_utility.at(t) = mean_utility_t;
    }
  }
  return(mean_utility);
}


// [[Rcpp::export]]
Eigen::MatrixXd vstack_rcpp(
    Rcpp::List L
) {
  int row_size = 0;
  int col_size = 0;
  for (int t = 0; t < L.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> L_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (L.at(t)));
    row_size = row_size + L_t.rows();
    col_size = L_t.cols();
  }
  Eigen::MatrixXd L_vec = Eigen::MatrixXd::Zero(row_size, col_size);
  row_size = 0;
  for (int t = 0; t < L.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> L_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (L.at(t)));
    L_vec.block(row_size, 0, L_t.rows(), L_t.cols()) = L_t;
    row_size = row_size + L_t.rows();
  }
  return(L_vec);
}

// # estimate the linear parameters
// estimate_linear_parameters <-
//   function(mean_utility, X, p, Z, W) {
//     mean_utility_vec <- 
//       mean_utility %>%
//       purrr::reduce(rbind)
//     X_vec <- 
//       X %>%
//       purrr::reduce(rbind)
//     p_vec <-
//       p %>%
//       purrr::reduce(rbind)
//     Z_vec <-
//       Z %>%
//       purrr::reduce(rbind)
//     XP <- cbind(X_vec, p_vec)
//     XZ <- cbind(X_vec, Z_vec)
//     A <- crossprod(XP, XZ) %*% solve(W, crossprod(XZ, XP))
//     b <- crossprod(XP, XZ) %*% solve(W, crossprod(XZ, mean_utility_vec))
//     theta_linear_hat <- solve(A, b)
//     return(theta_linear_hat)
//   }
// [[Rcpp::export]]
Eigen::MatrixXd estimate_linear_parameters_rcpp(
    Rcpp::List mean_utility, 
    Rcpp::List X, 
    Rcpp::List p, 
    Rcpp::List instruments,
    Eigen::MatrixXd W
) { 
  Eigen::MatrixXd mean_utility_vec = vstack_rcpp(mean_utility);
  Eigen::MatrixXd X_vec = vstack_rcpp(X);
  Eigen::MatrixXd p_vec = vstack_rcpp(p);
  Eigen::MatrixXd XP(X_vec.rows(), X_vec.cols() + p_vec.cols());
  XP << X_vec, p_vec;
  Eigen::MatrixXd instruments_vec = vstack_rcpp(instruments);
  Eigen::MatrixXd A = (XP.transpose() * instruments_vec) * (W.colPivHouseholderQr().solve(instruments_vec.transpose() * XP));
  Eigen::MatrixXd b = (XP.transpose() * instruments_vec) * W.colPivHouseholderQr().solve(instruments_vec.transpose() * mean_utility_vec);
  Eigen::MatrixXd theta_linear_hat = A.colPivHouseholderQr().solve(b);
  return(theta_linear_hat);
}

// # elicit xi
// elicit_xi <-
//   function(linear_parameters, mean_utility, X, p) {
//     beta_hat <- linear_parameters[1:(length(linear_parameters) - 1)]
//     alpha_hat <- linear_parameters[length(linear_parameters)]
//     xi_hat <-
//       foreach (t = 1:length(mean_utility)) %do% {
//         xi_hat_t <-
//           mean_utility[[t]] - X[[t]] %*% beta_hat - p[[t]] * alpha_hat
//         return(xi_hat_t)
//       }
//       return(xi_hat) 
//   }
// [[Rcpp::export]]
Rcpp::List elicit_xi_rcpp(
    Eigen::MatrixXd linear_parameters, 
    Rcpp::List mean_utility, 
    Rcpp::List X, 
    Rcpp::List p
) {
  Eigen::MatrixXd beta_hat = linear_parameters.block(0, 0, linear_parameters.rows() - 1, 1);
  double alpha_hat = linear_parameters(linear_parameters.rows() - 1, 0);
  Rcpp::List xi_hat;
  for (int t = 0; t < mean_utility.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> mean_utility_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (mean_utility.at(t)));
    Eigen::Map<Eigen::MatrixXd> X_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (X.at(t)));
    Eigen::Map<Eigen::MatrixXd> p_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (p.at(t)));
    Eigen::MatrixXd xi_hat_t = mean_utility_t - X_t * beta_hat - p_t * alpha_hat;
    xi_hat.push_back(xi_hat_t);
  }
  return(xi_hat);
}

// # compute derivatives of share with respect to delta
// compute_share_derivatives_wrt_mean_utility <-
//   function(individual_share) {
//     share_derivatives_wrt_mean_utility <-
//       foreach (t = 1:length(individual_share)) %do% {
//         s_t <- individual_share[[t]]
//         share_derivatives_wrt_mean_utility_t <-
//           foreach (i = 1:ncol(s_t)) %do% {
//             s_ti <- s_t[, i, drop = FALSE]
//             ss_ti <- diag(as.numeric(s_ti)) - tcrossprod(s_ti, s_ti)
//             return(ss_ti)
//           }
//           share_derivatives_wrt_mean_utility_t <-
//             share_derivatives_wrt_mean_utility_t %>%
//             purrr::reduce(., `+`)
//             share_derivatives_wrt_mean_utility_t <- 
//               share_derivatives_wrt_mean_utility_t / ncol(s_t)
//             return(share_derivatives_wrt_mean_utility_t)
//       }
//       return(share_derivatives_wrt_mean_utility)
//   }
// [[Rcpp::export]]
Rcpp::List compute_share_derivatives_wrt_mean_utility_rcpp(
  Rcpp::List individual_share
) {
  Rcpp::List share_derivatives_wrt_mean_utility;
  for (int t = 0; t < individual_share.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> s_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (individual_share.at(t)));
    Eigen::MatrixXd share_derivatives_wrt_mean_utility_t = Eigen::MatrixXd::Zero(s_t.rows(), s_t.rows());
    for (int i = 0; i < s_t.cols(); i++) {
      Eigen::VectorXd s_ti = s_t.col(i);
      Eigen::MatrixXd ss_ti = Eigen::MatrixXd(s_ti.asDiagonal()) - s_ti * s_ti.transpose();
      share_derivatives_wrt_mean_utility_t = share_derivatives_wrt_mean_utility_t + ss_ti;
    }
    share_derivatives_wrt_mean_utility_t = share_derivatives_wrt_mean_utility_t / s_t.cols();
    share_derivatives_wrt_mean_utility.push_back(share_derivatives_wrt_mean_utility_t);
  }
  return(share_derivatives_wrt_mean_utility);
}

// # compute derivatives of share with respect to non-linear parameters
// compute_share_derivatives_wrt_theta_nonlinear <-
//   function(individual_share, X, p, nu, upsilon) {
//     share_derivatives_wrt_theta_nonlinear <-
//       foreach (t = 1:length(individual_share)) %do% {
//         s_t <- individual_share[[t]]
//         X_t <- X[[t]]
//         p_t <- p[[t]]
//         nu_t <- nu[[t]]
//         upsilon_t <- upsilon[[t]]
//         Xp_t <- cbind(X_t, p_t)
//         taste_t <- rbind(nu_t, upsilon_t)
//         share_derivatives_wrt_theta_nonlinear_t <-
//           foreach (i = 1:ncol(s_t)) %do% {
//             s_ti <- s_t[, i, drop = FALSE]
//             taste_ti <- taste_t[, i, drop = FALSE]
//             xs <- crossprod(s_ti, Xp_t)
//             xs <- Xp_t - matrix(rep(1, nrow(Xp_t))) %*% xs
//             xs <- matrix(rep(1, nrow(Xp_t))) %*% t(taste_ti) * xs
//             xs <- s_ti %*% matrix(rep(1, ncol(Xp_t)), nrow = 1) * xs
//             return(xs)
//           }
//           share_derivatives_wrt_theta_nonlinear_t <-
//             share_derivatives_wrt_theta_nonlinear_t %>%
//             purrr::reduce(., `+`)
//             share_derivatives_wrt_theta_nonlinear_t <-
//               share_derivatives_wrt_theta_nonlinear_t / ncol(s_t)
//             return(share_derivatives_wrt_theta_nonlinear_t)
//       }
//       return(share_derivatives_wrt_theta_nonlinear)
//   }
// [[Rcpp::export]]
Rcpp::List compute_share_derivatives_wrt_theta_nonlinear_rcpp(
  Rcpp::List individual_share,
  Rcpp::List X,
  Rcpp::List p,
  Rcpp::List nu,
  Rcpp::List upsilon
) {
  Rcpp::List share_derivatives_wrt_theta_nonlinear;
  for (int t = 0; t < individual_share.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> s_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (individual_share.at(t)));
    Eigen::Map<Eigen::MatrixXd> X_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (X.at(t)));
    Eigen::Map<Eigen::MatrixXd> p_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (p.at(t)));
    Eigen::Map<Eigen::MatrixXd> nu_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (nu.at(t)));
    Eigen::Map<Eigen::MatrixXd> upsilon_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (upsilon.at(t)));
    Eigen::MatrixXd Xp_t(X_t.rows(), X_t.cols() + p_t.cols());
    Xp_t << X_t, p_t;
    Eigen::MatrixXd taste_t(nu_t.rows() + upsilon_t.rows(), nu_t.cols());
    taste_t << nu_t, upsilon_t;
    Eigen::MatrixXd share_derivatives_wrt_theta_nonlinear_t = Eigen::MatrixXd::Zero(Xp_t.rows(), Xp_t.cols());
    for (int i = 0; i < s_t.cols(); i++) {
      Eigen::MatrixXd s_ti = s_t.col(i);
      Eigen::MatrixXd taste_ti = taste_t.col(i);
      Eigen::MatrixXd xs = s_ti.transpose() * Xp_t;
      xs = Xp_t - Eigen::VectorXd::Ones(Xp_t.rows()) * xs;
      xs = ((Eigen::VectorXd::Ones(Xp_t.rows()) * taste_ti.transpose()).array() * xs.array()).matrix();
      xs = ((s_ti * Eigen::MatrixXd::Ones(1, Xp_t.cols())).array() * xs.array()).matrix();
      share_derivatives_wrt_theta_nonlinear_t = share_derivatives_wrt_theta_nonlinear_t + xs;
    }
    share_derivatives_wrt_theta_nonlinear_t = share_derivatives_wrt_theta_nonlinear_t / s_t.cols();
    share_derivatives_wrt_theta_nonlinear.push_back(share_derivatives_wrt_theta_nonlinear_t);
  }
  return(share_derivatives_wrt_theta_nonlinear);
}

// # compute derivatives of mean utility with respect to non-inear parameters
// compute_mean_utility_derivatives_wrt_theta_nonlinear <-
//   function(individual_share, X, p, nu, upsilon) {
// # compute derivatives of share with respect to delta
//     share_derivatives_wrt_mean_utility <- compute_share_derivatives_wrt_mean_utility(individual_share)
// # compute derivatives of share with respect to non-linear parameters
//     share_derivatives_wrt_theta_nonlinear <- compute_share_derivatives_wrt_theta_nonlinear(individual_share, X, p, nu, upsilon)
// # derivatives
//     mean_utility_derivatives_wrt_theta_nonlinear <-
//       purrr::map2(share_derivatives_wrt_mean_utility, share_derivatives_wrt_theta_nonlinear, ~ - qr.solve(.x, .y))
//     return(mean_utility_derivatives_wrt_theta_nonlinear)
//   }
// [[Rcpp::export]]
Rcpp::List compute_mean_utility_derivatives_wrt_theta_nonlinear_rcpp(
  Rcpp::List individual_share,
  Rcpp::List X,
  Rcpp::List p,
  Rcpp::List nu,
  Rcpp::List upsilon
) {
  // compute derivatives of share with respect to delta
  Rcpp::List share_derivatives_wrt_mean_utility = compute_share_derivatives_wrt_mean_utility_rcpp(individual_share);
  // compute derivatives of share with respect to non-linear parameters
  Rcpp::List share_derivatives_wrt_theta_nonlinear = compute_share_derivatives_wrt_theta_nonlinear_rcpp(individual_share, X, p, nu, upsilon);
  // derivatives
  Rcpp::List mean_utility_derivatives_wrt_theta_nonlinear;
  for (int t = 0; t < share_derivatives_wrt_mean_utility.size(); t++) {
    Eigen::Map<Eigen::MatrixXd> share_derivatives_wrt_mean_utility_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (share_derivatives_wrt_mean_utility.at(t)));
    Eigen::Map<Eigen::MatrixXd> share_derivatives_wrt_theta_nonlinear_t(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (share_derivatives_wrt_theta_nonlinear.at(t)));
    Eigen::MatrixXd mean_utility_derivatives_wrt_theta_nonlinear_t = - share_derivatives_wrt_mean_utility_t.colPivHouseholderQr().solve(share_derivatives_wrt_theta_nonlinear_t);
    mean_utility_derivatives_wrt_theta_nonlinear.push_back(mean_utility_derivatives_wrt_theta_nonlinear_t);
  }
  return(mean_utility_derivatives_wrt_theta_nonlinear);
}

























