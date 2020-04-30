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
//         denominator_t <- matrix(apply(1 + individual_share_t, 2, sum), nrow = 1)
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












