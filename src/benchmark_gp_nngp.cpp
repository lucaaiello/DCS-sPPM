// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <boost/math/special_functions/bessel.hpp>
#include <chrono>
#include <cmath>
#include <vector>

using namespace Rcpp;

double benchmark_matern_correlation(double distance, double range, double nu) {
  if (distance == 0.0) return 1.0;
  double argument = std::sqrt(2.0 * nu) * distance / range;
  return std::pow(2.0, 1.0 - nu) / R::gammafn(nu) *
    std::pow(argument, nu) * boost::math::cyl_bessel_k(nu, argument);
}

// [[Rcpp::export]]
List benchmark_gp_nngp_cpp(
    const arma::mat& locations,
    const arma::umat& neighbor_indices,
    const arma::vec& field,
    std::string method,
    double variance = 1.0,
    double range = 0.2,
    double nu = 0.5
) {
  const int n = static_cast<int>(locations.n_rows);
  if (locations.n_cols != 2 || field.n_elem != static_cast<arma::uword>(n)) {
    stop("Invalid locations or field dimensions.");
  }
  auto setup_start = std::chrono::steady_clock::now();
  double log_determinant = 0.0;
  double quadratic = 0.0;
  double storage_bytes = 0.0;

  if (method == "exact_gp") {
    arma::mat covariance(n, n, arma::fill::zeros);
    for (int row = 0; row < n; ++row) {
      covariance(row, row) = variance;
      for (int column = row + 1; column < n; ++column) {
        double distance = arma::norm(
          locations.row(row) - locations.row(column), 2
        );
        double value = variance * benchmark_matern_correlation(
          distance, range, nu
        );
        covariance(row, column) = value;
        covariance(column, row) = value;
      }
    }
    arma::mat cholesky;
    bool factorized = false;
    for (int attempt = 0; attempt < 7; ++attempt) {
      arma::mat candidate = covariance;
      if (attempt > 0) {
        candidate.diag() += 1e-10 * variance * std::pow(10.0, attempt - 1);
      }
      factorized = arma::chol(cholesky, candidate, "lower");
      if (factorized) break;
    }
    if (!factorized) stop("Exact GP covariance factorization failed.");
    auto setup_end = std::chrono::steady_clock::now();
    auto evaluation_start = setup_end;
    arma::vec standardized = arma::solve(
      arma::trimatl(cholesky), field, arma::solve_opts::fast
    );
    quadratic = arma::dot(standardized, standardized);
    log_determinant = 2.0 * arma::sum(arma::log(cholesky.diag()));
    auto evaluation_end = std::chrono::steady_clock::now();
    storage_bytes = 8.0 * static_cast<double>(n) * static_cast<double>(n) * 2.0;
    return List::create(
      Named("setup_seconds") = std::chrono::duration<double>(setup_end - setup_start).count(),
      Named("evaluation_seconds") = std::chrono::duration<double>(evaluation_end - evaluation_start).count(),
      Named("total_seconds") = std::chrono::duration<double>(evaluation_end - setup_start).count(),
      Named("log_density") = -0.5 * (n * std::log(2.0 * arma::datum::pi) + log_determinant + quadratic),
      Named("storage_mb") = storage_bytes / (1024.0 * 1024.0)
    );
  }

  if (method != "nngp") stop("Method must be `exact_gp` or `nngp`.");
  if (neighbor_indices.n_rows != static_cast<arma::uword>(n)) {
    stop("NNGP neighbor matrix must have one row per location.");
  }

  std::vector<arma::uvec> neighbors(n);
  std::vector<arma::vec> coefficients(n);
  arma::vec conditional_variances(n, arma::fill::value(variance));
  arma::uword coefficient_count = 0;
  for (int site = 0; site < n; ++site) {
    arma::uvec supplied = neighbor_indices.row(site).t();
    supplied = supplied.elem(arma::find(supplied > 0));
    if (supplied.n_elem == 0) continue;
    neighbors[site] = supplied - 1;
    int neighbor_count = static_cast<int>(neighbors[site].n_elem);
    arma::vec cross_covariance(neighbor_count);
    arma::mat neighbor_covariance(neighbor_count, neighbor_count);
    for (int first = 0; first < neighbor_count; ++first) {
      double distance = arma::norm(
        locations.row(site) - locations.row(neighbors[site](first)), 2
      );
      cross_covariance(first) = variance * benchmark_matern_correlation(
        distance, range, nu
      );
      for (int second = first; second < neighbor_count; ++second) {
        double neighbor_distance = arma::norm(
          locations.row(neighbors[site](first)) -
            locations.row(neighbors[site](second)),
          2
        );
        double value = variance * benchmark_matern_correlation(
          neighbor_distance, range, nu
        );
        neighbor_covariance(first, second) = value;
        neighbor_covariance(second, first) = value;
      }
    }
    bool solved = arma::solve(
      coefficients[site],
      neighbor_covariance,
      cross_covariance,
      arma::solve_opts::likely_sympd
    );
    if (!solved || !coefficients[site].is_finite()) {
      stop("NNGP neighbor solve failed at site %d.", site + 1);
    }
    conditional_variances(site) = variance -
      arma::dot(cross_covariance, coefficients[site]);
    if (!std::isfinite(conditional_variances(site)) ||
        conditional_variances(site) <= 0.0) {
      stop("Invalid NNGP conditional variance at site %d.", site + 1);
    }
    coefficient_count += coefficients[site].n_elem;
  }
  auto setup_end = std::chrono::steady_clock::now();
  auto evaluation_start = setup_end;
  for (int site = 0; site < n; ++site) {
    double innovation = field(site);
    if (neighbors[site].n_elem > 0) {
      innovation -= arma::dot(
        coefficients[site], field.elem(neighbors[site])
      );
    }
    quadratic += innovation * innovation / conditional_variances(site);
  }
  log_determinant = arma::sum(arma::log(conditional_variances));
  auto evaluation_end = std::chrono::steady_clock::now();
  storage_bytes = 8.0 * static_cast<double>(n + coefficient_count) +
    4.0 * static_cast<double>(coefficient_count);
  return List::create(
    Named("setup_seconds") = std::chrono::duration<double>(setup_end - setup_start).count(),
    Named("evaluation_seconds") = std::chrono::duration<double>(evaluation_end - evaluation_start).count(),
    Named("total_seconds") = std::chrono::duration<double>(evaluation_end - setup_start).count(),
    Named("log_density") = -0.5 * (n * std::log(2.0 * arma::datum::pi) + log_determinant + quadratic),
    Named("storage_mb") = storage_bytes / (1024.0 * 1024.0)
  );
}
