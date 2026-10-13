// DCS-sPPM sampler with global beta and memory-efficient spatial prediction.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Matrix)]]
// [[Rcpp::plugins(openmp)]]

// #include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <boost/math/special_functions/bessel.hpp>
#include <chrono>
#include <cctype>
#include <cmath>
#include <deque>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <RcppArmadillo.h>
#include <sstream>
#include <vector>
using namespace Rcpp;
#include <Rcpp/Benchmark/Timer.h>
#include <Matrix.h>
#include <cholmod.h>
#ifdef _OPENMP
#include <omp.h>
#endif
extern "C" {
#include <Matrix_stubs.c>
}

////////////////////////////////////////////////////////////////////////////////
// CHOLMOD sampling utilities
////////////////////////////////////////////////////////////////////////////////

/**
 * ==============================================================================
 * CLASS: CholmodSampler
 * ==============================================================================
 * PURPOSE: 
 * A high-performance, stateful wrapper for drawing multivariate normal 
 * samples using CHOLMOD's sparse Cholesky factorization.
 * * WHY THIS IS NECESSARY (THE OPTIMIZATION):
 * In a standard MCMC loop, the spatial/temporal network structure (the locations 
 * of the zeros and non-zeros in the precision matrix Q) never changes. Only the 
 * actual variance values change. 
 * * Sparse Cholesky factorization operates in two distinct phases:
 * 1. Symbolic Factorization (analyze): A highly expensive graph-theory 
 * algorithm that finds the optimal permutation to minimize memory fill-in.
 * 2. Numeric Factorization (factorize): The actual floating-point math.
 * * Previous inline approaches recalculated the Symbolic map and re-allocated the 
 * CHOLMOD environment on every single MCMC iteration, creating a massive bottleneck. 
 * This class uses Object-Oriented design to create Persistent State:
 * * HOW TO USE:
 * - Step 1: Instantiate the object and call `analyze_pattern(Q)` exactly ONCE 
 * before the MCMC loop begins. This permanently stores the optimal 
 * sparsity map in memory.
 * - Step 2: Pass the object by reference into your update functions. Call 
 * `.sample(Q, b)` inside the loop. This skips the symbolic phase 
 * entirely and executes blazing-fast numeric-only factorization.
 * ==============================================================================
 */

class CholmodSampler {
private:
  cholmod_common c;
  cholmod_factor* L;
  bool is_initialized;
  arma::vec mean_buffer;
  arma::vec noise_buffer;
  arma::vec randn_buffer;
  
public:
  // Constructor: Starts the environment
  CholmodSampler() {
    M_cholmod_start(&c);
    c.supernodal = CHOLMOD_SIMPLICIAL;
    c.final_asis = 0;
    c.final_super = 0;
    c.final_ll = 1;
    L = nullptr;
    is_initialized = false;
  }
  
  // Destructor: Safely cleans up memory when the object is destroyed
  ~CholmodSampler() {
    if (L != nullptr) {
      M_cholmod_free_factor(&L, &c);
    }
    M_cholmod_finish(&c);
  }
  
  // Step 1: Run this ONCE outside your MCMC loop to analyze the pattern
  void analyze_pattern(const arma::sp_mat& Q) {
    cholmod_sparse A;
    A.nrow = Q.n_rows; A.ncol = Q.n_cols; A.nzmax = Q.n_nonzero;
    A.p = (void*) Q.col_ptrs; A.i = (void*) Q.row_indices; A.x = (void*) Q.values;
    A.stype = 1; A.itype = CHOLMOD_INT; A.xtype = CHOLMOD_REAL; A.dtype = CHOLMOD_DOUBLE;
    A.sorted = 1; A.packed = 1;
    
    // Perform the expensive symbolic factorization only once
    if (L != nullptr) M_cholmod_free_factor(&L, &c); 
    L = M_cholmod_analyze(&A, &c);
    is_initialized = true;
  }
  
  // Step 2: Run this INSIDE your MCMC loop
  arma::vec sample(const arma::sp_mat& Q, const arma::vec& b) {
    int n = Q.n_rows;
    
    // Map Armadillo sparse matrix to CHOLMOD (zero copy)
    cholmod_sparse A;
    A.nrow = Q.n_rows; A.ncol = Q.n_cols; A.nzmax = Q.n_nonzero;
    A.p = (void*) Q.col_ptrs; A.i = (void*) Q.row_indices; A.x = (void*) Q.values;
    A.stype = 1; A.itype = CHOLMOD_INT; A.xtype = CHOLMOD_REAL; A.dtype = CHOLMOD_DOUBLE;
    A.sorted = 1; A.packed = 1;
    
    // Fast Numeric Factorization (re-uses the symbolic pattern in L)
    M_cholmod_factorize(&A, L, &c);
    
    // Setup dense vector b
    cholmod_dense b_dense;
    b_dense.nrow = n; b_dense.ncol = 1; b_dense.d = n;
    b_dense.x = (void*) b.memptr(); b_dense.xtype = CHOLMOD_REAL; b_dense.dtype = CHOLMOD_DOUBLE;
    
    // Solve for mu: Q * mu = b
    cholmod_dense *mu_ptr = M_cholmod_solve(CHOLMOD_A, L, &b_dense, &c);
    mean_buffer.set_size(n);
    std::memcpy(mean_buffer.memptr(), mu_ptr->x, n * sizeof(double));
    
    // Generate noise
    randn_buffer = arma::randn(n);
    cholmod_dense z_dense;
    z_dense.nrow = n; z_dense.ncol = 1; z_dense.d = n;
    z_dense.x = (void*) randn_buffer.memptr(); z_dense.xtype = CHOLMOD_REAL; z_dense.dtype = CHOLMOD_DOUBLE;
    
    // For a CHOLMOD factorization of P A P' = L L', a valid draw from A^{-1}
    // is P' * L'^{-1} z with z ~ N(0, I). The previous code used only the
    // L' solve and missed the permutation, which inflates the covariance.
    cholmod_dense *v_ptr = M_cholmod_solve(CHOLMOD_Lt, L, &z_dense, &c);
    cholmod_dense *noise_ptr = M_cholmod_solve(CHOLMOD_Pt, L, v_ptr, &c);
    noise_buffer.set_size(n);
    std::memcpy(noise_buffer.memptr(), noise_ptr->x, n * sizeof(double));
    
    // Free only the temporary dense vectors
    M_cholmod_free_dense(&mu_ptr, &c);
    M_cholmod_free_dense(&v_ptr, &c);
    M_cholmod_free_dense(&noise_ptr, &c);
    
    return mean_buffer + noise_buffer;
  }

  arma::vec sample_factorized(const arma::vec& b) {
    if (!is_initialized) {
      Rcpp::stop("CholmodSampler::analyze_pattern must be called before sample_factorized.");
    }

    int n = static_cast<int>(b.n_elem);

    cholmod_dense b_dense;
    b_dense.nrow = n; b_dense.ncol = 1; b_dense.d = n;
    b_dense.x = (void*) b.memptr(); b_dense.xtype = CHOLMOD_REAL; b_dense.dtype = CHOLMOD_DOUBLE;

    cholmod_dense *mu_ptr = M_cholmod_solve(CHOLMOD_A, L, &b_dense, &c);
    mean_buffer.set_size(n);
    std::memcpy(mean_buffer.memptr(), mu_ptr->x, n * sizeof(double));

    randn_buffer = arma::randn(n);
    cholmod_dense z_dense;
    z_dense.nrow = n; z_dense.ncol = 1; z_dense.d = n;
    z_dense.x = (void*) randn_buffer.memptr(); z_dense.xtype = CHOLMOD_REAL; z_dense.dtype = CHOLMOD_DOUBLE;

    cholmod_dense *v_ptr = M_cholmod_solve(CHOLMOD_Lt, L, &z_dense, &c);
    cholmod_dense *noise_ptr = M_cholmod_solve(CHOLMOD_Pt, L, v_ptr, &c);
    noise_buffer.set_size(n);
    std::memcpy(noise_buffer.memptr(), noise_ptr->x, n * sizeof(double));

    M_cholmod_free_dense(&mu_ptr, &c);
    M_cholmod_free_dense(&v_ptr, &c);
    M_cholmod_free_dense(&noise_ptr, &c);

    return mean_buffer + noise_buffer;
  }

  arma::mat sample_factorized_columns_parallel(
      const arma::mat& right_hand_sides,
      const arma::mat& standard_normals
  ) const {
    if (!is_initialized || L == nullptr) {
      Rcpp::stop(
        "CholmodSampler must be factorized before parallel column sampling."
      );
    }
    if (L->is_super || !L->is_ll || L->xtype != CHOLMOD_REAL ||
        L->itype != CHOLMOD_INT) {
      Rcpp::stop(
        "Parallel column sampling requires a simplicial real LL' factor."
      );
    }

    const int n = static_cast<int>(L->n);
    const int n_columns = static_cast<int>(right_hand_sides.n_cols);
    if (right_hand_sides.n_rows != static_cast<arma::uword>(n) ||
        standard_normals.n_rows != static_cast<arma::uword>(n) ||
        standard_normals.n_cols != right_hand_sides.n_cols) {
      Rcpp::stop("Invalid dimensions in parallel column sampling.");
    }

    const int* column_ptr = static_cast<const int*>(L->p);
    const int* row_index = static_cast<const int*>(L->i);
    const int* column_nnz = static_cast<const int*>(L->nz);
    const int* permutation = static_cast<const int*>(L->Perm);
    const double* factor_values = static_cast<const double*>(L->x);

    std::vector<int> column_end(n);
    std::vector<int> diagonal_position(n, -1);
    for (int column = 0; column < n; ++column) {
      int start = column_ptr[column];
      int end = column_nnz == nullptr ? column_ptr[column + 1] :
        start + column_nnz[column];
      column_end[column] = end;
      for (int position = start; position < end; ++position) {
        if (row_index[position] == column) {
          diagonal_position[column] = position;
          break;
        }
      }
      if (diagonal_position[column] < 0 ||
          !std::isfinite(factor_values[diagonal_position[column]]) ||
          factor_values[diagonal_position[column]] <= 0.0) {
        Rcpp::stop("Invalid diagonal in sparse Cholesky factor.");
      }
    }

    arma::mat draws(n, n_columns, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int draw_index = 0; draw_index < n_columns; ++draw_index) {
      std::vector<double> mean_work(n);
      std::vector<double> noise_work(n);

      for (int permuted_index = 0; permuted_index < n; ++permuted_index) {
        int original_index = permutation == nullptr ?
          permuted_index : permutation[permuted_index];
        mean_work[permuted_index] =
          right_hand_sides(original_index, draw_index);
        noise_work[permuted_index] =
          standard_normals(permuted_index, draw_index);
      }

      for (int column = 0; column < n; ++column) {
        int diagonal = diagonal_position[column];
        mean_work[column] /= factor_values[diagonal];
        for (int position = column_ptr[column];
             position < column_end[column]; ++position) {
          int row = row_index[position];
          if (row > column) {
            mean_work[row] -= factor_values[position] * mean_work[column];
          }
        }
      }

      for (int column = n - 1; column >= 0; --column) {
        double mean_value = mean_work[column];
        double noise_value = noise_work[column];
        for (int position = column_ptr[column];
             position < column_end[column]; ++position) {
          int row = row_index[position];
          if (row > column) {
            mean_value -= factor_values[position] * mean_work[row];
            noise_value -= factor_values[position] * noise_work[row];
          }
        }
        double diagonal = factor_values[diagonal_position[column]];
        mean_work[column] = mean_value / diagonal;
        noise_work[column] = noise_value / diagonal;
      }

      for (int permuted_index = 0; permuted_index < n; ++permuted_index) {
        int original_index = permutation == nullptr ?
          permuted_index : permutation[permuted_index];
        draws(original_index, draw_index) =
          mean_work[permuted_index] + noise_work[permuted_index];
      }
    }

    return draws;
  }

  arma::vec solve_mean(const arma::sp_mat& Q, const arma::vec& b) {
    int n = Q.n_rows;

    cholmod_sparse A;
    A.nrow = Q.n_rows; A.ncol = Q.n_cols; A.nzmax = Q.n_nonzero;
    A.p = (void*) Q.col_ptrs; A.i = (void*) Q.row_indices; A.x = (void*) Q.values;
    A.stype = 1; A.itype = CHOLMOD_INT; A.xtype = CHOLMOD_REAL; A.dtype = CHOLMOD_DOUBLE;
    A.sorted = 1; A.packed = 1;

    M_cholmod_factorize(&A, L, &c);

    cholmod_dense b_dense;
    b_dense.nrow = n; b_dense.ncol = 1; b_dense.d = n;
    b_dense.x = (void*) b.memptr(); b_dense.xtype = CHOLMOD_REAL; b_dense.dtype = CHOLMOD_DOUBLE;

    cholmod_dense *mu_ptr = M_cholmod_solve(CHOLMOD_A, L, &b_dense, &c);
    mean_buffer.set_size(n);
    std::memcpy(mean_buffer.memptr(), mu_ptr->x, n * sizeof(double));
    M_cholmod_free_dense(&mu_ptr, &c);

    return mean_buffer;
  }

  void factorize_numeric(const arma::sp_mat& Q) {
    if (!is_initialized) {
      Rcpp::stop("CholmodSampler::analyze_pattern must be called before factorize_numeric.");
    }

    cholmod_sparse A;
    A.nrow = Q.n_rows; A.ncol = Q.n_cols; A.nzmax = Q.n_nonzero;
    A.p = (void*) Q.col_ptrs; A.i = (void*) Q.row_indices; A.x = (void*) Q.values;
    A.stype = 1; A.itype = CHOLMOD_INT; A.xtype = CHOLMOD_REAL; A.dtype = CHOLMOD_DOUBLE;
    A.sorted = 1; A.packed = 1;

    M_cholmod_factorize(&A, L, &c);
  }

  arma::vec solve_mean_factorized(const arma::vec& b) {
    int n = static_cast<int>(b.n_elem);

    cholmod_dense b_dense;
    b_dense.nrow = n; b_dense.ncol = 1; b_dense.d = n;
    b_dense.x = (void*) b.memptr(); b_dense.xtype = CHOLMOD_REAL; b_dense.dtype = CHOLMOD_DOUBLE;

    cholmod_dense *mu_ptr = M_cholmod_solve(CHOLMOD_A, L, &b_dense, &c);
    mean_buffer.set_size(n);
    std::memcpy(mean_buffer.memptr(), mu_ptr->x, n * sizeof(double));
    M_cholmod_free_dense(&mu_ptr, &c);

    return mean_buffer;
  }

  arma::mat solve_means_factorized(const arma::mat& B) {
    int n = B.n_rows;
    int m = B.n_cols;

    cholmod_dense b_dense;
    b_dense.nrow = n; b_dense.ncol = m; b_dense.d = n;
    b_dense.x = (void*) B.memptr(); b_dense.xtype = CHOLMOD_REAL; b_dense.dtype = CHOLMOD_DOUBLE;

    cholmod_dense *mu_ptr = M_cholmod_solve(CHOLMOD_A, L, &b_dense, &c);
    arma::mat solutions(n, m);
    std::memcpy(solutions.memptr(), mu_ptr->x, static_cast<std::size_t>(n) * static_cast<std::size_t>(m) * sizeof(double));
    M_cholmod_free_dense(&mu_ptr, &c);

    return solutions;
  }

  double log_determinant_factorized() const {
    if (!is_initialized || L == nullptr) {
      Rcpp::stop("CholmodSampler must be factorized before requesting its log determinant.");
    }
    if (L->is_super || !L->is_ll || L->xtype != CHOLMOD_REAL ||
        L->itype != CHOLMOD_INT) {
      Rcpp::stop("Unsupported CHOLMOD factor representation for log determinant.");
    }

    const int* column_ptr = static_cast<const int*>(L->p);
    const int* row_index = static_cast<const int*>(L->i);
    const int* column_nnz = static_cast<const int*>(L->nz);
    const double* factor_values = static_cast<const double*>(L->x);
    double log_determinant = 0.0;

    for (std::size_t column = 0; column < L->n; ++column) {
      int start = column_ptr[column];
      int end = column_nnz == nullptr ? column_ptr[column + 1] :
        start + column_nnz[column];
      double diagonal = NA_REAL;
      for (int position = start; position < end; ++position) {
        if (row_index[position] == static_cast<int>(column)) {
          diagonal = factor_values[position];
          break;
        }
      }
      if (!std::isfinite(diagonal) || diagonal <= 0.0) {
        Rcpp::stop("Invalid diagonal in CHOLMOD factor.");
      }
      log_determinant += 2.0 * std::log(diagonal);
    }
    return log_determinant;
  }
};

////////////////////////////////////////////////////////////////////////////////
// Basic parameter updates
////////////////////////////////////////////////////////////////////////////////

double sigma2_update_global(const arma::vec& y, const arma::vec& w, const arma::vec& Zbeta,
                            int n, int T, double a_sigma, double b_sigma) {
  arma::vec residuals = y - Zbeta - w;
  double ss = arma::dot(residuals, residuals);
  double shape = a_sigma + 0.5 * static_cast<double>(n * T);
  double scale = b_sigma + 0.5 * ss;
  return 1.0 / arma::randg<double>(arma::distr_param(shape, 1.0 / scale));
}

bool use_collapsed_sigma2_update(std::string sigma2_update_mode) {
  std::transform(
    sigma2_update_mode.begin(), sigma2_update_mode.end(), sigma2_update_mode.begin(),
    [](unsigned char character) { return std::tolower(character); }
  );
  if (sigma2_update_mode == "collapsed") {
    return true;
  }
  if (sigma2_update_mode == "conditional") {
    return false;
  }
  Rcpp::stop("`sigma2_update_mode` must be either 'collapsed' or 'conditional'.");
  return false;
}

bool use_collapsed_spatial_update(std::string spatial_update_mode) {
  std::transform(
    spatial_update_mode.begin(), spatial_update_mode.end(), spatial_update_mode.begin(),
    [](unsigned char character) { return std::tolower(character); }
  );
  if (spatial_update_mode == "collapsed") {
    return true;
  }
  if (spatial_update_mode == "conditional") {
    return false;
  }
  Rcpp::stop("`spatial_update_mode` must be either 'collapsed' or 'conditional'.");
  return false;
}

int tau2_update_code(std::string tau2_update_mode) {
  std::transform(
    tau2_update_mode.begin(), tau2_update_mode.end(), tau2_update_mode.begin(),
    [](unsigned char character) { return std::tolower(character); }
  );
  if (tau2_update_mode == "collapsed") {
    return 1;
  }
  if (tau2_update_mode == "conditional") {
    return 0;
  }
  if (tau2_update_mode == "joint_collapsed") {
    return 2;
  }
  if (tau2_update_mode == "pairwise_collapsed") {
    return 3;
  }
  if (tau2_update_mode == "pairwise_interwoven") {
    return 4;
  }
  Rcpp::stop(
    "`tau2_update_mode` must be 'conditional', 'collapsed', "
    "'joint_collapsed', 'pairwise_collapsed', or 'pairwise_interwoven'."
  );
  return 0;
}

struct RunningCovariance {
  int count = 0;
  int proposals = 0;
  double log_scale = 0.0;
  arma::vec mean;
  arma::mat centered_sum;

  void update(const arma::vec& value) {
    if (count == 0) {
      mean = value;
      centered_sum.zeros(value.n_elem, value.n_elem);
      count = 1;
      return;
    }
    count += 1;
    arma::vec delta = value - mean;
    mean += delta / static_cast<double>(count);
    centered_sum += delta * (value - mean).t();
  }

  void adapt_scale(
      double acceptance_probability,
      double target,
      double minimum_log_scale = -5.0,
      double maximum_log_scale = 5.0
  ) {
    proposals += 1;
    double gain = std::min(
      0.05,
      std::pow(static_cast<double>(proposals + 10), -0.6)
    );
    log_scale += gain * (acceptance_probability - target);
    log_scale = std::max(
      minimum_log_scale,
      std::min(log_scale, maximum_log_scale)
    );
  }

  arma::mat proposal_covariance(
      double initial_scale,
      double diagonal_ridge = 1e-6
  ) const {
    arma::uword dimension = mean.n_elem;
    double scale_multiplier = std::exp(2.0 * log_scale);
    if (count <= static_cast<int>(2 * dimension + 20)) {
      return scale_multiplier * std::pow(initial_scale, 2.0) *
        arma::eye(dimension, dimension);
    }
    arma::mat empirical = centered_sum / static_cast<double>(count - 1);
    double optimal_scale = 2.38 * 2.38 / static_cast<double>(dimension);
    return scale_multiplier * optimal_scale * empirical +
      diagonal_ridge * arma::eye(dimension, dimension);
  }
};

double collapsed_sigma2_log_target_cpp(
    double log_sigma2,
    const arma::vec& y,
    const arma::vec& Xbeta,
    const arma::vec& u,
    const arma::vec& phi,
    const arma::vec& tau2,
    int n,
    int T,
    double a_sigma,
    double b_sigma
) {
  double sigma2 = std::exp(log_sigma2);
  if (!std::isfinite(sigma2) || sigma2 <= 0.0) {
    return -arma::datum::inf;
  }

  arma::mat y_mat(const_cast<double*>(y.memptr()), n, T, false, false);
  arma::mat mean_mat(const_cast<double*>(Xbeta.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);
  double loglik = 0.0;

  for (int site = 0; site < n; ++site) {
    double phi_i = phi(site);
    double tau2_i = tau2(site);
    double state_mean = 0.0;
    double state_variance = tau2_i / (1.0 - phi_i * phi_i);

    for (int t = 0; t < T; ++t) {
      double innovation_variance = state_variance + sigma2;
      double innovation = y_mat(site, t) - mean_mat(site, t) - state_mean;
      loglik += -0.5 * (
        std::log(2.0 * arma::datum::pi * innovation_variance) +
        innovation * innovation / innovation_variance
      );

      double kalman_gain = state_variance / innovation_variance;
      state_mean += kalman_gain * innovation;
      state_variance *= (1.0 - kalman_gain);

      if (t + 1 < T) {
        state_mean = phi_i * state_mean + u_mat(site, t + 1);
        state_variance = phi_i * phi_i * state_variance + tau2_i;
      }
    }
  }

  return loglik - a_sigma * log_sigma2 - b_sigma / sigma2;
}

double collapsed_tau2_log_target_cpp(
    double log_tau2,
    const arma::uvec& members,
    const arma::vec& y,
    const arma::vec& Xbeta,
    const arma::vec& u,
    const arma::vec& phi,
    double sigma2,
    int n,
    int T,
    double a_tau,
    double b_tau
) {
  double tau2_candidate = std::exp(log_tau2);
  if (!std::isfinite(tau2_candidate) || tau2_candidate <= 0.0) {
    return -arma::datum::inf;
  }

  arma::mat y_mat(const_cast<double*>(y.memptr()), n, T, false, false);
  arma::mat mean_mat(const_cast<double*>(Xbeta.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);
  double loglik = 0.0;

  for (arma::uword member = 0; member < members.n_elem; ++member) {
    int site = static_cast<int>(members(member));
    double phi_i = phi(site);
    double state_mean = 0.0;
    double state_variance = tau2_candidate / (1.0 - phi_i * phi_i);

    for (int t = 0; t < T; ++t) {
      double innovation_variance = state_variance + sigma2;
      double innovation = y_mat(site, t) - mean_mat(site, t) - state_mean;
      loglik += -0.5 * (
        std::log(2.0 * arma::datum::pi * innovation_variance) +
        innovation * innovation / innovation_variance
      );

      double kalman_gain = state_variance / innovation_variance;
      state_mean += kalman_gain * innovation;
      state_variance *= (1.0 - kalman_gain);

      if (t + 1 < T) {
        state_mean = phi_i * state_mean + u_mat(site, t + 1);
        state_variance = phi_i * phi_i * state_variance + tau2_candidate;
      }
    }
  }

  return loglik - a_tau * log_tau2 - b_tau / tau2_candidate;
}

// [[Rcpp::export]]
double joint_collapsed_variance_log_target_cpp(
    const arma::vec& transformed_variances,
    const arma::ivec& cluster,
    const arma::vec& y,
    const arma::vec& Xbeta,
    const arma::vec& u,
    const arma::vec& phi,
    int n,
    int T,
    double a_sigma,
    double b_sigma,
    double a_tau,
    double b_tau
) {
  int n_clusters = transformed_variances.n_elem - 1;
  if (n_clusters < 1 || cluster.n_elem != static_cast<arma::uword>(n)) {
    return -arma::datum::inf;
  }

  double log_sigma2 = transformed_variances(0);
  double sigma2 = std::exp(log_sigma2);
  arma::vec tau2_star(n_clusters);
  for (int k = 0; k < n_clusters; ++k) {
    tau2_star(k) = std::exp(log_sigma2 + transformed_variances(k + 1));
  }
  if (!std::isfinite(sigma2) || sigma2 <= 0.0 || !tau2_star.is_finite() ||
      arma::any(tau2_star <= 0.0)) {
    return -arma::datum::inf;
  }

  arma::mat y_mat(const_cast<double*>(y.memptr()), n, T, false, false);
  arma::mat mean_mat(const_cast<double*>(Xbeta.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);
  double log_target = -a_sigma * log_sigma2 - b_sigma / sigma2;
  for (int k = 0; k < n_clusters; ++k) {
    double log_tau2 = log_sigma2 + transformed_variances(k + 1);
    log_target += -a_tau * log_tau2 - b_tau / tau2_star(k);
  }

  for (int site = 0; site < n; ++site) {
    int cluster_index = cluster(site) - 1;
    if (cluster_index < 0 || cluster_index >= n_clusters) {
      return -arma::datum::inf;
    }
    double phi_i = phi(site);
    double tau2_i = tau2_star(cluster_index);
    double state_mean = 0.0;
    double state_variance = tau2_i / (1.0 - phi_i * phi_i);

    for (int t = 0; t < T; ++t) {
      double innovation_variance = state_variance + sigma2;
      double innovation = y_mat(site, t) - mean_mat(site, t) - state_mean;
      log_target += -0.5 * (
        std::log(2.0 * arma::datum::pi * innovation_variance) +
        innovation * innovation / innovation_variance
      );
      double kalman_gain = state_variance / innovation_variance;
      state_mean += kalman_gain * innovation;
      state_variance *= (1.0 - kalman_gain);
      if (t + 1 < T) {
        state_mean = phi_i * state_mean + u_mat(site, t + 1);
        state_variance = phi_i * phi_i * state_variance + tau2_i;
      }
    }
  }
  return log_target;
}

double collapsed_variance_loglik_parallel_cpp(
    double sigma2,
    const arma::vec& tau2_star,
    const arma::ivec& cluster,
    const arma::vec& y,
    const arma::vec& Xbeta,
    const arma::vec& u,
    const arma::vec& phi,
    int n,
    int T
) {
  arma::mat y_mat(const_cast<double*>(y.memptr()), n, T, false, false);
  arma::mat mean_mat(const_cast<double*>(Xbeta.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);
  double loglik = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:loglik) schedule(static)
#endif
  for (int site = 0; site < n; ++site) {
    int site_cluster = cluster(site) - 1;
    double phi_i = phi(site);
    double tau2_i = tau2_star(site_cluster);
    double state_mean = 0.0;
    double state_variance = tau2_i / (1.0 - phi_i * phi_i);
    double site_loglik = 0.0;
    for (int t = 0; t < T; ++t) {
      double innovation_variance = state_variance + sigma2;
      double innovation = y_mat(site, t) - mean_mat(site, t) - state_mean;
      site_loglik += -0.5 * (
        std::log(2.0 * arma::datum::pi * innovation_variance) +
        innovation * innovation / innovation_variance
      );
      double kalman_gain = state_variance / innovation_variance;
      state_mean += kalman_gain * innovation;
      state_variance *= (1.0 - kalman_gain);
      if (t + 1 < T) {
        state_mean = phi_i * state_mean + u_mat(site, t + 1);
        state_variance = phi_i * phi_i * state_variance + tau2_i;
      }
    }
    loglik += site_loglik;
  }
  return loglik;
}

// [[Rcpp::export]]
double pairwise_collapsed_variance_log_target_cpp(
    const arma::vec& transformed_pair,
    int cluster_to_update,
    const arma::vec& tau2_star_current,
    const arma::ivec& cluster,
    const arma::vec& y,
    const arma::vec& Xbeta,
    const arma::vec& u,
    const arma::vec& phi,
    int n,
    int T,
    double a_sigma,
    double b_sigma,
    double a_tau,
    double b_tau
) {
  int cluster_index = cluster_to_update - 1;
  if (transformed_pair.n_elem != 2 || cluster_index < 0 ||
      cluster_index >= static_cast<int>(tau2_star_current.n_elem) ||
      cluster.n_elem != static_cast<arma::uword>(n)) {
    return -arma::datum::inf;
  }

  double log_sigma2 = transformed_pair(0);
  double sigma2 = std::exp(log_sigma2);
  double log_tau2_candidate = log_sigma2 + transformed_pair(1);
  double tau2_candidate = std::exp(log_tau2_candidate);
  if (!std::isfinite(sigma2) || sigma2 <= 0.0 ||
      !std::isfinite(tau2_candidate) || tau2_candidate <= 0.0 ||
      !tau2_star_current.is_finite() || arma::any(tau2_star_current <= 0.0)) {
    return -arma::datum::inf;
  }
  if (arma::any(cluster < 1) ||
      arma::any(cluster > static_cast<int>(tau2_star_current.n_elem))) {
    return -arma::datum::inf;
  }

  arma::vec tau2_star = tau2_star_current;
  tau2_star(cluster_index) = tau2_candidate;
  double log_target =
    -a_sigma * log_sigma2 - b_sigma / sigma2 -
    a_tau * log_tau2_candidate - b_tau / tau2_candidate;
  return log_target + collapsed_variance_loglik_parallel_cpp(
    sigma2, tau2_star, cluster, y, Xbeta, u, phi, n, T
  );
}

void update_missing_y_cpp(
    arma::vec& y,
    const arma::uvec& missing_indices,
    const arma::vec& w,
    const arma::vec& Zbeta,
    double sigma2
) {
  if (missing_indices.n_elem == 0) {
    return;
  }

  y.elem(missing_indices) =
    Zbeta.elem(missing_indices) +
    w.elem(missing_indices) +
    std::sqrt(sigma2) * arma::randn<arma::vec>(missing_indices.n_elem);
}

double observed_gaussian_loglik_cpp(
    const arma::vec& y_observed,
    const arma::uvec& observed_indices,
    const arma::vec& w,
    const arma::vec& Zbeta,
    double sigma2
) {
  if (observed_indices.n_elem == 0) {
    return 0.0;
  }

  arma::vec residuals =
    y_observed.elem(observed_indices) -
    Zbeta.elem(observed_indices) -
    w.elem(observed_indices);
  double count = static_cast<double>(observed_indices.n_elem);
  return -0.5 * count * std::log(2.0 * arma::datum::pi * sigma2) -
    0.5 * arma::dot(residuals, residuals) / sigma2;
}

////////////////////////////////////////////////////////////////////////////////
// Likelihood and NNGP cache helpers
////////////////////////////////////////////////////////////////////////////////

arma::rowvec compute_w_tilde_row_cpp(
    const arma::vec& w,
    const arma::vec& phi,
    arma::uword row_idx,
    int n,
    int T,
    double phi_override = arma::datum::nan
) {
  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::rowvec row_tilde(T);
  double phi_row = std::isnan(phi_override) ? phi(row_idx) : phi_override;

  row_tilde(0) = std::sqrt(1.0 - phi_row * phi_row) * w_mat(row_idx, 0);
  for (int t = 1; t < T; ++t) {
    row_tilde(t) = w_mat(row_idx, t) - phi_row * w_mat(row_idx, t - 1);
  }

  return row_tilde;
}

void build_w_tilde_cache_cpp(
    const arma::vec& w,
    const arma::vec& phi,
    int n,
    int T,
    arma::mat& w_tilde
) {
  w_tilde.set_size(n, T);
  for (int i = 0; i < n; ++i) {
    w_tilde.row(i) = compute_w_tilde_row_cpp(w, phi, i, n, T);
  }
}

double compute_log_phi_term_cpp(const arma::vec& phi) {
  return 0.5 * arma::sum(arma::log(1.0 - arma::square(phi)));
}

arma::rowvec compute_innovation_row_cpp(
    arma::uword row_idx,
    const arma::mat& w_tilde,
    const List& precomp,
    const arma::vec& values,
    const arma::rowvec* updated_row_tilde = nullptr,
    arma::uword updated_idx = 0
) {
  arma::uvec row_offsets = as<arma::uvec>(precomp["row_offsets"]);
  List neighbors_list = as<List>(precomp["neighbors"]);
  arma::uvec neighbors = as<arma::uvec>(neighbors_list[row_idx]);

  arma::rowvec innovation = (updated_row_tilde != nullptr && row_idx == updated_idx)
    ? *updated_row_tilde
    : w_tilde.row(row_idx);

  arma::uword start = row_offsets(row_idx);
  for (arma::uword j = 0; j < neighbors.n_elem; ++j) {
    const arma::rowvec& neighbor_row = (updated_row_tilde != nullptr && neighbors(j) == updated_idx)
      ? *updated_row_tilde
      : w_tilde.row(neighbors(j));
    innovation -= values(start + j) * neighbor_row;
  }

  return innovation;
}

void build_innovation_cache_cpp(
    const arma::mat& w_tilde,
    const List& precomp,
    const arma::vec& values,
    arma::mat& innovations
) {
  int n = w_tilde.n_rows;
  innovations.set_size(w_tilde.n_rows, w_tilde.n_cols);
  for (int i = 0; i < n; ++i) {
    innovations.row(i) = compute_innovation_row_cpp(i, w_tilde, precomp, values);
  }
}

double compute_quad_from_innovations_cpp(const arma::mat& innovations, const arma::vec& f) {
  double quad = 0.0;
  for (arma::uword i = 0; i < innovations.n_rows; ++i) {
    quad += arma::dot(innovations.row(i), innovations.row(i)) / f(i);
  }
  return quad;
}

double matern_cov_scalar(double d, double nu, double rho, double sigma_i, double sigma_j) {
  if (d == 0.0) return sigma_i * sigma_j;
  if (rho <= 0.0 || nu <= 0.0 || sigma_i < 0.0 || sigma_j < 0.0) {
    Rcpp::stop("Invalid Matern covariance parameters.");
  }
  
  double arg = sqrt(2.0 * nu) * d / rho;
  double term1 = sigma_i * sigma_j * (std::pow(2.0, 1.0 - nu) / R::gammafn(nu));
  double bessel_k_val = 0.0;
  try {
    bessel_k_val = boost::math::cyl_bessel_k(nu, arg);
  } catch (...) {
    Rcpp::stop("Invalid Matern Bessel evaluation at nonzero distance.");
  }
  double term2 = std::pow(arg, nu) * bessel_k_val;
  double covariance = term1 * term2;

  if (!std::isfinite(covariance)) {
    Rcpp::stop("Invalid Matern covariance evaluation at nonzero distance.");
  }
  return covariance;
}

List precompute_NNGP_cpp(const arma::mat& coords, const arma::umat& NN_ind) {
  int n = coords.n_rows;
  List D_list(n);
  List neighbors_list(n);
  
  arma::umat padded_NN_ind = NN_ind;
  if (padded_NN_ind.n_rows < static_cast<arma::uword>(n)) {
    arma::umat zeros(static_cast<arma::uword>(n) - padded_NN_ind.n_rows, padded_NN_ind.n_cols, arma::fill::zeros);
    padded_NN_ind = arma::join_cols(zeros, padded_NN_ind);
  }

  arma::uvec row_offsets(n + 1, arma::fill::zeros);
  std::vector<std::vector<unsigned int>> dependents_vec(n);
  int nnz_est = 0;
  
  for (int i = 0; i < n; ++i) {
    arma::uvec row_ind = padded_NN_ind.row(i).t();
    arma::uvec neighbors_1based = row_ind(arma::find(row_ind > 0));
    int m = neighbors_1based.n_elem;
    row_offsets(i) = nnz_est;
    nnz_est += m;
    dependents_vec[i].push_back(i);
    
    if (m == 0) {
      D_list[i] = R_NilValue;
      neighbors_list[i] = arma::uvec();
    } else {
      arma::uvec neighbors = neighbors_1based - 1; 
      neighbors_list[i] = neighbors;
      for (int j = 0; j < m; ++j) {
        dependents_vec[neighbors(j)].push_back(i);
      }
      
      arma::mat coords_sub(m + 1, coords.n_cols);
      coords_sub.row(0) = coords.row(i);
      for (int j = 0; j < m; ++j) {
        coords_sub.row(j + 1) = coords.row(neighbors(j));
      }
      
      arma::mat D(m + 1, m + 1, arma::fill::zeros);
      for (int r = 0; r < m + 1; ++r) {
        for (int c = r + 1; c < m + 1; ++c) {
          double dist = arma::norm(coords_sub.row(r) - coords_sub.row(c), 2);
          D(r, c) = dist;
          D(c, r) = dist; 
        }
      }
      D_list[i] = D;
    }
  }
  row_offsets(n) = nnz_est;

  arma::umat locations(2, nnz_est);
  int counter = 0;
  for (int i = 0; i < n; ++i) {
    arma::uvec row_ind = padded_NN_ind.row(i).t();
    arma::uvec neighbors_1based = row_ind(arma::find(row_ind > 0));
    arma::uvec neighbors = neighbors_1based - 1;
    for (arma::uword j = 0; j < neighbors.n_elem; ++j) {
      locations(0, counter) = i;
      locations(1, counter) = neighbors(j);
      counter++;
    }
  }

  List dependents(n);
  for (int i = 0; i < n; ++i) {
    dependents[i] = arma::conv_to<arma::uvec>::from(dependents_vec[i]);
  }
  
  return List::create(Named("NN_ind") = padded_NN_ind,
                      Named("D_list") = D_list,
                      Named("neighbors") = neighbors_list,
                      Named("row_offsets") = row_offsets,
                      Named("locations") = locations,
                      Named("dependents") = dependents,
                      Named("n") = n);
}

void compute_nngp_row_factor_cpp(
    const List& neighbors_list,
    const List& D_list,
    const arma::vec& sigma_vec,
    arma::uword row_idx,
    double rho,
    double nu,
    double& F_out,
    arma::vec& B_out
) {
  arma::uvec neighbors = as<arma::uvec>(neighbors_list[row_idx]);
  int m = neighbors.n_elem;

  if (m == 0) {
    double sigma_i = sigma_vec(row_idx);
    F_out = sigma_i * sigma_i;
    B_out.reset();
    return;
  }

  arma::mat D = as<arma::mat>(D_list[row_idx]);

  double sigma_i = sigma_vec(row_idx);
  arma::vec c_i(m);
  for (int j = 0; j < m; ++j) {
    c_i(j) = matern_cov_scalar(D(0, j + 1), nu, rho, sigma_i, sigma_vec(neighbors(j)));
  }

  arma::mat C_N(m, m);
  for (int r = 0; r < m; ++r) {
    for (int c = r; c < m; ++c) {
      if (r == c) {
        C_N(r, c) = sigma_vec(neighbors(r)) * sigma_vec(neighbors(r));
      } else {
        double cov_val = matern_cov_scalar(
          D(r + 1, c + 1), nu, rho, sigma_vec(neighbors(r)), sigma_vec(neighbors(c))
        );
        C_N(r, c) = cov_val;
        C_N(c, r) = cov_val;
      }
    }
  }

  double marginal_var = sigma_i * sigma_i;
  double base_jitter = std::max(1e-8 * marginal_var, 1e-12);
  double minimum_conditional_variance = 0.5 * base_jitter;
  double F_i = arma::datum::nan;
  arma::vec B_candidate;

  for (int attempt = 0; attempt < 7; ++attempt) {
    arma::mat C_N_attempt = C_N;
    if (attempt > 0) {
      C_N_attempt.diag() += base_jitter * std::pow(10.0, attempt - 1);
    }

    bool success = arma::solve(
      B_candidate,
      C_N_attempt,
      c_i,
      arma::solve_opts::fast + arma::solve_opts::no_approx
    );
    if (!success || !B_candidate.is_finite()) {
      continue;
    }

    F_i = marginal_var - arma::dot(c_i, B_candidate);
    if (std::isfinite(F_i) && F_i >= minimum_conditional_variance) {
      B_out = B_candidate;
      F_out = F_i;
      return;
    }
  }

  Rcpp::stop(
    "Invalid NNGP conditional variance at row %d: F_i=%g, rho=%g, sigma2=%g.",
    static_cast<int>(row_idx) + 1,
    F_i,
    rho,
    marginal_var
  );
}

arma::sp_mat build_Q_from_nngp_state_cpp(
    const List& precomp,
    const arma::vec& values,
    const arma::vec& F_vec
) {
  int n = as<int>(precomp["n"]);
  arma::umat locations = as<arma::umat>(precomp["locations"]);
  arma::sp_mat A(locations, values, n, n);
  arma::sp_mat I_n(n, n);
  I_n.eye();
  arma::sp_mat ImA = I_n - A;

  arma::sp_mat Dinv(n, n);
  for (int i = 0; i < n; ++i) Dinv(i, i) = 1.0 / F_vec(i);

  arma::sp_mat Q = ImA.t() * Dinv * ImA;
  return 0.5 * (Q + Q.t());
}

arma::sp_mat compute_Q_fast_cpp(
    const List& precomp,
    const arma::vec& sigma2_vec,
    double rho,
    double nu,
    arma::vec& F_vec,
    arma::vec& values
) {
  int n = as<int>(precomp["n"]);
  List neighbors_list = as<List>(precomp["neighbors"]);
  List D_list = as<List>(precomp["D_list"]);
  arma::uvec row_offsets = as<arma::uvec>(precomp["row_offsets"]);
  arma::vec sigma_sd = arma::sqrt(sigma2_vec);

  F_vec.set_size(n);
  values.set_size(row_offsets(n));

  for (int i = 0; i < n; ++i) {
    arma::vec B_i;
    compute_nngp_row_factor_cpp(neighbors_list, D_list, sigma_sd, i, rho, nu, F_vec(i), B_i);

    arma::uword start = row_offsets(i);
    if (B_i.n_elem > 0) {
      values.subvec(start, start + B_i.n_elem - 1) = B_i;
    }
  }

  return build_Q_from_nngp_state_cpp(precomp, values, F_vec);
}

void update_Q_fast_single_row_cpp(
    const List& precomp,
    const arma::vec& sigma2_vec,
    arma::uword update_idx,
    double rho,
    double nu,
    arma::vec& F_vec,
    arma::vec& values
) {
  List neighbors_list = as<List>(precomp["neighbors"]);
  List D_list = as<List>(precomp["D_list"]);
  arma::uvec row_offsets = as<arma::uvec>(precomp["row_offsets"]);
  List dependents = as<List>(precomp["dependents"]);
  arma::uvec affected_rows = as<arma::uvec>(dependents[update_idx]);
  arma::vec sigma_sd = arma::sqrt(sigma2_vec);

  for (arma::uword idx = 0; idx < affected_rows.n_elem; ++idx) {
    arma::uword row = affected_rows(idx);
    arma::vec B_i;
    compute_nngp_row_factor_cpp(neighbors_list, D_list, sigma_sd, row, rho, nu, F_vec(row), B_i);

    arma::uword start = row_offsets(row);
    if (B_i.n_elem > 0) {
      values.subvec(start, start + B_i.n_elem - 1) = B_i;
    }
  }
}

arma::sp_mat compute_global_nngp_state_cpp(
    const List& precomp,
    int n,
    double sigma2_u,
    double rho,
    double nu,
    arma::vec& F_u,
    arma::vec& values_u
) {
  arma::vec sigma2_u_vec(n, arma::fill::value(sigma2_u));
  return compute_Q_fast_cpp(precomp, sigma2_u_vec, rho, nu, F_u, values_u);
}

arma::sp_mat build_Tauinv_u_cpp(
    const List& precomp,
    int n,
    double sigma2_u,
    double rho,
    double nu,
    arma::vec& F_u,
    arma::vec& values_u
) {
  return compute_global_nngp_state_cpp(precomp, n, sigma2_u, rho, nu, F_u, values_u);
}

arma::sp_mat build_exact_gp_state_cpp(
    const arma::mat& locations,
    double sigma2_u,
    double rho,
    double nu,
    arma::vec& F_u,
    arma::vec& values_u
) {
  arma::uword n = locations.n_rows;
  arma::mat covariance(n, n, arma::fill::zeros);
  double marginal_sd = std::sqrt(sigma2_u);
  for (arma::uword row = 0; row < n; ++row) {
    covariance(row, row) = sigma2_u;
    for (arma::uword column = row + 1; column < n; ++column) {
      double distance = arma::norm(
        locations.row(row) - locations.row(column),
        2
      );
      double covariance_value = matern_cov_scalar(
        distance, nu, rho, marginal_sd, marginal_sd
      );
      covariance(row, column) = covariance_value;
      covariance(column, row) = covariance_value;
    }
  }

  arma::mat covariance_cholesky;
  arma::mat stabilized_covariance;
  bool factorized = false;
  double base_jitter = std::max(1e-10 * sigma2_u, 1e-12);
  for (int attempt = 0; attempt < 7; ++attempt) {
    stabilized_covariance = covariance;
    if (attempt > 0) {
      stabilized_covariance.diag() +=
        base_jitter * std::pow(10.0, attempt - 1);
    }
    factorized = arma::chol(
      covariance_cholesky,
      stabilized_covariance,
      "lower"
    );
    if (factorized) {
      break;
    }
  }
  if (!factorized) {
    Rcpp::stop("Exact GP covariance is not positive definite.");
  }

  arma::mat precision;
  if (!arma::inv_sympd(precision, stabilized_covariance) ||
      !precision.is_finite()) {
    Rcpp::stop("Unable to invert the exact GP covariance.");
  }
  F_u = arma::square(covariance_cholesky.diag());
  values_u.reset();
  arma::mat symmetric_precision = 0.5 * (precision + precision.t());
  return arma::sp_mat(symmetric_precision);
}

List precompute_Q_ar1_diag_structure_cpp(int n, int T) {
  arma::uword total_nnz = static_cast<arma::uword>(n) * static_cast<arma::uword>(T) +
    2u * static_cast<arma::uword>(n) * static_cast<arma::uword>(std::max(T - 1, 0));
  arma::umat q_locs(2, total_nnz);
  arma::uvec block_type(total_nnz, arma::fill::zeros);
  arma::uvec time_index(total_nnz, arma::fill::zeros);
  arma::uvec site_index(total_nnz, arma::fill::zeros);
  arma::uvec diag_positions(static_cast<arma::uword>(n * T), arma::fill::zeros);

  arma::uword counter = 0;
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < n; ++i) {
      arma::uword idx = static_cast<arma::uword>(t * n + i);
      q_locs(0, counter) = idx;
      q_locs(1, counter) = idx;
      block_type(counter) = 0;
      time_index(counter) = static_cast<arma::uword>(t);
      site_index(counter) = static_cast<arma::uword>(i);
      diag_positions(idx) = counter;
      ++counter;
    }
  }

  for (int t = 0; t < T - 1; ++t) {
    for (int i = 0; i < n; ++i) {
      arma::uword idx_curr = static_cast<arma::uword>(t * n + i);
      arma::uword idx_next = static_cast<arma::uword>((t + 1) * n + i);

      q_locs(0, counter) = idx_curr;
      q_locs(1, counter) = idx_next;
      block_type(counter) = 1;
      time_index(counter) = static_cast<arma::uword>(t);
      site_index(counter) = static_cast<arma::uword>(i);
      ++counter;

      q_locs(0, counter) = idx_next;
      q_locs(1, counter) = idx_curr;
      block_type(counter) = 2;
      time_index(counter) = static_cast<arma::uword>(t);
      site_index(counter) = static_cast<arma::uword>(i);
      ++counter;
    }
  }

  return List::create(
    Named("q_locs") = q_locs,
    Named("block_type") = block_type,
    Named("time_index") = time_index,
    Named("site_index") = site_index,
    Named("diag_positions") = diag_positions,
    Named("n") = n,
    Named("T") = T
  );
}

arma::sp_mat build_Q_ar1_diag_fast_cpp(
    const List& q_precomp,
    const arma::vec& phi,
    const arma::vec& tau2,
    arma::vec& q_vals
) {
  arma::umat q_locs = as<arma::umat>(q_precomp["q_locs"]);
  arma::uvec block_type = as<arma::uvec>(q_precomp["block_type"]);
  arma::uvec time_index = as<arma::uvec>(q_precomp["time_index"]);
  arma::uvec site_index = as<arma::uvec>(q_precomp["site_index"]);
  int n = as<int>(q_precomp["n"]);
  int T = as<int>(q_precomp["T"]);

  q_vals.set_size(q_locs.n_cols);

  for (arma::uword idx = 0; idx < q_locs.n_cols; ++idx) {
    arma::uword site = site_index(idx);
    arma::uword time = time_index(idx);
    double phi_i = phi(site);
    double inv_tau2_i = 1.0 / tau2(site);

    if (block_type(idx) == 0) {
      double val = 0.0;
      if (T == 1) {
        val = (1.0 - phi_i * phi_i) * inv_tau2_i;
      } else if (time == 0 || time == static_cast<arma::uword>(T - 1)) {
        val = inv_tau2_i;
      } else {
        val = (1.0 + phi_i * phi_i) * inv_tau2_i;
      }
      q_vals(idx) = val;
    } else {
      q_vals(idx) = -phi_i * inv_tau2_i;
    }
  }

  return arma::sp_mat(q_locs, q_vals, n * T, n * T);
}

arma::vec build_w_prior_canonical_cpp(
    const arma::vec& u,
    const arma::vec& phi,
    const arma::vec& tau2,
    int n,
    int T
) {
  arma::vec rhs(n * T, arma::fill::zeros);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);

  for (int t = 1; t < T; ++t) {
    for (int i = 0; i < n; ++i) {
      double scaled_u = u_mat(i, t) / tau2(i);
      rhs(t * n + i) += scaled_u;
      rhs((t - 1) * n + i) -= phi(i) * scaled_u;
    }
  }

  return rhs;
}

struct TransitionSufficientStatistics {
  arma::vec constant;
  arma::vec linear;
  arma::vec quadratic;
  int n_times;
};

TransitionSufficientStatistics build_transition_statistics_cpp(
    const arma::vec& w,
    const arma::vec& u,
    int n,
    int T
) {
  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);
  TransitionSufficientStatistics statistics{
    arma::vec(n, arma::fill::zeros),
    arma::vec(n, arma::fill::zeros),
    arma::vec(n, arma::fill::zeros),
    T
  };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int site = 0; site < n; ++site) {
    double w_initial = w_mat(site, 0);
    double constant = w_initial * w_initial;
    double linear = 0.0;
    double quadratic = -w_initial * w_initial;
    for (int time = 1; time < T; ++time) {
      double adjusted_state = w_mat(site, time) - u_mat(site, time);
      double previous_state = w_mat(site, time - 1);
      constant += adjusted_state * adjusted_state;
      linear -= 2.0 * adjusted_state * previous_state;
      quadratic += previous_state * previous_state;
    }
    statistics.constant(site) = constant;
    statistics.linear(site) = linear;
    statistics.quadratic(site) = quadratic;
  }

  return statistics;
}

double transition_sse_from_statistics_cpp(
    const TransitionSufficientStatistics& statistics,
    arma::uword site,
    double phi_candidate
) {
  return statistics.constant(site) +
    phi_candidate * statistics.linear(site) +
    phi_candidate * phi_candidate * statistics.quadratic(site);
}

double evaluate_candidate_transition_loglik_cached_cpp(
    const TransitionSufficientStatistics& statistics,
    arma::uword site,
    double phi_candidate,
    double tau2_candidate
) {
  if (phi_candidate <= 0.0 || phi_candidate >= 1.0 ||
      tau2_candidate <= 0.0) {
    return -arma::datum::inf;
  }
  double one_minus_phi2 = 1.0 - phi_candidate * phi_candidate;
  if (one_minus_phi2 <= 0.0) {
    return -arma::datum::inf;
  }
  double sse = transition_sse_from_statistics_cpp(
    statistics,
    site,
    phi_candidate
  );
  return 0.5 * std::log(one_minus_phi2) -
    0.5 * static_cast<double>(statistics.n_times) * std::log(tau2_candidate) -
    0.5 * sse / tau2_candidate;
}

double compute_transition_loglik_cached_cpp(
    const TransitionSufficientStatistics& statistics,
    const arma::vec& phi,
    const arma::vec& tau2
) {
  double loglik = 0.0;
  for (arma::uword site = 0; site < phi.n_elem; ++site) {
    loglik += evaluate_candidate_transition_loglik_cached_cpp(
      statistics,
      site,
      phi(site),
      tau2(site)
    );
  }
  return loglik;
}

double compute_cluster_transition_loglik_cached_cpp(
    const TransitionSufficientStatistics& statistics,
    const arma::uvec& members,
    double phi_candidate,
    double tau2_candidate
) {
  double loglik = 0.0;
  for (arma::uword index = 0; index < members.n_elem; ++index) {
    loglik += evaluate_candidate_transition_loglik_cached_cpp(
      statistics,
      members(index),
      phi_candidate,
      tau2_candidate
    );
  }
  return loglik;
}

double compute_cluster_tau2_sse_cached_cpp(
    const TransitionSufficientStatistics& statistics,
    const arma::uvec& members,
    double phi_candidate
) {
  double sse = 0.0;
  for (arma::uword index = 0; index < members.n_elem; ++index) {
    sse += transition_sse_from_statistics_cpp(
      statistics,
      members(index),
      phi_candidate
    );
  }
  return sse;
}

double evaluate_candidate_local_transition_loglik(
    const arma::vec& w,
    const arma::vec& u,
    arma::uword i,
    double phi_candidate,
    double tau2_candidate,
    int n,
    int T
) {
  if (phi_candidate <= 0.0 || phi_candidate >= 1.0 || tau2_candidate <= 0.0) {
    return -arma::datum::inf;
  }

  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);

  double one_minus_phi2 = 1.0 - phi_candidate * phi_candidate;
  if (one_minus_phi2 <= 0.0) {
    return -arma::datum::inf;
  }

  double loglik =
    0.5 * std::log(one_minus_phi2) -
    0.5 * std::log(tau2_candidate) -
    0.5 * one_minus_phi2 * std::pow(w_mat(i, 0), 2.0) / tau2_candidate;

  for (int t = 1; t < T; ++t) {
    double resid = w_mat(i, t) - phi_candidate * w_mat(i, t - 1) - u_mat(i, t);
    loglik += -0.5 * std::log(tau2_candidate) - 0.5 * resid * resid / tau2_candidate;
  }

  return loglik;
}

double compute_transition_loglik_w_given_u_cpp(
    const arma::vec& w,
    const arma::vec& u,
    const arma::vec& phi,
    const arma::vec& tau2,
    int n,
    int T
) {
  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);

  double loglik = 0.0;
  for (int i = 0; i < n; ++i) {
    double phi_i = phi(i);
    double tau2_i = tau2(i);
    if (phi_i <= 0.0 || phi_i >= 1.0 || tau2_i <= 0.0) {
      return -arma::datum::inf;
    }

    double one_minus_phi2 = 1.0 - phi_i * phi_i;
    if (one_minus_phi2 <= 0.0) {
      return -arma::datum::inf;
    }

    double w1 = w_mat(i, 0);
    loglik +=
      0.5 * std::log(one_minus_phi2) -
      0.5 * std::log(tau2_i) -
      0.5 * one_minus_phi2 * w1 * w1 / tau2_i;

    for (int t = 1; t < T; ++t) {
      double resid = w_mat(i, t) - phi_i * w_mat(i, t - 1) - u_mat(i, t);
      loglik += -0.5 * std::log(tau2_i) - 0.5 * resid * resid / tau2_i;
    }
  }

  return loglik;
}

double compute_cluster_tau2_sse_cpp(
    const arma::vec& w,
    const arma::vec& u,
    const arma::uvec& members,
    double phi_k,
    int n,
    int T
) {
  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);

  double sse = 0.0;
  double one_minus_phi2 = 1.0 - phi_k * phi_k;
  for (arma::uword idx = 0; idx < members.n_elem; ++idx) {
    arma::uword i = members(idx);
    sse += one_minus_phi2 * std::pow(w_mat(i, 0), 2.0);
    for (int t = 1; t < T; ++t) {
      double resid = w_mat(i, t) - phi_k * w_mat(i, t - 1) - u_mat(i, t);
      sse += resid * resid;
    }
  }

  return sse;
}

double compute_nngp_loglik_u_cpp(
    const arma::vec& u,
    const arma::vec& F_u,
    const arma::vec& values_u,
    const List& precomp,
    int n,
    int T
) {
  if (T <= 1) {
    return 0.0;
  }

  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);
  arma::mat u_active = u_mat.cols(1, T - 1);
  arma::mat innovations;
  build_innovation_cache_cpp(u_active, precomp, values_u, innovations);
  double quad = compute_quad_from_innovations_cpp(innovations, F_u);
  double log_det = -0.5 * static_cast<double>(T - 1) * arma::sum(arma::log(F_u));
  return log_det - 0.5 * quad;
}

double compute_sigma2_u_sse_cpp(
    const arma::vec& u,
    const arma::vec& F_u,
    const arma::vec& values_u,
    const List& precomp,
    int n,
    int T,
    double sigma2_u
) {
  if (T <= 1) {
    return 0.0;
  }

  arma::mat u_mat(const_cast<double*>(u.memptr()), n, T, false, false);
  arma::mat u_active = u_mat.cols(1, T - 1);
  arma::mat innovations;
  build_innovation_cache_cpp(u_active, precomp, values_u, innovations);
  double quad = compute_quad_from_innovations_cpp(innovations, F_u);
  return sigma2_u * quad;
}

arma::sp_mat build_u_precision_cpp(const arma::sp_mat& Tauinv_u, const arma::vec& tau2) {
  arma::sp_mat precision_u = Tauinv_u;
  for (arma::uword i = 0; i < tau2.n_elem; ++i) {
    precision_u(i, i) += 1.0 / tau2(i);
  }
  return precision_u;
}

double collapsed_spatial_log_target_from_state_cpp(
    double log_rho,
    double log_sigma2_u,
    const arma::vec& w,
    const arma::vec& phi,
    const arma::vec& tau2,
    const arma::sp_mat& Tauinv_u,
    const arma::vec& F_u,
    CholmodSampler& sampler,
    int n,
    int T,
    double a_rho,
    double b_rho,
    double a_sigma2_u,
    double b_sigma2_u
) {
  double rho = std::exp(log_rho);
  double sigma2_u = std::exp(log_sigma2_u);
  if (!std::isfinite(rho) || rho <= 0.0 ||
      !std::isfinite(sigma2_u) || sigma2_u <= 0.0 ||
      !F_u.is_finite() || arma::any(F_u <= 0.0) ||
      !tau2.is_finite() || arma::any(tau2 <= 0.0)) {
    return -arma::datum::inf;
  }

  double log_target =
    a_rho * log_rho - b_rho * rho -
    a_sigma2_u * log_sigma2_u - b_sigma2_u / sigma2_u;
  if (T <= 1) {
    return log_target;
  }

  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::mat state_residuals(n, T - 1);
  for (int t = 1; t < T; ++t) {
    state_residuals.col(t - 1) = w_mat.col(t) - phi % w_mat.col(t - 1);
  }

  arma::vec inv_tau2 = 1.0 / tau2;
  arma::mat rhs = state_residuals.each_col() % inv_tau2;
  arma::sp_mat precision = build_u_precision_cpp(Tauinv_u, tau2);
  sampler.factorize_numeric(precision);
  double log_det_precision = sampler.log_determinant_factorized();
  arma::mat solved_rhs = sampler.solve_means_factorized(rhs);

  double direct_quadratic = arma::accu(state_residuals % rhs);
  double correction = arma::accu(rhs % solved_rhs);
  double quadratic = direct_quadratic - correction;
  double tolerance = 1e-8 * std::max(1.0, std::abs(direct_quadratic));
  if (quadratic < -tolerance || !std::isfinite(quadratic)) {
    return -arma::datum::inf;
  }
  quadratic = std::max(0.0, quadratic);

  double log_det_marginal =
    arma::sum(arma::log(F_u)) +
    arma::sum(arma::log(tau2)) +
    log_det_precision;
  log_target +=
    -0.5 * static_cast<double>(T - 1) * log_det_marginal -
    0.5 * quadratic;
  return log_target;
}

// [[Rcpp::export]]
double collapsed_spatial_log_target_cpp(
    const arma::vec& transformed_spatial,
    const arma::vec& w,
    const arma::vec& phi,
    const arma::vec& tau2,
    const arma::mat& locations,
    const arma::umat& NN_ind,
    int T,
    double nu,
    double a_kappa,
    double b_kappa,
    double a_sigma2_u,
    double b_sigma2_u
) {
  int n = static_cast<int>(phi.n_elem);
  if (transformed_spatial.n_elem != 2 || tau2.n_elem != phi.n_elem ||
      w.n_elem != static_cast<arma::uword>(n * T) ||
      locations.n_rows != static_cast<arma::uword>(n)) {
    return -arma::datum::inf;
  }

  double rho = std::exp(transformed_spatial(0));
  double sigma2_u = std::exp(transformed_spatial(1));
  List precomp = precompute_NNGP_cpp(locations, NN_ind);
  arma::vec F_u;
  arma::vec values_u;
  arma::sp_mat Tauinv_u = build_Tauinv_u_cpp(
    precomp, n, sigma2_u, rho, nu, F_u, values_u
  );
  CholmodSampler sampler;
  sampler.analyze_pattern(build_u_precision_cpp(Tauinv_u, tau2));
  return collapsed_spatial_log_target_from_state_cpp(
    transformed_spatial(0), transformed_spatial(1), w, phi, tau2,
    Tauinv_u, F_u, sampler, n, T, a_kappa, b_kappa,
    a_sigma2_u, b_sigma2_u
  );
}

// [[Rcpp::export]]
List global_nngp_state_test_cpp(
    const arma::mat& locations,
    const arma::umat& NN_ind,
    double sigma2_u,
    double kappa,
    double nu
) {
  int n = static_cast<int>(locations.n_rows);
  List precomp = precompute_NNGP_cpp(locations, NN_ind);
  arma::vec F_u;
  arma::vec values_u;
  arma::sp_mat precision = build_Tauinv_u_cpp(
    precomp, n, sigma2_u, kappa, nu, F_u, values_u
  );
  return List::create(
    Named("precision") = precision,
    Named("F") = F_u
  );
}

arma::vec update_u_cpp(
    CholmodSampler& sampler,
    const arma::vec& w,
    const arma::vec& phi,
    const arma::vec& tau2,
    const arma::sp_mat& Tauinv_u,
    int n,
    int T
) {
  arma::vec u_new(n * T, arma::fill::zeros);
  if (T <= 1) {
    return u_new;
  }

#ifndef _OPENMP
  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::sp_mat precision_u = build_u_precision_cpp(Tauinv_u, tau2);
  sampler.factorize_numeric(precision_u);

  for (int t = 1; t < T; ++t) {
    arma::vec rhs_t(n, arma::fill::zeros);
    for (int i = 0; i < n; ++i) {
      double a_ti = w_mat(i, t) - phi(i) * w_mat(i, t - 1);
      rhs_t(i) = a_ti / tau2(i);
    }
    u_new.subvec(t * n, (t + 1) * n - 1) = sampler.sample_factorized(rhs_t);
  }

  return u_new;
#else
  if (omp_get_max_threads() <= 1) {
    arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
    arma::sp_mat precision_u = build_u_precision_cpp(Tauinv_u, tau2);
    sampler.factorize_numeric(precision_u);

    for (int t = 1; t < T; ++t) {
      arma::vec rhs_t = (1.0 / tau2) %
        (w_mat.col(t) - phi % w_mat.col(t - 1));
      u_new.subvec(t * n, (t + 1) * n - 1) =
        sampler.sample_factorized(rhs_t);
    }
    return u_new;
  }

  arma::mat w_mat(const_cast<double*>(w.memptr()), n, T, false, false);
  arma::sp_mat precision_u = build_u_precision_cpp(Tauinv_u, tau2);
  sampler.factorize_numeric(precision_u);
  arma::vec inv_tau2 = 1.0 / tau2;
  arma::mat right_hand_sides(n, T - 1, arma::fill::zeros);

#pragma omp parallel for schedule(static)
  for (int t = 1; t < T; ++t) {
    right_hand_sides.col(t - 1) = inv_tau2 %
      (w_mat.col(t) - phi % w_mat.col(t - 1));
  }

  arma::mat standard_normals = arma::randn<arma::mat>(n, T - 1);
  arma::mat active_draws = sampler.sample_factorized_columns_parallel(
    right_hand_sides,
    standard_normals
  );
  arma::mat u_new_mat(u_new.memptr(), n, T, false, false);
  u_new_mat.cols(1, T - 1) = active_draws;

  return u_new;
#endif
}

// [[Rcpp::export]]
double validate_parallel_sparse_solves_cpp(
    const arma::sp_mat& precision,
    const arma::mat& right_hand_sides
) {
  CholmodSampler sampler;
  sampler.analyze_pattern(precision);
  sampler.factorize_numeric(precision);
  arma::mat reference = sampler.solve_means_factorized(right_hand_sides);
  arma::mat candidate = sampler.sample_factorized_columns_parallel(
    right_hand_sides,
    arma::mat(right_hand_sides.n_rows, right_hand_sides.n_cols,
              arma::fill::zeros)
  );
  return arma::abs(reference - candidate).max();
}

// CHECKED
List scaleandperiods(const arma::mat& data, bool scale) {
  // Get the dimensions of the data
  //int n = data.n_rows;
  int m = data.n_cols;
  
  // Initialize maxima and minima matrices
  arma::vec maxima(m);
  arma::vec minima(m);
  
  // Calculate maxima and minima for each column (time series)
  for (int i = 0; i < m; i++) {
    maxima(i) = data.col(i).max();
    minima(i) = data.col(i).min();
  }
  
  // Identify columns with constant time series
  arma::uvec cts = arma::find(maxima == minima);
  
  // If there are constant time series, remove them
  arma::mat mydata = data;
  if (!cts.is_empty()) {
    Rcpp::Rcout << "Removing series ";
    for (unsigned int i = 0; i < cts.n_elem; i++) {
      Rcpp::Rcout << (cts(i) + 1) << " ";
    }
    Rcpp::Rcout << "since they are constant." << std::endl;
    
    // Remove constant columns
    mydata.shed_cols(cts);
    m = mydata.n_cols;
    
    // Recompute maxima and minima for the remaining data
    maxima.resize(m);
    minima.resize(m);
    for (int i = 0; i < m; i++) {
      maxima(i) = mydata.col(i).max();
      minima(i) = mydata.col(i).min();
    }
  }
  
  // Scale data if required
  if (scale) {
    for (int j = 0; j < m; j++) {
      double range = maxima(j) - minima(j);
      if (range != 0) {
        mydata.col(j) = 1 + (1.0 / range) * (mydata.col(j) - maxima(j));
      }
    }
  }
  
  // Output the list
  return List::create(
    Named("mydata") = mydata,
    Named("cts") = cts
  );
}

////////////////////////////////////////////////////////////////////////////////
// Fixed sparse-structure builders
////////////////////////////////////////////////////////////////////////////////

List precompute_Q_ar1_nngp_structure_cpp(const arma::sp_mat& Tauinv, int n, int T) {
  arma::umat tau_locs(2, Tauinv.n_nonzero);
  arma::uword nz = 0;
  for (arma::sp_mat::const_iterator it = Tauinv.begin(); it != Tauinv.end(); ++it) {
    tau_locs(0, nz) = it.row();
    tau_locs(1, nz) = it.col();
    ++nz;
  }

  arma::uword tau_nnz = Tauinv.n_nonzero;
  arma::uword total_nnz = static_cast<arma::uword>(T) * tau_nnz +
    static_cast<arma::uword>(std::max(T - 1, 0)) * 2 * tau_nnz;
  arma::umat q_locs(2, total_nnz);
  arma::uvec source_tau_idx(total_nnz);
  arma::uvec block_type(total_nnz);
  arma::uvec diag_positions(n * T, arma::fill::zeros);
  arma::uword q_counter = 0;

  for (int t = 0; t < T; ++t) {
    for (arma::uword idx = 0; idx < tau_nnz; ++idx) {
      arma::uword row = tau_locs(0, idx);
      arma::uword col = tau_locs(1, idx);
      q_locs(0, q_counter) = t * n + row;
      q_locs(1, q_counter) = t * n + col;
      source_tau_idx(q_counter) = idx;
      block_type(q_counter) = 0;
      if (row == col) {
        diag_positions(t * n + row) = q_counter;
      }
      ++q_counter;
    }

    if (t < T - 1) {
      for (arma::uword idx = 0; idx < tau_nnz; ++idx) {
        arma::uword row = tau_locs(0, idx);
        arma::uword col = tau_locs(1, idx);

        q_locs(0, q_counter) = t * n + row;
        q_locs(1, q_counter) = (t + 1) * n + col;
        source_tau_idx(q_counter) = idx;
        block_type(q_counter) = 1;
        ++q_counter;

        q_locs(0, q_counter) = (t + 1) * n + row;
        q_locs(1, q_counter) = t * n + col;
        source_tau_idx(q_counter) = idx;
        block_type(q_counter) = 2;
        ++q_counter;
      }
    }
  }

  return List::create(
    Named("tau_locs") = tau_locs,
    Named("q_locs") = q_locs,
    Named("source_tau_idx") = source_tau_idx,
    Named("block_type") = block_type,
    Named("diag_positions") = diag_positions,
    Named("tau_nnz") = tau_nnz,
    Named("n") = n,
    Named("T") = T
  );
}

arma::sp_mat build_Q_ar1_nngp_fast_cpp(
    const List& q_precomp,
    const arma::sp_mat& Tauinv,
    const arma::vec& phi,
    arma::vec& q_vals
) {
  arma::umat tau_locs = as<arma::umat>(q_precomp["tau_locs"]);
  arma::umat q_locs = as<arma::umat>(q_precomp["q_locs"]);
  arma::uvec source_tau_idx = as<arma::uvec>(q_precomp["source_tau_idx"]);
  arma::uvec block_type = as<arma::uvec>(q_precomp["block_type"]);
  int n = as<int>(q_precomp["n"]);
  int T = as<int>(q_precomp["T"]);

  q_vals.set_size(q_locs.n_cols);
  arma::vec phi_root = arma::sqrt(1.0 - arma::square(phi));

  for (arma::uword idx = 0; idx < q_locs.n_cols; ++idx) {
    arma::uword tau_idx = source_tau_idx(idx);
    arma::uword row = tau_locs(0, tau_idx);
    arma::uword col = tau_locs(1, tau_idx);
    arma::uword block = q_locs(0, idx) / static_cast<arma::uword>(n);
    double tau_val = Tauinv.values[tau_idx];

    if (block_type(idx) == 0) {
      double val = tau_val;
      if (block == 0) {
        val = phi_root(row) * tau_val * phi_root(col);
        if (T > 1) val += phi(row) * tau_val * phi(col);
      } else if (block < static_cast<arma::uword>(T - 1)) {
        val += phi(row) * tau_val * phi(col);
      }
      q_vals(idx) = val;
    } else if (block_type(idx) == 1) {
      q_vals(idx) = -phi(row) * tau_val;
    } else {
      q_vals(idx) = -tau_val * phi(col);
    }
  }

  return arma::sp_mat(q_locs, q_vals, n * T, n * T);
}

arma::sp_mat build_Qw_from_Q_fast_cpp(
    const List& q_precomp,
    const arma::vec& q_vals,
    const arma::vec& sigma_inv_diag,
    arma::vec& q_w_vals
) {
  arma::umat q_locs = as<arma::umat>(q_precomp["q_locs"]);
  arma::uvec diag_positions = as<arma::uvec>(q_precomp["diag_positions"]);
  int n = as<int>(q_precomp["n"]);
  int T = as<int>(q_precomp["T"]);

  q_w_vals = q_vals;
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < n; ++i) {
      q_w_vals(diag_positions(t * n + i)) += sigma_inv_diag(t * n + i);
    }
  }

  return arma::sp_mat(q_locs, q_w_vals, n * T, n * T);
}

////////////////////////////////////////////////////////////////////////////////
// Gibbs update helpers
////////////////////////////////////////////////////////////////////////////////

arma::mat build_site_Q_dense_cpp(double phi_i, double tau2_i, int T) {
  arma::mat q_site(T, T, arma::fill::zeros);
  double inv_tau2_i = 1.0 / tau2_i;

  if (T == 1) {
    q_site(0, 0) = (1.0 - phi_i * phi_i) * inv_tau2_i;
    return q_site;
  }

  q_site.diag().fill((1.0 + phi_i * phi_i) * inv_tau2_i);
  q_site(0, 0) = inv_tau2_i;
  q_site(T - 1, T - 1) = inv_tau2_i;
  q_site.diag(1).fill(-phi_i * inv_tau2_i);
  q_site.diag(-1).fill(-phi_i * inv_tau2_i);
  return q_site;
}

bool use_collapsed_beta_update(std::string beta_update_mode) {
  std::transform(beta_update_mode.begin(), beta_update_mode.end(), beta_update_mode.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (beta_update_mode == "collapsed") {
    return true;
  }
  if (beta_update_mode == "conditional") {
    return false;
  }

  Rcpp::stop("`beta_update_mode` must be either 'collapsed' or 'conditional'.");
}

arma::mat extract_site_design_global_cpp(
    const arma::mat& X,
    int site,
    int n,
    int T
) {
  arma::mat X_site(T, X.n_cols, arma::fill::zeros);
  for (int t = 0; t < T; ++t) {
    X_site.row(t) = X.row(static_cast<arma::uword>(t * n + site));
  }
  return X_site;
}

arma::vec compute_Xbeta_global_cpp(
    const arma::mat& X,
    const arma::vec& beta
) {
  return X * beta;
}

arma::vec xi2_update_global_cpp(
    const arma::vec& beta,
    double a_xi,
    double b_xi
) {
  arma::vec xi2(beta.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    double shape = a_xi + 0.5;
    double scale = b_xi + 0.5 * beta(j) * beta(j);
    xi2(j) = 1.0 / arma::randg<double>(
      arma::distr_param(shape, 1.0 / scale)
    );
  }
  return xi2;
}

arma::vec draw_gaussian_from_precision_cpp(
    const arma::mat& precision,
    const arma::vec& rhs,
    const std::string& update_name
) {
  arma::vec mean;
  if (!arma::solve(
        mean,
        precision,
        rhs,
        arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
      )) {
    Rcpp::stop(update_name + " mean solve failed.");
  }

  arma::mat chol_precision;
  if (!arma::chol(chol_precision, precision)) {
    Rcpp::stop(update_name + " precision is not positive definite.");
  }

  arma::vec noise = arma::solve(
    arma::trimatu(chol_precision),
    arma::randn<arma::vec>(rhs.n_elem)
  );
  return mean + noise;
}

arma::vec beta_update_global_conditional_cpp(
    const arma::mat& X,
    const arma::vec& xi2_inv,
    double sigma_inv,
    const arma::vec& y,
    const arma::vec& w
) {
  arma::mat precision = arma::symmatu(sigma_inv * X.t() * X);
  precision.diag() += xi2_inv;
  arma::vec rhs = sigma_inv * X.t() * (y - w);
  return draw_gaussian_from_precision_cpp(
    precision,
    rhs,
    "Global conditional beta"
  );
}

arma::vec beta_update_global_collapsed_cpp(
    const arma::mat& X,
    int n,
    int T,
    const arma::vec& xi2_inv,
    double sigma_inv,
    const arma::vec& phi,
    const arma::vec& tau2,
    const arma::vec& prior_rhs_w,
    const arma::vec& y
) {
  int p = X.n_cols;
  arma::mat y_mat(const_cast<double*>(y.memptr()), n, T, false, false);
  arma::mat prior_rhs_mat(
    const_cast<double*>(prior_rhs_w.memptr()), n, T, false, false
  );
  int thread_count = 1;
#ifdef _OPENMP
  thread_count = omp_get_max_threads();
#endif
  std::vector<arma::mat> precision_by_thread(
    thread_count,
    arma::mat(p, p, arma::fill::zeros)
  );
  std::vector<arma::vec> rhs_by_thread(
    thread_count,
    arma::vec(p, arma::fill::zeros)
  );
  std::atomic<int> error_code(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int site = 0; site < n; ++site) {
    int thread_index = 0;
#ifdef _OPENMP
    thread_index = omp_get_thread_num();
#endif
    arma::mat X_site = extract_site_design_global_cpp(X, site, n, T);
    arma::vec y_site = y_mat.row(site).t();
    arma::vec h_site = prior_rhs_mat.row(site).t();

    arma::mat solve_rhs(T, p + 2, arma::fill::zeros);
    solve_rhs.col(0) = sigma_inv * y_site;
    solve_rhs.col(1) = h_site;
    solve_rhs.cols(2, p + 1) = sigma_inv * X_site;

    arma::vec diagonal(T, arma::fill::zeros);
    arma::vec cholesky_diagonal(T, arma::fill::zeros);
    arma::vec cholesky_subdiagonal(std::max(0, T - 1), arma::fill::zeros);
    double inv_tau2 = 1.0 / tau2(site);
    double off_diagonal = -phi(site) * inv_tau2;
    if (T == 1) {
      diagonal(0) = (1.0 - phi(site) * phi(site)) * inv_tau2 + sigma_inv;
    } else {
      diagonal.fill((1.0 + phi(site) * phi(site)) * inv_tau2 + sigma_inv);
      diagonal(0) = inv_tau2 + sigma_inv;
      diagonal(T - 1) = inv_tau2 + sigma_inv;
    }

    bool factor_ok = std::isfinite(diagonal(0)) && diagonal(0) > 0.0;
    if (factor_ok) {
      cholesky_diagonal(0) = std::sqrt(diagonal(0));
      for (int t = 1; t < T; ++t) {
        cholesky_subdiagonal(t - 1) =
          off_diagonal / cholesky_diagonal(t - 1);
        double pivot = diagonal(t) -
          cholesky_subdiagonal(t - 1) * cholesky_subdiagonal(t - 1);
        if (!std::isfinite(pivot) || pivot <= 0.0) {
          factor_ok = false;
          break;
        }
        cholesky_diagonal(t) = std::sqrt(pivot);
      }
    }
    if (!factor_ok) {
      error_code.store(1, std::memory_order_relaxed);
      continue;
    }

    arma::mat solved(T, p + 2, arma::fill::zeros);
    solved.row(0) = solve_rhs.row(0) / cholesky_diagonal(0);
    for (int t = 1; t < T; ++t) {
      solved.row(t) =
        (solve_rhs.row(t) - cholesky_subdiagonal(t - 1) * solved.row(t - 1)) /
        cholesky_diagonal(t);
    }
    solved.row(T - 1) /= cholesky_diagonal(T - 1);
    for (int t = T - 2; t >= 0; --t) {
      solved.row(t) =
        (solved.row(t) - cholesky_subdiagonal(t) * solved.row(t + 1)) /
        cholesky_diagonal(t);
    }

    arma::vec Sinv_y = sigma_inv * y_site - sigma_inv * solved.col(0);
    arma::vec Sinv_m = sigma_inv * solved.col(1);
    arma::mat Sinv_X = sigma_inv * X_site -
      sigma_inv * solved.cols(2, p + 1);

    precision_by_thread[thread_index] += X_site.t() * Sinv_X;
    rhs_by_thread[thread_index] += X_site.t() * (Sinv_y - Sinv_m);
  }

  if (error_code.load(std::memory_order_relaxed) != 0) {
    Rcpp::stop("Global collapsed beta tridiagonal factorization failed.");
  }

  arma::mat precision = arma::diagmat(xi2_inv);
  arma::vec rhs(p, arma::fill::zeros);
  for (int thread_index = 0; thread_index < thread_count; ++thread_index) {
    precision += precision_by_thread[thread_index];
    rhs += rhs_by_thread[thread_index];
  }

  precision = arma::symmatu(precision);
  return draw_gaussian_from_precision_cpp(
    precision,
    rhs,
    "Global collapsed beta"
  );
}

arma::vec beta_update_global_dispatch_cpp(
    bool use_collapsed,
    const arma::mat& X,
    int n,
    int T,
    const arma::vec& xi2_inv,
    double sigma_inv,
    const arma::vec& phi,
    const arma::vec& tau2,
    const arma::vec& prior_rhs_w,
    const arma::vec& y,
    const arma::vec& w
) {
  if (use_collapsed) {
    return beta_update_global_collapsed_cpp(
      X, n, T, xi2_inv, sigma_inv, phi, tau2, prior_rhs_w, y
    );
  }
  return beta_update_global_conditional_cpp(
    X, xi2_inv, sigma_inv, y, w
  );
}

arma::vec w_update(
    const arma::vec& sigma_inv_site,
    const arma::vec& phi,
    const arma::vec& tau2,
    const arma::vec& y,
    const arma::vec& Zbeta,
    const arma::vec& prior_rhs_w
) {
  int n = sigma_inv_site.n_elem;
  int T = y.n_elem / n;
  arma::mat y_mat(const_cast<double*>(y.memptr()), n, T, false, false);
  arma::mat zbeta_mat(const_cast<double*>(Zbeta.memptr()), n, T, false, false);
  arma::mat prior_rhs_mat(const_cast<double*>(prior_rhs_w.memptr()), n, T, false, false);
  arma::vec w_new(n * T, arma::fill::zeros);
  arma::mat w_std_normals = arma::randn<arma::mat>(T, n);
  std::atomic<int> error_code(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int site = 0; site < n; ++site) {
    arma::vec rhs = sigma_inv_site(site) * (y_mat.row(site).t() - zbeta_mat.row(site).t()) +
      prior_rhs_mat.row(site).t();
    arma::vec diagonal(T, arma::fill::zeros);
    arma::vec cholesky_diagonal(T, arma::fill::zeros);
    arma::vec cholesky_subdiagonal(std::max(0, T - 1), arma::fill::zeros);
    arma::vec forward_mean(T, arma::fill::zeros);
    arma::vec w_mean(T, arma::fill::zeros);
    arma::vec w_noise(T, arma::fill::zeros);
    double inv_tau2 = 1.0 / tau2(site);
    double off_diagonal = -phi(site) * inv_tau2;

    if (T == 1) {
      diagonal(0) = (1.0 - phi(site) * phi(site)) * inv_tau2 +
        sigma_inv_site(site);
    } else {
      diagonal.fill((1.0 + phi(site) * phi(site)) * inv_tau2 +
        sigma_inv_site(site));
      diagonal(0) = inv_tau2 + sigma_inv_site(site);
      diagonal(T - 1) = inv_tau2 + sigma_inv_site(site);
    }

    bool factor_ok = std::isfinite(diagonal(0)) && diagonal(0) > 0.0;
    if (factor_ok) {
      cholesky_diagonal(0) = std::sqrt(diagonal(0));
      for (int t = 1; t < T; ++t) {
        cholesky_subdiagonal(t - 1) =
          off_diagonal / cholesky_diagonal(t - 1);
        double pivot = diagonal(t) -
          cholesky_subdiagonal(t - 1) * cholesky_subdiagonal(t - 1);
        if (!std::isfinite(pivot) || pivot <= 0.0) {
          factor_ok = false;
          break;
        }
        cholesky_diagonal(t) = std::sqrt(pivot);
      }
    }

    if (!factor_ok) {
      if (error_code.load(std::memory_order_relaxed) == 0) {
#ifdef _OPENMP
#pragma omp critical(w_update_error)
#endif
        {
          if (error_code.load(std::memory_order_relaxed) == 0) {
            error_code.store(1, std::memory_order_relaxed);
          }
        }
      }
      continue;
    }

    forward_mean(0) = rhs(0) / cholesky_diagonal(0);
    for (int t = 1; t < T; ++t) {
      forward_mean(t) =
        (rhs(t) - cholesky_subdiagonal(t - 1) * forward_mean(t - 1)) /
        cholesky_diagonal(t);
    }

    w_mean(T - 1) = forward_mean(T - 1) / cholesky_diagonal(T - 1);
    w_noise(T - 1) = w_std_normals(T - 1, site) /
      cholesky_diagonal(T - 1);
    for (int t = T - 2; t >= 0; --t) {
      w_mean(t) =
        (forward_mean(t) - cholesky_subdiagonal(t) * w_mean(t + 1)) /
        cholesky_diagonal(t);
      w_noise(t) =
        (w_std_normals(t, site) -
         cholesky_subdiagonal(t) * w_noise(t + 1)) /
        cholesky_diagonal(t);
    }

    arma::vec w_site = w_mean + w_noise;
    for (int t = 0; t < T; ++t) {
      w_new(t * n + site) = w_site(t);
    }
  }

  if (error_code.load(std::memory_order_relaxed) == 1) {
    Rcpp::stop("Sitewise w precision is not positive definite.");
  }

  return w_new;
}

////////////////////////////////////////////////////////////////////////////////

int sample_cpp(arma::uword n, int size, const arma::vec& probs, bool replace = true) {
  if (probs.n_elem != n) {
    Rcpp::stop("Length of probabilities vector must match n.");
  }
  // Convert arma::vec to Rcpp::NumericVector
  Rcpp::NumericVector rcpp_probs(probs.begin(), probs.end());
  
  // Perform sampling
  return Rcpp::sample(n, size, replace, rcpp_probs)[0];
}

////////////////////////////////////////////////////////////////////////////////
// Spatial cohesion and clustering helpers
////////////////////////////////////////////////////////////////////////////////

// Log-density for a diagonal-covariance Gaussian vector.
double ldmvnorm_diag(const arma::vec& x, const arma::vec& mean, const arma::vec& diag){
  int n = x.n_elem;  // Dimensionality of the multivariate normal
  double log_det = arma::sum(arma::trunc_log(diag));  // Log of the determinant of diagonal covariance
  double log_density = - 0.5 * n * std::log(2 * M_PI) - 0.5 * log_det;
  
  // Compute the exponent term (squared Mahalanobis distance)
  arma::vec diff = x - mean;
  double exponent = - 0.5 * arma::sum(arma::square(diff) / diag);
  
  // Combine the terms to get the log density
  log_density += exponent;
  
  // Return the density in the original scale by exponentiating
  return log_density;
}

////////////////////////////////////////////////////////////////////////////////

// Cohesion 1

double Cohesion1(const arma::vec& s1, const arma::vec& s2, double epsilon, int dim, int lg){
  
  double cent1 = 0.0; 
  double cent2 = 0.0;
  double out = 0.0;
  double sdist, maxdist, dist;
  
  for(int ii = 0; ii < dim; ii++){
    cent1 = cent1 + s1(ii)/(double) dim;
    cent2 = cent2 + s2(ii)/(double) dim;
  }
  
  maxdist = sqrt((s1(0) - cent1)*(s1(0) - cent1) +
    (s2(0) - cent2)*(s2(0) - cent2));
  
  sdist = 0;
  
  for(int ii = 0; ii < dim; ii++){
    dist = sqrt((s1(ii) - cent1)*(s1(ii) - cent1) +
      (s2(ii) - cent2)*(s2(ii) - cent2));
    sdist = sdist + dist;
    if(maxdist < dist) {
      maxdist = dist;
    } 
  }
  
  if (sdist >= 1) {
    out = -lgamma(sdist * epsilon);
  } else if (sdist != 0) {
    // out = -std::log(sdist);
    out = -arma::trunc_log(sdist);
  } else {
    out = log(1);
  }
  if (!lg) {
    // out = std::exp(out);
    out = arma::trunc_exp(out);
  } 
  return(out);
}

// Cohesion 2

double Cohesion2(const arma::vec& s1, const arma::vec& s2, double a, int dim, int lg){
  
  double dist = 0.0;
  double out = 0.0;
  
  for(int i = 0; i < dim; i++){
    for(int ii = 0; ii < dim; ++ii){
      dist = sqrt((s1(i) - s1(ii))*(s1(i) - s1(ii)) +
        (s2(i) - s2(ii))*(s2(i) - s2(ii)));
      if (dist < a) {
        out = 1.0;
      } else {
        out = 0.0;
        break;
      }
    }
    if (out < 1) {
      break;
    } 
  }
  
  if(lg) { 
    // out = std::log(out);
    out = arma::trunc_log(out);
  } 
  
  return(out);
  
}


// Cohesion 3 and 4

double G2a(double a, int lg){
  
  double out;
  
  out = log(M_PI) + lgamma(a) + lgamma(a - 0.5);
  
  if(!lg) {
    out = arma::trunc_exp(out);
  } 
  
  return(out);
}

double Cohesion3_4(const arma::vec& s1, const arma::vec& s2, const arma::vec& mu0, double k0, double v0, const arma::mat& L0,
                   int dim, int Cohesion, int lg){
  
  
  double kn, knn, vn, vnn, sbar1 = 0.0, sbar2 = 0.0, dL0, dLn, dLnn;
  double out = 0.0;
  double s_sbar1 = 0.0;
  double s_sbar2 = 0.0;
  double sbar_mu01, sbar_mu02, sbar_mun1, sbar_mun2, mun1, mun2;
  double Vs1 = 0.0;
  double Vs2 = 0.0;
  double Vs3 = 0.0;
  double Vs4 = 0.0;
  double Vsbarmu01, Vsbarmu02, Vsbarmu03, Vsbarmu04;
  double Vsbarmun1, Vsbarmun2, Vsbarmun3, Vsbarmun4;
  double Ln1, Ln2, Ln3, Ln4, Lnn1, Lnn2, Lnn3, Lnn4;
  
  for(int ii = 0; ii < dim; ii++){
    
    sbar1 = sbar1 + s1(ii)/(double) dim;
    sbar2 = sbar2 + s2(ii)/(double) dim;
    
  }
  
  for(int ii = 0; ii < dim; ii++){
    
    s_sbar1 = s1(ii) - sbar1;
    s_sbar2 = s2(ii) - sbar2;
    
    
    Vs1 = Vs1 + s_sbar1 * s_sbar1;
    Vs2 = Vs2 + s_sbar1 * s_sbar2;
    Vs3 = Vs3 + s_sbar2 * s_sbar1;
    Vs4 = Vs4 + s_sbar2 * s_sbar2;
    
  }
  
  kn = k0 + dim; vn = v0 + dim;
  knn = kn + dim; vnn = vn + dim;
  
  mun1 = k0 / (k0 + dim) * mu0(0) + dim / (k0 + dim) * sbar1;
  mun2 = k0 / (k0 + dim) * mu0(1) + dim / (k0 + dim) * sbar2;
  
  sbar_mu01 = sbar1 - mu0(0);
  sbar_mu02 = sbar2 - mu0(1);
  sbar_mun1 = sbar1 - mun1;
  sbar_mun2 = sbar2 - mun2;
  
  Vsbarmu01 = sbar_mu01 * sbar_mu01, Vsbarmun1 = sbar_mun1 * sbar_mun1;
  Vsbarmu02 = sbar_mu01 * sbar_mu02, Vsbarmun2 = sbar_mun1 * sbar_mun2;
  Vsbarmu03 = sbar_mu02 * sbar_mu01, Vsbarmun3 = sbar_mun2 * sbar_mun1;
  Vsbarmu04 = sbar_mu02 * sbar_mu02, Vsbarmun4 = sbar_mun2 * sbar_mun2;
  
  Ln1 = L0(0, 0) + Vs1 + k0 * dim / (k0 + dim) * Vsbarmu01;
  Ln2 = L0(0, 1) + Vs2 + k0 * dim / (k0 + dim) * Vsbarmu02;
  Ln3 = L0(1, 0) + Vs3 + k0 * dim / (k0 + dim) * Vsbarmu03;
  Ln4 = L0(1, 1) + Vs4 + k0 * dim / (k0 + dim) * Vsbarmu04;
  
  Lnn1 = Ln1 + Vs1 + kn * dim / (kn + dim) * Vsbarmun1;
  Lnn2 = Ln2 + Vs2 + kn * dim / (kn + dim) * Vsbarmun2;
  Lnn3 = Ln3 + Vs3 + kn * dim / (kn + dim) * Vsbarmun3;
  Lnn4 = Ln4 + Vs4 + kn * dim / (kn + dim) * Vsbarmun4;
  
  dL0 = L0(0, 0) * L0(1, 1) - L0(0, 1) * L0(1, 0);
  dLn = Ln1 * Ln4 - Ln2 * Ln3;
  dLnn = Lnn1 * Lnn4 - Lnn2 * Lnn3;
  
  if(Cohesion == 3){
    out = - dim * log(M_PI) +
      (G2a(0.5 * vn, 1) - G2a(0.5 * v0, 1)) +
      (0.5 * v0 * arma::trunc_log(dL0) - 0.5 * vn * arma::trunc_log(dLn)) +
      (arma::trunc_log(k0) - arma::trunc_log(kn));
  }
  if(Cohesion == 4){
    out = - dim * log(M_PI) +
      (G2a(0.5 * vnn, 1) - G2a(0.5 * vn, 1)) +
      (0.5 * vn * arma::trunc_log(dLn) - 0.5 * vnn * arma::trunc_log(dLnn)) +
      (arma::trunc_log(kn) - arma::trunc_log(knn));
  }
  if(!lg) {
    out = arma::trunc_exp(out);  
  }
  
  return(out);
  
}

double compute_cohesion_log_value(
    int cohesion,
    const arma::vec& s1,
    const arma::vec& s2,
    const arma::vec& mu0,
    double k0,
    double v0,
    const arma::mat& L0,
    double eC1,
    double aC2,
    int dim
) {
  if (cohesion == 1) return Cohesion1(s1, s2, eC1, dim, 1);
  if (cohesion == 2) return Cohesion2(s1, s2, aC2, dim, 1);
  if (cohesion == 3) return Cohesion3_4(s1, s2, mu0, k0, v0, L0, dim, 3, 1);
  if (cohesion == 4) return Cohesion3_4(s1, s2, mu0, k0, v0, L0, dim, 4, 1);
  return 0.0;
}

struct SpatialClusterStatistics {
  int count;
  double sum1;
  double sum2;
  double sum11;
  double sum12;
  double sum22;
};

SpatialClusterStatistics empty_spatial_cluster_statistics_cpp() {
  return SpatialClusterStatistics{0, 0.0, 0.0, 0.0, 0.0, 0.0};
}

void update_spatial_cluster_statistics_cpp(
    SpatialClusterStatistics& statistics,
    double coordinate1,
    double coordinate2,
    int direction
) {
  statistics.count += direction;
  statistics.sum1 += direction * coordinate1;
  statistics.sum2 += direction * coordinate2;
  statistics.sum11 += direction * coordinate1 * coordinate1;
  statistics.sum12 += direction * coordinate1 * coordinate2;
  statistics.sum22 += direction * coordinate2 * coordinate2;
  if (statistics.count < 0) {
    Rcpp::stop("Invalid cached spatial-cluster count.");
  }
}

double cohesion3_4_from_statistics_cpp(
    const SpatialClusterStatistics& statistics,
    const arma::vec& mu0,
    double k0,
    double v0,
    const arma::mat& L0,
    int cohesion
) {
  const int dim = statistics.count;
  if (dim <= 0) {
    return 0.0;
  }

  double sbar1 = statistics.sum1 / static_cast<double>(dim);
  double sbar2 = statistics.sum2 / static_cast<double>(dim);
  double Vs1 = std::max(
    0.0,
    statistics.sum11 - statistics.sum1 * statistics.sum1 /
      static_cast<double>(dim)
  );
  double Vs2 = statistics.sum12 - statistics.sum1 * statistics.sum2 /
    static_cast<double>(dim);
  double Vs4 = std::max(
    0.0,
    statistics.sum22 - statistics.sum2 * statistics.sum2 /
      static_cast<double>(dim)
  );

  double kn = k0 + dim;
  double vn = v0 + dim;
  double knn = kn + dim;
  double vnn = vn + dim;
  double mun1 = k0 / (k0 + dim) * mu0(0) +
    dim / (k0 + dim) * sbar1;
  double mun2 = k0 / (k0 + dim) * mu0(1) +
    dim / (k0 + dim) * sbar2;
  double sbar_mu01 = sbar1 - mu0(0);
  double sbar_mu02 = sbar2 - mu0(1);
  double sbar_mun1 = sbar1 - mun1;
  double sbar_mun2 = sbar2 - mun2;

  double Ln1 = L0(0, 0) + Vs1 +
    k0 * dim / (k0 + dim) * sbar_mu01 * sbar_mu01;
  double Ln2 = L0(0, 1) + Vs2 +
    k0 * dim / (k0 + dim) * sbar_mu01 * sbar_mu02;
  double Ln3 = L0(1, 0) + Vs2 +
    k0 * dim / (k0 + dim) * sbar_mu02 * sbar_mu01;
  double Ln4 = L0(1, 1) + Vs4 +
    k0 * dim / (k0 + dim) * sbar_mu02 * sbar_mu02;

  double Lnn1 = Ln1 + Vs1 +
    kn * dim / (kn + dim) * sbar_mun1 * sbar_mun1;
  double Lnn2 = Ln2 + Vs2 +
    kn * dim / (kn + dim) * sbar_mun1 * sbar_mun2;
  double Lnn3 = Ln3 + Vs2 +
    kn * dim / (kn + dim) * sbar_mun2 * sbar_mun1;
  double Lnn4 = Ln4 + Vs4 +
    kn * dim / (kn + dim) * sbar_mun2 * sbar_mun2;

  double dL0 = L0(0, 0) * L0(1, 1) - L0(0, 1) * L0(1, 0);
  double dLn = Ln1 * Ln4 - Ln2 * Ln3;
  double dLnn = Lnn1 * Lnn4 - Lnn2 * Lnn3;
  if (dL0 <= 0.0 || dLn <= 0.0 || dLnn <= 0.0) {
    return -arma::datum::inf;
  }

  if (cohesion == 3) {
    return -dim * std::log(M_PI) +
      (G2a(0.5 * vn, 1) - G2a(0.5 * v0, 1)) +
      (0.5 * v0 * std::log(dL0) - 0.5 * vn * std::log(dLn)) +
      (std::log(k0) - std::log(kn));
  }
  return -dim * std::log(M_PI) +
    (G2a(0.5 * vnn, 1) - G2a(0.5 * vn, 1)) +
    (0.5 * vn * std::log(dLn) - 0.5 * vnn * std::log(dLnn)) +
    (std::log(kn) - std::log(knn));
}

std::vector<SpatialClusterStatistics> build_spatial_cluster_statistics_cpp(
    const arma::ivec& cluster,
    int n_clusters,
    const arma::vec& coordinate1,
    const arma::vec& coordinate2
) {
  std::vector<SpatialClusterStatistics> statistics(
    n_clusters,
    empty_spatial_cluster_statistics_cpp()
  );
  for (arma::uword site = 0; site < cluster.n_elem; ++site) {
    int cluster_index = cluster(site) - 1;
    if (cluster_index < 0 || cluster_index >= n_clusters) {
      Rcpp::stop("Invalid cluster label while caching spatial statistics.");
    }
    update_spatial_cluster_statistics_cpp(
      statistics[cluster_index],
      coordinate1(site),
      coordinate2(site),
      1
    );
  }
  return statistics;
}

double cached_cohesion_log_value_cpp(
    int cohesion,
    const SpatialClusterStatistics& statistics,
    const arma::vec& mu0,
    double k0,
    double v0,
    const arma::mat& L0
) {
  if (cohesion == 3 || cohesion == 4) {
    return cohesion3_4_from_statistics_cpp(
      statistics,
      mu0,
      k0,
      v0,
      L0,
      cohesion
    );
  }
  Rcpp::stop("Cached spatial statistics are available only for cohesion 3 or 4.");
}

struct PredictionNeighborGeometry {
  arma::uvec neighbors;
  arma::vec cross_distances;
  arma::mat neighbor_distances;
};

std::vector<PredictionNeighborGeometry> build_prediction_geometry_cpp(
    const arma::mat& locations,
    const arma::mat& prediction_locations,
    const arma::umat& prediction_neighbors
) {
  std::vector<PredictionNeighborGeometry> geometry(prediction_locations.n_rows);
  for (arma::uword prediction = 0;
       prediction < prediction_locations.n_rows;
       ++prediction) {
    arma::uvec neighbors_1based = prediction_neighbors.row(prediction).t();
    neighbors_1based = neighbors_1based.elem(arma::find(neighbors_1based > 0));
    if (neighbors_1based.n_elem == 0) {
      Rcpp::stop("Every prediction location requires at least one observed neighbor.");
    }
    if (arma::any(neighbors_1based > locations.n_rows)) {
      Rcpp::stop("Prediction-neighbor indices must refer to observed locations.");
    }
    arma::uvec neighbors = neighbors_1based - 1;
    PredictionNeighborGeometry item;
    item.neighbors = neighbors;
    item.cross_distances.set_size(neighbors.n_elem);
    item.neighbor_distances.zeros(neighbors.n_elem, neighbors.n_elem);
    for (arma::uword first = 0; first < neighbors.n_elem; ++first) {
      item.cross_distances(first) = arma::norm(
        prediction_locations.row(prediction) - locations.row(neighbors(first)),
        2
      );
      for (arma::uword second = first + 1;
           second < neighbors.n_elem;
           ++second) {
        double distance = arma::norm(
          locations.row(neighbors(first)) - locations.row(neighbors(second)),
          2
        );
        item.neighbor_distances(first, second) = distance;
        item.neighbor_distances(second, first) = distance;
      }
    }
    geometry[prediction] = item;
  }
  return geometry;
}

void prediction_nngp_factor_cpp(
    const PredictionNeighborGeometry& geometry,
    double sigma2_u,
    double rho,
    double nu,
    arma::vec& coefficients,
    double& conditional_variance
) {
  arma::uword neighbors = geometry.neighbors.n_elem;
  arma::mat covariance(neighbors, neighbors, arma::fill::zeros);
  arma::vec cross_covariance(neighbors, arma::fill::zeros);
  double marginal_sd = std::sqrt(sigma2_u);
  for (arma::uword first = 0; first < neighbors; ++first) {
    covariance(first, first) = sigma2_u;
    cross_covariance(first) = matern_cov_scalar(
      geometry.cross_distances(first), nu, rho, marginal_sd, marginal_sd
    );
    for (arma::uword second = first + 1; second < neighbors; ++second) {
      double value = matern_cov_scalar(
        geometry.neighbor_distances(first, second),
        nu,
        rho,
        marginal_sd,
        marginal_sd
      );
      covariance(first, second) = value;
      covariance(second, first) = value;
    }
  }

  double base_jitter = std::max(1e-8 * sigma2_u, 1e-12);
  double minimum_conditional_variance = 0.5 * base_jitter;
  for (int attempt = 0; attempt < 7; ++attempt) {
    arma::mat candidate = covariance;
    if (attempt > 0) {
      candidate.diag() += base_jitter * std::pow(10.0, attempt - 1);
    }
    bool success = arma::solve(
      coefficients,
      candidate,
      cross_covariance,
      arma::solve_opts::fast + arma::solve_opts::no_approx
    );
    if (!success || !coefficients.is_finite()) {
      continue;
    }
    conditional_variance = sigma2_u - arma::dot(
      cross_covariance,
      coefficients
    );
    if (std::isfinite(conditional_variance) &&
        conditional_variance >= minimum_conditional_variance) {
      return;
    }
  }
  Rcpp::stop("Unable to construct a positive prediction conditional variance.");
}

double empirical_quantile_sorted_cpp(
    const std::vector<float>& sorted,
    double probability
) {
  if (sorted.empty()) return NA_REAL;
  if (sorted.size() == 1) return sorted[0];
  double index = probability * static_cast<double>(sorted.size() - 1);
  std::size_t lower = static_cast<std::size_t>(std::floor(index));
  std::size_t upper = static_cast<std::size_t>(std::ceil(index));
  double fraction = index - static_cast<double>(lower);
  return (1.0 - fraction) * sorted[lower] + fraction * sorted[upper];
}

List summarize_prediction_draws_cpp(
    const arma::fmat& draws,
    int prediction_times,
    int prediction_locations,
    SEXP truth_sexp
) {
  int retained = static_cast<int>(draws.n_rows);
  int targets = prediction_times * prediction_locations;
  arma::mat mean(prediction_times, prediction_locations, arma::fill::zeros);
  arma::mat sd(prediction_times, prediction_locations, arma::fill::zeros);
  arma::mat q025(prediction_times, prediction_locations, arma::fill::zeros);
  arma::mat q975(prediction_times, prediction_locations, arma::fill::zeros);
  bool has_truth = !Rf_isNull(truth_sexp);
  arma::mat truth;
  if (has_truth) {
    truth = Rcpp::as<arma::mat>(truth_sexp);
    if (truth.n_rows != static_cast<arma::uword>(prediction_times) ||
        truth.n_cols != static_cast<arma::uword>(prediction_locations) ||
        !truth.is_finite()) {
      Rcpp::stop("Prediction truth matrices must be finite T x n_prediction matrices.");
    }
  }
  arma::mat squared_error;
  arma::mat absolute_error;
  arma::mat coverage;
  arma::mat interval_width;
  arma::mat crps;
  if (has_truth) {
    squared_error.zeros(prediction_times, prediction_locations);
    absolute_error.zeros(prediction_times, prediction_locations);
    coverage.zeros(prediction_times, prediction_locations);
    interval_width.zeros(prediction_times, prediction_locations);
    crps.zeros(prediction_times, prediction_locations);
  }

  std::vector<float> values(retained);
  for (int target = 0; target < targets; ++target) {
    double sum = 0.0;
    double sum_squares = 0.0;
    for (int draw = 0; draw < retained; ++draw) {
      float value = draws(draw, target);
      values[draw] = value;
      sum += value;
      sum_squares += static_cast<double>(value) * value;
    }
    std::sort(values.begin(), values.end());
    int time = target / prediction_locations;
    int location = target % prediction_locations;
    double posterior_mean = sum / static_cast<double>(retained);
    double variance = retained > 1 ?
      (sum_squares - sum * sum / retained) / (retained - 1.0) : 0.0;
    mean(time, location) = posterior_mean;
    sd(time, location) = std::sqrt(std::max(0.0, variance));
    q025(time, location) = empirical_quantile_sorted_cpp(values, 0.025);
    q975(time, location) = empirical_quantile_sorted_cpp(values, 0.975);

    if (has_truth) {
      double observed = truth(time, location);
      double error = posterior_mean - observed;
      squared_error(time, location) = error * error;
      absolute_error(time, location) = std::abs(error);
      coverage(time, location) = observed >= q025(time, location) &&
        observed <= q975(time, location);
      interval_width(time, location) = q975(time, location) - q025(time, location);
      double first_term = 0.0;
      double second_term = 0.0;
      for (int draw = 0; draw < retained; ++draw) {
        first_term += std::abs(static_cast<double>(values[draw]) - observed);
        second_term += (2.0 * (draw + 1) - retained - 1.0) * values[draw];
      }
      crps(time, location) =
        first_term / retained - second_term / (retained * retained);
    }
  }

  List result = List::create(
    Named("mean") = mean,
    Named("sd") = sd,
    Named("q025") = q025,
    Named("q975") = q975,
    Named("retained_draws") = retained
  );
  if (has_truth) {
    result["squared_error"] = squared_error;
    result["absolute_error"] = absolute_error;
    result["coverage"] = coverage;
    result["interval_width"] = interval_width;
    result["crps"] = crps;
  }
  return result;
}

double evaluate_candidate_loglik(
    const List& precomp,
    const arma::vec& w,
    const arma::vec& phi,
    const arma::vec& tau2,
    const arma::vec& current_F,
    const arma::vec& current_values,
    const arma::mat& current_w_tilde,
    const arma::mat& current_innovations,
    double current_quad_form,
    double current_log_phi_term,
    arma::uword update_idx,
    double phi_new,
    double tau2_new,
    double rho,
    double nu,
    int n,
    int T
) {
  arma::vec tau2_cand = tau2;
  tau2_cand(update_idx) = tau2_new;
  arma::vec f_cand = current_F;
  arma::vec values_cand = current_values;
  update_Q_fast_single_row_cpp(precomp, tau2_cand, update_idx, rho, nu, f_cand, values_cand);

  List dependents = as<List>(precomp["dependents"]);
  arma::uvec affected_rows = as<arma::uvec>(dependents[update_idx]);
  arma::rowvec updated_row_tilde = compute_w_tilde_row_cpp(w, phi, update_idx, n, T, phi_new);

  double quad_form = current_quad_form;
  for (arma::uword idx = 0; idx < affected_rows.n_elem; ++idx) {
    arma::uword row = affected_rows(idx);
    quad_form -= arma::dot(current_innovations.row(row), current_innovations.row(row)) / current_F(row);
    arma::rowvec innovation_row = compute_innovation_row_cpp(
      row, current_w_tilde, precomp, values_cand, &updated_row_tilde, update_idx
    );
    quad_form += arma::dot(innovation_row, innovation_row) / f_cand(row);
  }

  double log_phi_term = current_log_phi_term -
    0.5 * std::log(1.0 - phi(update_idx) * phi(update_idx)) +
    0.5 * std::log(1.0 - phi_new * phi_new);

  double log_det_Q = -arma::sum(arma::log(f_cand));
  return log_phi_term + (T / 2.0) * log_det_Q - 0.5 * quad_form;
}

double compute_exact_loglik_w_cpp(
    const arma::vec& w,
    const List& precomp,
    const arma::vec& phi,
    const arma::vec& tau2,
    double rho,
    double nu,
    int n,
    int T
) {
  arma::vec f_prop;
  arma::vec values_prop;
  compute_Q_fast_cpp(precomp, tau2, rho, nu, f_prop, values_prop);

  arma::mat w_tilde_prop;
  arma::mat innovation_prop;
  build_w_tilde_cache_cpp(w, phi, n, T, w_tilde_prop);
  build_innovation_cache_cpp(w_tilde_prop, precomp, values_prop, innovation_prop);

  double log_phi_term = compute_log_phi_term_cpp(phi);
  double log_det_term = -static_cast<double>(T) * 0.5 * arma::sum(arma::log(f_prop));
  double quad_form = compute_quad_from_innovations_cpp(innovation_prop, f_prop);

  return log_phi_term + log_det_term - 0.5 * quad_form;
}

List unique_columns_info(const arma::mat& gamma) {
  int n = gamma.n_cols; // Number of columns in gamma
  
  // Containers for unique column indices, cardinalities, and labels
  arma::uvec unique_indices;
  arma::uvec cardinalities; // Cardinality of each unique column
  arma::uvec labels(n, arma::fill::zeros); // Group labels for each column
  
  // Iterate over all columns to find unique ones
  for (int j = 0; j < n; j++) {
    // Extract the current column
    arma::vec current_col = gamma.col(j);
    bool is_unique = true;
    
    // Check if this column matches any of the already identified unique columns
    for (size_t k = 0; k < unique_indices.n_elem; k++) {
      if (arma::all(gamma.col(unique_indices(k)) == current_col)) {
        is_unique = false;
        labels(j) = k + 1; // Assign the group label (1-based index)
        ++cardinalities(k); // Increment the cardinality of the matched group
        break;
      }
    }
    
    // If the column is unique, add it to the list of unique columns
    if (is_unique) {
      unique_indices.insert_rows(unique_indices.n_elem, arma::uvec({(unsigned int)j}));
      cardinalities.insert_rows(cardinalities.n_elem, arma::uvec({1})); // Initialize cardinality to 1
      labels(j) = unique_indices.n_elem; // Assign a new group label
    }
  }
  
  // Extract the unique columns using the indices
  arma::mat unique_columns = gamma.cols(unique_indices);
  
  // Return the results as a list
  return List::create(
    Named("unique_columns") = unique_columns,
    Named("cardinalities") = cardinalities,
    Named("labels") = labels
  );
}


arma::vec solve_projection_coefficients_cpp(
    const arma::mat& crossproduct,
    const arma::vec& rhs,
    const char* label
) {
  arma::vec coefficients;
  bool solved = arma::solve(
    coefficients,
    crossproduct,
    rhs,
    arma::solve_opts::likely_sympd
  );
  if (!solved || !coefficients.is_finite()) {
    coefficients = arma::pinv(crossproduct) * rhs;
  }
  if (!coefficients.is_finite()) {
    Rcpp::stop(std::string(label) + " projection failed.");
  }
  return coefficients;
}

void compute_latent_design_projections_cpp(
    const arma::mat& X,
    const arma::mat& design_crossproduct,
    const arma::vec& w,
    const arma::vec& u,
    const arma::vec& phi,
    int n,
    int T,
    arma::vec& w_projection,
    arma::vec& u_transition_projection
) {
  int p = X.n_cols;
  w_projection = solve_projection_coefficients_cpp(
    design_crossproduct,
    X.t() * w,
    "Latent-state"
  );
  arma::mat transition_crossproduct(p, p, arma::fill::zeros);
  arma::vec transition_rhs(p, arma::fill::zeros);
  for (int time = 1; time < T; ++time) {
    arma::mat lagged_design = X.rows((time - 1) * n, time * n - 1);
    lagged_design.each_col() %= phi;
    arma::mat transition_design =
      X.rows(time * n, (time + 1) * n - 1) - lagged_design;
    transition_crossproduct += transition_design.t() * transition_design;
    transition_rhs += transition_design.t() *
      u.subvec(time * n, (time + 1) * n - 1);
  }
  u_transition_projection = solve_projection_coefficients_cpp(
    arma::symmatu(transition_crossproduct),
    transition_rhs,
    "Common-shock transition"
  );
}

////////////////////////////////////////////////////////////////////////////////
// Main sampler entry point
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List dcs_sppm_global_beta_cpp(arma::mat Y, arma::mat locations, arma::mat X,
                                         int runs, int burn, int thin,
                                         int cohesion,
                                         arma::vec init_beta,
                                         arma::vec init_xi2, arma::vec init_sigma2,
                                         arma::vec init_phi, arma::vec init_tau2, double init_M,
                                         double init_kappa, double nu,
                                         arma::umat NN_ind,
                                         double k0 = 1.0, double v0 = 2.0,
                                         double eC1 = 1.0, double aC2 = 1.5,
                                         double a_sigma = 1.0, double b_sigma = 1.0, double a_xi = 1.0, double b_xi = 1.0,
                                         double a_tau = 1.0, double b_tau = 1.0, double a_phi = 1.0, double b_phi = 1.0,
                                         double a_kappa = 1.0, double b_kappa = 1.0,
                                         double eta_phi = 1.0, double eta_kappa = 1.0,
                                         int CC = 3,
                                         std::string beta_update_mode = "collapsed",
                                         double init_sigma2_u = 1.0,
                                         double a_sigma2_u = 1.0, double b_sigma2_u = 1.0,
                                         bool debug_checks = false,
                                         SEXP locations_similarity = R_NilValue,
                                         bool save_latent = true,
                                         bool save_w = true,
                                         SEXP fixed_u = R_NilValue,
                                         std::string sigma2_update_mode = "conditional",
                                         double eta_sigma2 = 0.15,
                                         std::string tau2_update_mode = "conditional",
                                         double eta_tau2 = 0.15,
                                         double joint_variance_initial_scale = 0.15,
                                         std::string spatial_update_mode = "conditional",
                                         double spatial_initial_scale = 0.08,
                                         std::string spatial_process_model = "nngp",
                                         std::string progress_file = "",
                                         bool show_progress = true,
                                         SEXP prediction_locations = R_NilValue,
                                         SEXP prediction_locations_similarity = R_NilValue,
                                         SEXP prediction_neighbors = R_NilValue,
                                         SEXP prediction_X = R_NilValue,
                                         int prediction_thin = 10,
                                         int prediction_seed = 20400001,
                                         SEXP prediction_truth_u = R_NilValue,
                                         SEXP prediction_truth_w = R_NilValue,
                                         SEXP prediction_truth_y = R_NilValue,
                                         SEXP prediction_true_cluster = R_NilValue,
                                         SEXP prediction_cluster_phi = R_NilValue,
                                         SEXP prediction_cluster_tau2 = R_NilValue) {

  int T = Y.n_rows;
  int n = Y.n_cols;
  int p = X.n_cols;

  if (T < 2 || n < 2 || p < 1) {
    Rcpp::stop("Require T >= 2, n >= 2, and p >= 1.");
  }
  if (locations.n_rows != static_cast<arma::uword>(n) ||
      locations.n_cols != 2) {
    Rcpp::stop("`locations` must be an n x 2 matrix.");
  }
  if (X.n_rows != static_cast<arma::uword>(n * T)) {
    Rcpp::stop("`X` must have n * T rows ordered by time blocks.");
  }
  if (!X.is_finite()) {
    Rcpp::stop("`X` must contain only finite values.");
  }
  bool condition_on_fixed_u = !Rf_isNull(fixed_u);
  arma::vec fixed_u_values;
  if (condition_on_fixed_u) {
    if (!Rf_isMatrix(fixed_u)) {
      Rcpp::stop("`fixed_u` must be a numeric T x n matrix when supplied.");
    }
    arma::mat fixed_u_matrix = Rcpp::as<arma::mat>(fixed_u);
    if (fixed_u_matrix.n_rows != static_cast<arma::uword>(T) ||
        fixed_u_matrix.n_cols != static_cast<arma::uword>(n) ||
        !fixed_u_matrix.is_finite()) {
      Rcpp::stop("`fixed_u` must be a finite T x n matrix when supplied.");
    }
    fixed_u_values = arma::vectorise(fixed_u_matrix.t());
  }
  if (init_beta.n_elem != static_cast<arma::uword>(p) ||
      init_xi2.n_elem != static_cast<arma::uword>(p)) {
    Rcpp::stop("`init_beta` and `init_xi2` must both have length p.");
  }
  if (init_phi.n_elem != static_cast<arma::uword>(n) ||
      init_tau2.n_elem != static_cast<arma::uword>(n)) {
    Rcpp::stop("`init_phi` and `init_tau2` must both have length n.");
  }
  if (runs <= burn || burn < 0 || thin <= 0) {
    Rcpp::stop("Require runs > burn >= 0 and thin > 0.");
  }
  if (prediction_thin <= 0) {
    Rcpp::stop("`prediction_thin` must be positive.");
  }
  if (!std::isfinite(eta_sigma2) || eta_sigma2 <= 0.0 ||
      !std::isfinite(eta_tau2) || eta_tau2 <= 0.0 ||
      !std::isfinite(joint_variance_initial_scale) ||
      joint_variance_initial_scale <= 0.0 ||
      !std::isfinite(spatial_initial_scale) || spatial_initial_scale <= 0.0) {
    Rcpp::stop("All variance proposal scales must be positive and finite.");
  }
  if (init_M <= 0.0 || init_kappa <= 0.0 || nu <= 0.0 ||
      init_sigma2_u <= 0.0) {
    Rcpp::stop("Initial mass, kappa, nu, and sigma2_u must be positive.");
  }

  arma::mat similarity_locations = locations;
  if (!Rf_isNull(locations_similarity)) {
    if (!Rf_isMatrix(locations_similarity)) {
      Rcpp::stop("`locations_similarity` must be a numeric matrix.");
    }
    similarity_locations = Rcpp::as<arma::mat>(locations_similarity);
    if (similarity_locations.n_rows != locations.n_rows ||
        similarity_locations.n_cols != locations.n_cols) {
      Rcpp::stop("`locations_similarity` must have the same dimensions as `locations`.");
    }
  }

  arma::vec mu0 = trans(arma::mean(similarity_locations, 0));
  arma::mat L0 = arma::eye(similarity_locations.n_cols, similarity_locations.n_cols);
  double Mdp = init_M;

  bool prediction_enabled = !Rf_isNull(prediction_locations);
  arma::mat prediction_location_matrix;
  arma::mat prediction_similarity_matrix;
  arma::umat prediction_neighbor_matrix;
  arma::mat prediction_design;
  int n_prediction = 0;
  std::vector<PredictionNeighborGeometry> prediction_geometry;
  if (prediction_enabled) {
    if (Rf_isNull(prediction_locations_similarity) ||
        Rf_isNull(prediction_neighbors) || Rf_isNull(prediction_X)) {
      Rcpp::stop(
        "Prediction requires locations, similarity locations, neighbors, and X."
      );
    }
    prediction_location_matrix = Rcpp::as<arma::mat>(prediction_locations);
    prediction_similarity_matrix = Rcpp::as<arma::mat>(
      prediction_locations_similarity
    );
    prediction_neighbor_matrix = Rcpp::as<arma::umat>(prediction_neighbors);
    prediction_design = Rcpp::as<arma::mat>(prediction_X);
    n_prediction = static_cast<int>(prediction_location_matrix.n_rows);
    if (n_prediction < 1 || prediction_location_matrix.n_cols != 2 ||
        prediction_similarity_matrix.n_rows !=
          prediction_location_matrix.n_rows ||
        prediction_similarity_matrix.n_cols != 2 ||
        prediction_neighbor_matrix.n_rows != prediction_location_matrix.n_rows ||
        prediction_design.n_rows !=
          static_cast<arma::uword>(T * n_prediction) ||
        prediction_design.n_cols != static_cast<arma::uword>(p) ||
        !prediction_location_matrix.is_finite() ||
        !prediction_similarity_matrix.is_finite() ||
        !prediction_design.is_finite()) {
      Rcpp::stop("Invalid dimensions or values in prediction inputs.");
    }
    prediction_geometry = build_prediction_geometry_cpp(
      locations,
      prediction_location_matrix,
      prediction_neighbor_matrix
    );
    if (cohesion != 3 && cohesion != 4) {
      Rcpp::stop("The prediction implementation currently requires cohesion 3 or 4.");
    }
  } else if (!Rf_isNull(prediction_locations_similarity) ||
             !Rf_isNull(prediction_neighbors) || !Rf_isNull(prediction_X)) {
    Rcpp::stop("Supply all prediction inputs or none of them.");
  }

  arma::vec y = arma::vectorise(Y.t());
  arma::vec y_observed = y;
  arma::uvec missing_indices(y.n_elem);
  arma::uvec observed_indices(y.n_elem);
  arma::uword n_missing = 0;
  arma::uword n_observed = 0;

  for (arma::uword idx = 0; idx < y.n_elem; ++idx) {
    if (std::isnan(y(idx))) {
      missing_indices(n_missing++) = idx;
      y(idx) = 0.0;
    } else if (!std::isfinite(y(idx))) {
      Rcpp::stop("`Y` contains an infinite value.");
    } else {
      observed_indices(n_observed++) = idx;
    }
  }
  missing_indices.resize(n_missing);
  observed_indices.resize(n_observed);
  if (n_observed == 0) {
    Rcpp::stop("`Y` must contain at least one observed value.");
  }

  arma::vec s1 = similarity_locations.col(1);
  arma::vec s2 = similarity_locations.col(0);
  arma::vec s1o(n, arma::fill::zeros);
  arma::vec s2o(n, arma::fill::zeros);
  arma::vec s1n(n, arma::fill::zeros);
  arma::vec s2n(n, arma::fill::zeros);

  if (init_sigma2.n_elem < 1) {
    Rcpp::stop("`init_sigma2` must contain at least one initial value.");
  }
  double sigma2 = init_sigma2(0);
  if (!std::isfinite(sigma2) || sigma2 <= 0.0 ||
      arma::any(init_xi2 <= 0.0) || arma::any(init_tau2 <= 0.0) ||
      arma::any(init_phi <= 0.0) || arma::any(init_phi >= 1.0)) {
    Rcpp::stop("Initial variances must be positive and initial phi values must lie in (0, 1).");
  }
  arma::vec sigma_inv_diag(n, arma::fill::value(1.0 / sigma2));
  arma::vec xi2 = init_xi2;
  arma::vec xi2_inv = 1.0 / xi2;
  arma::vec beta = init_beta;

  std::transform(
    spatial_process_model.begin(),
    spatial_process_model.end(),
    spatial_process_model.begin(),
    [](unsigned char character) { return std::tolower(character); }
  );
  bool use_exact_spatial_process = spatial_process_model == "exact_gp";
  if (!use_exact_spatial_process && spatial_process_model != "nngp") {
    Rcpp::stop("`spatial_process_model` must be 'nngp' or 'exact_gp'.");
  }
  List precomp = use_exact_spatial_process ?
    List::create() : precompute_NNGP_cpp(locations, NN_ind);

  arma::vec phi = init_phi;
  arma::vec tau2 = init_tau2;
  double rho = init_kappa;
  double a_rho = a_kappa;
  double b_rho = b_kappa;
  double eta_rho = eta_kappa;
  double sigma2_u = init_sigma2_u;
  arma::vec F_u;
  arma::vec values_u;
  auto build_spatial_precision = [&](double variance, double range,
                                     arma::vec& conditional_variances,
                                     arma::vec& conditional_coefficients) {
    if (use_exact_spatial_process) {
      return build_exact_gp_state_cpp(
        locations, variance, range, nu,
        conditional_variances, conditional_coefficients
      );
    }
    return build_Tauinv_u_cpp(
      precomp, n, variance, range, nu,
      conditional_variances, conditional_coefficients
    );
  };
  arma::sp_mat Tauinv_u = build_spatial_precision(
    sigma2_u, rho, F_u, values_u
  );

  arma::vec phi_star;
  arma::vec tau2_star;
  int k_current = 0;
  arma::ivec nj_current;
  arma::ivec clu_current;

  arma::mat gamma(2, n);
  gamma.row(0) = phi.t();
  gamma.row(1) = tau2.t();
  List summary_gamma = unique_columns_info(gamma);
  arma::mat gamma_star = summary_gamma["unique_columns"];
  phi_star = gamma_star.row(0).t();
  tau2_star = gamma_star.row(1).t();
  k_current = gamma_star.n_cols;
  arma::uvec nj_current_u = Rcpp::as<arma::uvec>(summary_gamma["cardinalities"]);
  arma::uvec clu_current_u = Rcpp::as<arma::uvec>(summary_gamma["labels"]);
  nj_current = arma::conv_to<arma::ivec>::from(nj_current_u);
  clu_current = arma::conv_to<arma::ivec>::from(clu_current_u);

  CholmodSampler u_sampler;
  u_sampler.analyze_pattern(build_u_precision_cpp(Tauinv_u, tau2));
  CholmodSampler collapsed_spatial_sampler;
  collapsed_spatial_sampler.analyze_pattern(build_u_precision_cpp(Tauinv_u, tau2));

  bool use_collapsed_beta = use_collapsed_beta_update(beta_update_mode);
  bool use_collapsed_sigma2 = use_collapsed_sigma2_update(sigma2_update_mode);
  bool use_collapsed_spatial = use_collapsed_spatial_update(spatial_update_mode);
  int tau2_update = tau2_update_code(tau2_update_mode);
  bool use_collapsed_tau2 = tau2_update == 1;
  bool use_joint_collapsed_variance = tau2_update == 2;
  bool use_pairwise_interweaving = tau2_update == 4;
  bool use_pairwise_collapsed_variance =
    tau2_update == 3 || use_pairwise_interweaving;
  if ((use_joint_collapsed_variance || use_pairwise_collapsed_variance) &&
      !use_collapsed_sigma2) {
    Rcpp::stop("Blocked collapsing requires `sigma2_update_mode = 'collapsed'`.");
  }
  if ((use_collapsed_sigma2 || use_collapsed_tau2 || use_joint_collapsed_variance ||
       use_pairwise_collapsed_variance) &&
      n_missing > 0) {
    Rcpp::stop("The partially collapsed validation updates currently require complete Y.");
  }
  if (use_collapsed_spatial && condition_on_fixed_u) {
    Rcpp::stop("The collapsed spatial update cannot be used with `fixed_u`.");
  }
  if (use_exact_spatial_process && !use_collapsed_spatial) {
    Rcpp::stop(
      "The exact GP backend currently requires `spatial_update_mode = 'collapsed'`."
    );
  }

  arma::vec u = condition_on_fixed_u
    ? fixed_u_values
    : arma::vec(n * T, arma::fill::zeros);
  arma::vec w_prior_rhs = build_w_prior_canonical_cpp(u, phi, tau2, n, T);
  arma::vec Xbeta = compute_Xbeta_global_cpp(X, beta);
  arma::vec w = w_update(sigma_inv_diag, phi, tau2, y, Xbeta, w_prior_rhs);

  double lCo = 0.0;
  double lCn = 0.0;
  double lCn_1 = 0.0;
  int CL = 1 + (runs - burn - 1) / thin;
  int prediction_draw_count = prediction_enabled ?
    1 + (CL - 1) / prediction_thin : 0;
  int prediction_draw_index = 0;
  int prediction_targets = T * n_prediction;
  arma::fmat prediction_u_draws;
  arma::fmat prediction_w_draws;
  arma::fmat prediction_y_draws;
  arma::fmat prediction_phi_draws;
  arma::fmat prediction_tau2_draws;
  if (prediction_enabled) {
    prediction_u_draws.set_size(prediction_draw_count, prediction_targets);
    prediction_w_draws.set_size(prediction_draw_count, prediction_targets);
    prediction_y_draws.set_size(prediction_draw_count, prediction_targets);
    prediction_phi_draws.set_size(prediction_draw_count, n_prediction);
    prediction_tau2_draws.set_size(prediction_draw_count, n_prediction);
  }
  std::mt19937 prediction_rng(static_cast<std::mt19937::result_type>(
    prediction_seed
  ));
  std::normal_distribution<double> prediction_standard_normal(0.0, 1.0);
  std::uniform_real_distribution<double> prediction_uniform(0.0, 1.0);
  double prediction_seconds = 0.0;

  arma::mat betasample(CL, p, arma::fill::zeros);
  arma::mat w_design_projection(CL, p, arma::fill::zeros);
  arma::mat u_transition_design_projection(CL, p, arma::fill::zeros);
  arma::mat latent_spatial_means(CL, 2, arma::fill::zeros);
  arma::mat design_crossproduct = arma::symmatu(X.t() * X);
  arma::mat wsample;
  arma::mat usample;
  if (save_w) {
    wsample.zeros(CL, n * T);
  }
  if (save_latent) {
    usample.zeros(CL, n * T);
  }
  arma::vec w_sum(n * T, arma::fill::zeros);
  arma::vec u_sum(n * T, arma::fill::zeros);
  arma::mat sigma2sample(CL, 1, arma::fill::zeros);
  arma::mat xi2sample(CL, p, arma::fill::zeros);
  arma::mat phisample(CL, n, arma::fill::zeros);
  arma::mat tau2sample(CL, n, arma::fill::zeros);
  arma::imat clustmemory(CL, n, arma::fill::zeros);
  arma::ivec k_clust(CL, arma::fill::zeros);
  arma::vec M_dp(CL, arma::fill::zeros);
  arma::vec loglik(CL, arma::fill::zeros);
  arma::vec rhosample(CL, arma::fill::zeros);
  arma::vec sigma2_usample(CL, arma::fill::zeros);
  arma::vec y_missing_sum(n_missing, arma::fill::zeros);

  List nj(CL);
  List phistarsample(CL);
  List tau2starsample(CL);

  int acc_phi = 0;
  int acc_tau = 0;
  int acc_rho = 0;
  int prop_phi = 0;
  int prop_tau = 0;
  int prop_rho = 0;
  int gibbs_tau_draws = 0;
  int gibbs_sigma2_u_draws = 0;
  int prop_sigma2 = 0;
  int acc_sigma2 = 0;
  int prop_joint_variance = 0;
  int acc_joint_variance = 0;
  int prop_pairwise_variance = 0;
  int acc_pairwise_variance = 0;
  int prop_spatial = 0;
  int acc_spatial = 0;
  int prop_pairwise_variance_warmup = 0;
  int acc_pairwise_variance_warmup = 0;
  int prop_pairwise_variance_sampling = 0;
  int acc_pairwise_variance_sampling = 0;
  int prop_spatial_warmup = 0;
  int acc_spatial_warmup = 0;
  int prop_spatial_sampling = 0;
  int acc_spatial_sampling = 0;
  int pairwise_window_prop = 0;
  int pairwise_window_acc = 0;
  int spatial_window_prop = 0;
  int spatial_window_acc = 0;
  const double phi_adapt_target = 0.44;
  const double eta_phi_min = 1e-6;
  const double eta_phi_max = 2.0;
  double log_eta_phi = std::log(eta_phi);
  const double rho_adapt_target = 0.44;
  const double eta_rho_min = 1e-6;
  const double eta_rho_max = 2.0;
  double log_eta_rho = std::log(eta_rho);
  const double sigma2_adapt_target = 0.44;
  const double eta_sigma2_min = 1e-6;
  const double eta_sigma2_max = 2.0;
  double log_eta_sigma2 = std::log(eta_sigma2);
  const double tau2_adapt_target = 0.44;
  const double eta_tau2_min = 1e-6;
  const double eta_tau2_max = 2.0;
  double log_eta_tau2 = std::log(eta_tau2);
  std::map<int, RunningCovariance> joint_variance_adaptation;
  std::map<int, RunningCovariance> pairwise_variance_adaptation;
  RunningCovariance spatial_adaptation;
  std::vector<int> adaptation_diagnostic_iteration;
  std::vector<double> adaptation_diagnostic_pairwise_acceptance;
  std::vector<double> adaptation_diagnostic_spatial_acceptance;
  std::vector<double> adaptation_diagnostic_pairwise_scale_min;
  std::vector<double> adaptation_diagnostic_pairwise_scale_median;
  std::vector<double> adaptation_diagnostic_pairwise_scale_max;
  std::vector<double> adaptation_diagnostic_spatial_scale;
  std::vector<double> adaptation_diagnostic_spatial_sd_log_kappa;
  std::vector<double> adaptation_diagnostic_spatial_sd_log_sigma2_u;
  std::deque<int> recent_k_history;
  double recent_k_sum = 0.0;
  std::size_t last_progress_width = 0;
  int next_progress_tenth_percent = 0;
  if (!progress_file.empty()) {
    std::ofstream progress_output(progress_file.c_str(), std::ios::trunc);
    progress_output << 0.0 << std::endl;
  }
  auto chain_start = std::chrono::steady_clock::now();
  double timing_u_seconds = 0.0;
  double timing_beta_seconds = 0.0;
  double timing_w_measurement_seconds = 0.0;
  double timing_allocation_seconds = 0.0;
  double timing_regime_seconds = 0.0;
  double timing_regime_phi_seconds = 0.0;
  double timing_regime_variance_seconds = 0.0;
  double timing_regime_w_redraw_seconds = 0.0;
  double timing_regime_interweave_seconds = 0.0;
  double timing_spatial_seconds = 0.0;
  double timing_storage_diagnostics_seconds = 0.0;
  double timing_iteration_seconds = 0.0;

  for (int iter = 0; iter < runs; ++iter) {
    auto iteration_start = std::chrono::steady_clock::now();
    auto block_start = iteration_start;
    if (!condition_on_fixed_u && !use_collapsed_spatial) {
      u = update_u_cpp(u_sampler, w, phi, tau2, Tauinv_u, n, T);
    }
    timing_u_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - block_start
    ).count();
    w_prior_rhs = build_w_prior_canonical_cpp(u, phi, tau2, n, T);

    update_missing_y_cpp(y, missing_indices, w, Xbeta, sigma2);

    block_start = std::chrono::steady_clock::now();
    beta = beta_update_global_dispatch_cpp(
      use_collapsed_beta, X, n, T,
      xi2_inv, 1.0 / sigma2, phi, tau2, w_prior_rhs, y, w
    );
    Xbeta = compute_Xbeta_global_cpp(X, beta);
    timing_beta_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - block_start
    ).count();

    if (use_pairwise_collapsed_variance) {
      block_start = std::chrono::steady_clock::now();
      std::vector<arma::uvec> pre_w_cluster_members(k_current);
      for (int k = 0; k < k_current; ++k) {
        pre_w_cluster_members[k] = arma::find(clu_current == (k + 1));
      }
      double current_collapsed_loglik =
        collapsed_variance_loglik_parallel_cpp(
          sigma2, tau2_star, clu_current, y_observed, Xbeta, u, phi, n, T
        );
      for (int k = 0; k < k_current; ++k) {
        arma::vec transformed_pair(2);
        transformed_pair(0) = std::log(sigma2);
        transformed_pair(1) = std::log(tau2_star(k) / sigma2);
        int adaptation_key = 0;
        int cluster_size_for_adaptation = static_cast<int>(
          pre_w_cluster_members[k].n_elem
        );
        for (int remaining_size = cluster_size_for_adaptation;
             remaining_size > 1;
             remaining_size /= 2) {
          adaptation_key += 1;
        }
        RunningCovariance& adaptation = pairwise_variance_adaptation[adaptation_key];
        if ((iter + 1) <= burn) {
          adaptation.update(transformed_pair);
        } else if (adaptation.count == 0) {
          adaptation.update(transformed_pair);
        }
        arma::mat proposal_covariance = adaptation.proposal_covariance(
          joint_variance_initial_scale
        );
        arma::mat proposal_cholesky;
        if (!arma::chol(proposal_cholesky, proposal_covariance, "lower")) {
          proposal_covariance.diag() += 1e-6;
          if (!arma::chol(proposal_cholesky, proposal_covariance, "lower")) {
            Rcpp::stop("Pairwise variance proposal covariance is not positive definite.");
          }
        }
        arma::vec proposed_pair = transformed_pair +
          proposal_cholesky * arma::randn<arma::vec>(2);
        double sigma2_prop = std::exp(proposed_pair(0));
        double tau2_prop = std::exp(proposed_pair(0) + proposed_pair(1));
        arma::vec tau2_star_prop = tau2_star;
        tau2_star_prop(k) = tau2_prop;
        double proposed_collapsed_loglik =
          collapsed_variance_loglik_parallel_cpp(
            sigma2_prop, tau2_star_prop, clu_current, y_observed,
            Xbeta, u, phi, n, T
          );
        double log_target_curr = current_collapsed_loglik -
          a_sigma * transformed_pair(0) - b_sigma / sigma2 -
          a_tau * std::log(tau2_star(k)) - b_tau / tau2_star(k);
        double log_target_prop = proposed_collapsed_loglik -
          a_sigma * proposed_pair(0) - b_sigma / sigma2_prop -
          a_tau * std::log(tau2_prop) - b_tau / tau2_prop;
        double log_acceptance = log_target_prop - log_target_curr;
        double alpha_pairwise = std::min(
          1.0, std::exp(std::min(0.0, log_acceptance))
        );
        if (std::log(arma::randu<double>()) < std::min(0.0, log_acceptance)) {
          sigma2 = sigma2_prop;
          tau2_star(k) = tau2_prop;
          tau2.elem(pre_w_cluster_members[k]).fill(tau2_prop);
          current_collapsed_loglik = proposed_collapsed_loglik;
          acc_pairwise_variance += 1;
          pairwise_window_acc += 1;
          if ((iter + 1) <= burn) {
            acc_pairwise_variance_warmup += 1;
          } else {
            acc_pairwise_variance_sampling += 1;
          }
        }
        prop_pairwise_variance += 1;
        pairwise_window_prop += 1;
        if ((iter + 1) <= burn) {
          prop_pairwise_variance_warmup += 1;
        } else {
          prop_pairwise_variance_sampling += 1;
        }
        if ((iter + 1) <= burn) {
          adaptation.adapt_scale(alpha_pairwise, 0.234);
        }
      }
      double pre_w_variance_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - block_start
      ).count();
      timing_regime_variance_seconds += pre_w_variance_seconds;
      timing_regime_seconds += pre_w_variance_seconds;
      sigma_inv_diag.fill(1.0 / sigma2);
      w_prior_rhs = build_w_prior_canonical_cpp(u, phi, tau2, n, T);
    }

    block_start = std::chrono::steady_clock::now();
    if (use_joint_collapsed_variance || use_pairwise_collapsed_variance) {
      w = w_update(sigma_inv_diag, phi, tau2, y, Xbeta, w_prior_rhs);
    } else if (use_collapsed_sigma2) {
      double log_sigma2_curr = std::log(sigma2);
      double log_sigma2_prop = log_sigma2_curr + eta_sigma2 * arma::randn<double>();
      double log_target_curr = collapsed_sigma2_log_target_cpp(
        log_sigma2_curr, y_observed, Xbeta, u, phi, tau2,
        n, T, a_sigma, b_sigma
      );
      double log_target_prop = collapsed_sigma2_log_target_cpp(
        log_sigma2_prop, y_observed, Xbeta, u, phi, tau2,
        n, T, a_sigma, b_sigma
      );
      double log_acceptance = log_target_prop - log_target_curr;
      double alpha_sigma2 = std::min(1.0, std::exp(std::min(0.0, log_acceptance)));
      if (arma::randu<double>() < alpha_sigma2) {
        sigma2 = std::exp(log_sigma2_prop);
        acc_sigma2 += 1;
      }
      prop_sigma2 += 1;
      if ((iter + 1) <= burn) {
        double sigma2_gamma = std::pow(static_cast<double>(prop_sigma2), -0.6);
        log_eta_sigma2 += sigma2_gamma * (alpha_sigma2 - sigma2_adapt_target);
        eta_sigma2 = std::exp(log_eta_sigma2);
        eta_sigma2 = std::max(
          eta_sigma2_min, std::min(eta_sigma2, eta_sigma2_max)
        );
        log_eta_sigma2 = std::log(eta_sigma2);
      }
      sigma_inv_diag.fill(1.0 / sigma2);
      w = w_update(sigma_inv_diag, phi, tau2, y, Xbeta, w_prior_rhs);
    } else {
      w = w_update(sigma_inv_diag, phi, tau2, y, Xbeta, w_prior_rhs);
      sigma2 = sigma2_update_global(y, w, Xbeta, n, T, a_sigma, b_sigma);
      sigma_inv_diag.fill(1.0 / sigma2);
    }

    xi2 = xi2_update_global_cpp(beta, a_xi, b_xi);
    xi2_inv = 1.0 / xi2;
    timing_w_measurement_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - block_start
    ).count();

    block_start = std::chrono::steady_clock::now();
    double debug_alloc_max_sum = 0.0;
    double debug_alloc_second_sum = 0.0;
    TransitionSufficientStatistics transition_statistics =
      build_transition_statistics_cpp(w, u, n, T);
    std::vector<SpatialClusterStatistics> spatial_cluster_statistics =
      build_spatial_cluster_statistics_cpp(
        clu_current,
        k_current,
        s1,
        s2
      );
    const bool use_cached_spatial_cohesion = cohesion == 3 || cohesion == 4;
    if (debug_checks) {
      double direct_transition_loglik = compute_transition_loglik_w_given_u_cpp(
        w, u, phi, tau2, n, T
      );
      double cached_transition_loglik = compute_transition_loglik_cached_cpp(
        transition_statistics, phi, tau2
      );
      double transition_tolerance = 1e-8 *
        std::max(1.0, std::abs(direct_transition_loglik));
      if (std::abs(direct_transition_loglik - cached_transition_loglik) >
          transition_tolerance) {
        Rcpp::stop("Cached transition statistics failed validation.");
      }

      if (use_cached_spatial_cohesion) {
        for (int cluster_index = 0; cluster_index < k_current; ++cluster_index) {
          arma::uvec members = arma::find(clu_current == cluster_index + 1);
          arma::vec cluster_s1 = s1.elem(members);
          arma::vec cluster_s2 = s2.elem(members);
          double direct_cohesion = compute_cohesion_log_value(
            cohesion, cluster_s1, cluster_s2, mu0, k0, v0, L0,
            eC1, aC2, static_cast<int>(members.n_elem)
          );
          double cached_cohesion = cached_cohesion_log_value_cpp(
            cohesion, spatial_cluster_statistics[cluster_index],
            mu0, k0, v0, L0
          );
          double cohesion_tolerance = 1e-8 *
            std::max(1.0, std::abs(direct_cohesion));
          if (std::abs(direct_cohesion - cached_cohesion) >
              cohesion_tolerance) {
            Rcpp::stop("Cached spatial cohesion failed validation.");
          }
        }
      }
    }

    for (int i = 0; i < n; ++i) {
      arma::vec phi_empty = Rcpp::rbeta(CC, a_phi, b_phi);
      arma::vec tau2_empty = 1.0 / arma::randg(
        CC,
        arma::distr_param(a_tau, 1.0 / b_tau)
      );
      SpatialClusterStatistics singleton_statistics =
        empty_spatial_cluster_statistics_cpp();
      update_spatial_cluster_statistics_cpp(
        singleton_statistics,
        s1(i),
        s2(i),
        1
      );
      if (use_cached_spatial_cohesion) {
        lCn_1 = cached_cohesion_log_value_cpp(
          cohesion,
          singleton_statistics,
          mu0,
          k0,
          v0,
          L0
        );
      } else {
        s1o(0) = s1(i);
        s2o(0) = s2(i);
        lCn_1 = 0.0;
      }

      int previous_cluster = clu_current(i);
      int previous_index = previous_cluster - 1;
      update_spatial_cluster_statistics_cpp(
        spatial_cluster_statistics[previous_index],
        s1(i),
        s2(i),
        -1
      );

      if (nj_current(previous_index) > 1) {
        nj_current(previous_index) -= 1;
      } else {
        int retained_auxiliary = Rcpp::sample(CC, 1, true)[0] - 1;
        phi_empty(retained_auxiliary) = phi_star(previous_index);
        tau2_empty(retained_auxiliary) = tau2_star(previous_index);

        nj_current.shed_row(previous_index);
        phi_star.shed_row(previous_index);
        tau2_star.shed_row(previous_index);
        spatial_cluster_statistics.erase(
          spatial_cluster_statistics.begin() + previous_index
        );
        k_current -= 1;

        for (int other_site = 0; other_site < n; ++other_site) {
          if (clu_current(other_site) > previous_cluster) {
            clu_current(other_site) -= 1;
          }
        }
      }

      arma::vec existing_weights(k_current, arma::fill::zeros);
      for (int cluster_index = 0; cluster_index < k_current; ++cluster_index) {
        if (use_cached_spatial_cohesion) {
          lCo = cached_cohesion_log_value_cpp(
            cohesion,
            spatial_cluster_statistics[cluster_index],
            mu0,
            k0,
            v0,
            L0
          );
          SpatialClusterStatistics proposed_statistics =
            spatial_cluster_statistics[cluster_index];
          update_spatial_cluster_statistics_cpp(
            proposed_statistics,
            s1(i),
            s2(i),
            1
          );
          lCn = cached_cohesion_log_value_cpp(
            cohesion,
            proposed_statistics,
            mu0,
            k0,
            v0,
            L0
          );
        } else {
          int member_index = 0;
          for (int other_site = 0; other_site < n; ++other_site) {
            if (other_site != i &&
                clu_current(other_site) == cluster_index + 1) {
              s1o(member_index) = s1(other_site);
              s2o(member_index) = s2(other_site);
              s1n(member_index) = s1(other_site);
              s2n(member_index) = s2(other_site);
              member_index += 1;
            }
          }
          s1n(member_index) = s1(i);
          s2n(member_index) = s2(i);
          lCo = compute_cohesion_log_value(
            cohesion, s1o, s2o, mu0, k0, v0, L0,
            eC1, aC2, member_index
          );
          lCn = compute_cohesion_log_value(
            cohesion, s1n, s2n, mu0, k0, v0, L0,
            eC1, aC2, member_index + 1
          );
        }

        double data_loglik = evaluate_candidate_transition_loglik_cached_cpp(
          transition_statistics,
          i,
          phi_star(cluster_index),
          tau2_star(cluster_index)
        );
        existing_weights(cluster_index) = data_loglik + lCn - lCo +
          std::log(nj_current(cluster_index));
      }

      arma::vec auxiliary_weights(CC, arma::fill::zeros);
      for (int auxiliary = 0; auxiliary < CC; ++auxiliary) {
        double data_loglik = evaluate_candidate_transition_loglik_cached_cpp(
          transition_statistics,
          i,
          phi_empty(auxiliary),
          tau2_empty(auxiliary)
        );
        auxiliary_weights(auxiliary) = data_loglik + lCn_1 +
          std::log(Mdp) - std::log(CC);
      }

      arma::vec weights = arma::join_cols(existing_weights, auxiliary_weights);
      double max_weight = arma::max(weights);
      weights = arma::exp(weights - max_weight);
      weights /= arma::sum(weights);

      if (debug_checks) {
        arma::vec sorted_weights = arma::sort(weights, "descend");
        debug_alloc_max_sum += sorted_weights(0);
        debug_alloc_second_sum += sorted_weights.n_elem > 1 ?
          sorted_weights(1) : 0.0;
      }

      int new_cluster = sample_cpp(k_current + CC, 1, weights, true);
      if (new_cluster <= k_current) {
        clu_current(i) = new_cluster;
        nj_current(new_cluster - 1) += 1;
        update_spatial_cluster_statistics_cpp(
          spatial_cluster_statistics[new_cluster - 1],
          s1(i),
          s2(i),
          1
        );
      } else {
        int selected_auxiliary = new_cluster - k_current - 1;
        clu_current(i) = k_current + 1;
        k_current += 1;
        nj_current.resize(k_current);
        nj_current(k_current - 1) = 1;
        phi_star.resize(k_current);
        phi_star(k_current - 1) = phi_empty(selected_auxiliary);
        tau2_star.resize(k_current);
        tau2_star(k_current - 1) = tau2_empty(selected_auxiliary);
        spatial_cluster_statistics.push_back(singleton_statistics);
      }

      phi(i) = phi_star(clu_current(i) - 1);
      tau2(i) = tau2_star(clu_current(i) - 1);
    }
    timing_allocation_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - block_start
    ).count();

    block_start = std::chrono::steady_clock::now();
    if (debug_checks) {
      if (spatial_cluster_statistics.size() !=
          static_cast<std::size_t>(k_current)) {
        Rcpp::stop("Cached spatial-cluster state has invalid size.");
      }
      for (int i = 0; i < n; ++i) {
        int cluster_idx = clu_current(i) - 1;
        if (cluster_idx < 0 || cluster_idx >= k_current) {
          Rcpp::stop("Invalid cluster index after allocation update.");
        }
        if (std::abs(phi(i) - phi_star(cluster_idx)) > 1e-10 ||
            std::abs(tau2(i) - tau2_star(cluster_idx)) > 1e-10) {
          Rcpp::stop("Cluster-specific phi/tau2 do not match the current allocation.");
        }
      }
      for (int cluster_index = 0; cluster_index < k_current; ++cluster_index) {
        if (spatial_cluster_statistics[cluster_index].count !=
            nj_current(cluster_index)) {
          Rcpp::stop("Cached spatial-cluster count does not match allocation.");
        }
      }
      if (arma::any(phi <= 0.0) || arma::any(phi >= 1.0) || arma::any(tau2 <= 0.0) || sigma2_u <= 0.0) {
        Rcpp::stop("Invalid parameter state detected after allocation update.");
      }
      if (((iter + 1) % 100 == 0) || iter == 0) {
        Rcpp::Rcout << "[debug] iter " << (iter + 1)
                    << " avg max alloc prob=" << (debug_alloc_max_sum / static_cast<double>(n))
                    << ", avg second alloc prob=" << (debug_alloc_second_sum / static_cast<double>(n))
                    << ", global NNGP uses sigma2_u only." << std::endl;
      }
    }

    std::vector<arma::uvec> cluster_members(k_current);
    for (int k = 0; k < k_current; ++k) {
      cluster_members[k] = arma::find(clu_current == (k + 1));
    }

    double transition_loglik_curr = compute_transition_loglik_cached_cpp(
      transition_statistics,
      phi,
      tau2
    );

    auto regime_subblock_start = std::chrono::steady_clock::now();
    for (int k = 0; k < k_current; ++k) {
      const arma::uvec& members = cluster_members[k];
      if (members.n_elem == 0) {
        continue;
      }
      prop_phi += 1;
      double eta_phi_sd = std::sqrt(eta_phi);

      double phi_curr = phi_star(k);
      double trans_phi_curr = arma::trunc_log(phi_curr) - std::log1p(-phi_curr);
      double trans_phi_prop = R::rnorm(trans_phi_curr, eta_phi_sd);
      double phi_prop = arma::trunc_exp(trans_phi_prop) / (1.0 + arma::trunc_exp(trans_phi_prop));

      double cluster_loglik_curr =
        compute_cluster_transition_loglik_cached_cpp(
          transition_statistics,
          members,
          phi_curr,
          tau2_star(k)
        );
      double cluster_loglik_prop =
        compute_cluster_transition_loglik_cached_cpp(
          transition_statistics,
          members,
          phi_prop,
          tau2_star(k)
        );
      double transition_loglik_prop = transition_loglik_curr -
        cluster_loglik_curr + cluster_loglik_prop;

      double logprior_curr =
        (a_phi - 1.0) * arma::trunc_log(phi_curr) +
        (b_phi - 1.0) * std::log1p(-phi_curr) +
        trans_phi_curr - 2.0 * std::log1p(arma::trunc_exp(trans_phi_curr));

      double logprior_prop =
        (a_phi - 1.0) * arma::trunc_log(phi_prop) +
        (b_phi - 1.0) * std::log1p(-phi_prop) +
        trans_phi_prop - 2.0 * std::log1p(arma::trunc_exp(trans_phi_prop));

      double log_alpha = transition_loglik_prop + logprior_prop - transition_loglik_curr - logprior_curr;
      double alpha_phi = std::exp(std::min(0.0, log_alpha));
      if (std::log(arma::randu()) < std::min(0.0, log_alpha)) {
        phi_star(k) = phi_prop;
        phi.elem(members).fill(phi_prop);
        transition_loglik_curr = transition_loglik_prop;
        acc_phi += 1;
      }

      if ((iter + 1) <= burn) {
        double phi_gamma = std::pow(static_cast<double>(prop_phi), -0.6);
        log_eta_phi += phi_gamma * (alpha_phi - phi_adapt_target);
        eta_phi = std::exp(log_eta_phi);
        eta_phi = std::max(eta_phi_min, std::min(eta_phi, eta_phi_max));
        log_eta_phi = std::log(eta_phi);
      }
    }
    timing_regime_phi_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - regime_subblock_start
    ).count();

    regime_subblock_start = std::chrono::steady_clock::now();
    for (int k = 0;
         k < k_current && !use_joint_collapsed_variance &&
           !use_pairwise_collapsed_variance;
         ++k) {
      const arma::uvec& members = cluster_members[k];
      if (members.n_elem == 0) {
        continue;
      }

      if (use_collapsed_tau2) {
        double log_tau2_curr = std::log(tau2_star(k));
        double log_tau2_prop = log_tau2_curr + eta_tau2 * arma::randn<double>();
        double log_target_curr = collapsed_tau2_log_target_cpp(
          log_tau2_curr, members, y_observed, Xbeta, u, phi, sigma2,
          n, T, a_tau, b_tau
        );
        double log_target_prop = collapsed_tau2_log_target_cpp(
          log_tau2_prop, members, y_observed, Xbeta, u, phi, sigma2,
          n, T, a_tau, b_tau
        );
        double log_acceptance = log_target_prop - log_target_curr;
        double alpha_tau2 = std::min(1.0, std::exp(std::min(0.0, log_acceptance)));
        if (arma::randu<double>() < alpha_tau2) {
          double tau2_draw = std::exp(log_tau2_prop);
          tau2_star(k) = tau2_draw;
          tau2.elem(members).fill(tau2_draw);
          acc_tau += 1;
        }
        prop_tau += 1;
        if ((iter + 1) <= burn) {
          double tau2_gamma = std::pow(static_cast<double>(prop_tau), -0.6);
          log_eta_tau2 += tau2_gamma * (alpha_tau2 - tau2_adapt_target);
          eta_tau2 = std::exp(log_eta_tau2);
          eta_tau2 = std::max(eta_tau2_min, std::min(eta_tau2, eta_tau2_max));
          log_eta_tau2 = std::log(eta_tau2);
        }
      } else {
        double phi_k = phi_star(k);
        double sse_k = compute_cluster_tau2_sse_cached_cpp(
          transition_statistics,
          members,
          phi_k
        );
        double shape_post = a_tau + 0.5 * static_cast<double>(members.n_elem * T);
        double scale_post = b_tau + 0.5 * sse_k;
        double tau2_draw = 1.0 / arma::randg<double>(
          arma::distr_param(shape_post, 1.0 / scale_post)
        );
        tau2_star(k) = tau2_draw;
        tau2.elem(members).fill(tau2_draw);
        gibbs_tau_draws += 1;
      }
    }
    timing_regime_variance_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - regime_subblock_start
    ).count();

    if (use_pairwise_collapsed_variance) {
      if (use_pairwise_interweaving) {
        regime_subblock_start = std::chrono::steady_clock::now();
        for (int k = 0; k < k_current; ++k) {
          const arma::uvec& members = cluster_members[k];
          if (members.n_elem == 0) {
            continue;
          }
          double sse_k = compute_cluster_tau2_sse_cpp(
            w, u, members, phi_star(k), n, T
          );
          double shape_post =
            a_tau + 0.5 * static_cast<double>(members.n_elem * T);
          double scale_post = b_tau + 0.5 * sse_k;
          double tau2_draw = 1.0 / arma::randg<double>(
            arma::distr_param(shape_post, 1.0 / scale_post)
          );
          tau2_star(k) = tau2_draw;
          tau2.elem(members).fill(tau2_draw);
          gibbs_tau_draws += 1;
        }
        timing_regime_interweave_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - regime_subblock_start
        ).count();
      }
    } else if (use_joint_collapsed_variance) {
      arma::vec transformed_variances(k_current + 1);
      transformed_variances(0) = std::log(sigma2);
      for (int k = 0; k < k_current; ++k) {
        transformed_variances(k + 1) = std::log(tau2_star(k) / sigma2);
      }

      RunningCovariance& adaptation = joint_variance_adaptation[k_current + 1];
      if ((iter + 1) <= burn) {
        adaptation.update(transformed_variances);
      } else if (adaptation.count == 0) {
        adaptation.update(transformed_variances);
      }
      arma::mat proposal_covariance = adaptation.proposal_covariance(
        joint_variance_initial_scale
      );
      arma::mat proposal_cholesky;
      if (!arma::chol(proposal_cholesky, proposal_covariance, "lower")) {
        proposal_covariance.diag() += 1e-6;
        if (!arma::chol(proposal_cholesky, proposal_covariance, "lower")) {
          Rcpp::stop("Joint variance proposal covariance is not positive definite.");
        }
      }
      arma::vec proposed_variances = transformed_variances +
        proposal_cholesky * arma::randn<arma::vec>(k_current + 1);
      double log_target_curr = joint_collapsed_variance_log_target_cpp(
        transformed_variances, clu_current, y_observed, Xbeta, u, phi,
        n, T, a_sigma, b_sigma, a_tau, b_tau
      );
      double log_target_prop = joint_collapsed_variance_log_target_cpp(
        proposed_variances, clu_current, y_observed, Xbeta, u, phi,
        n, T, a_sigma, b_sigma, a_tau, b_tau
      );
      double log_acceptance = log_target_prop - log_target_curr;
      double alpha_joint_variance = std::min(
        1.0, std::exp(std::min(0.0, log_acceptance))
      );
      if (std::log(arma::randu<double>()) < std::min(0.0, log_acceptance)) {
        sigma2 = std::exp(proposed_variances(0));
        for (int k = 0; k < k_current; ++k) {
          double tau2_draw = std::exp(
            proposed_variances(0) + proposed_variances(k + 1)
          );
          tau2_star(k) = tau2_draw;
          tau2.elem(cluster_members[k]).fill(tau2_draw);
        }
        acc_joint_variance += 1;
      }
      prop_joint_variance += 1;
      if ((iter + 1) <= burn) {
        adaptation.adapt_scale(alpha_joint_variance, 0.234);
      }
      sigma_inv_diag.fill(1.0 / sigma2);
      w_prior_rhs = build_w_prior_canonical_cpp(u, phi, tau2, n, T);
      w = w_update(sigma_inv_diag, phi, tau2, y, Xbeta, w_prior_rhs);
    } else if (use_collapsed_tau2) {
      w_prior_rhs = build_w_prior_canonical_cpp(u, phi, tau2, n, T);
      w = w_update(sigma_inv_diag, phi, tau2, y, Xbeta, w_prior_rhs);
    }
    timing_regime_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - block_start
    ).count();

    block_start = std::chrono::steady_clock::now();
    if (use_collapsed_spatial) {
      arma::vec transformed_spatial(2);
      transformed_spatial(0) = std::log(rho);
      transformed_spatial(1) = std::log(sigma2_u);
      if ((iter + 1) <= burn || spatial_adaptation.count == 0) {
        spatial_adaptation.update(transformed_spatial);
      }
      arma::mat proposal_covariance = spatial_adaptation.proposal_covariance(
        spatial_initial_scale,
        1e-10
      );
      arma::mat proposal_cholesky;
      if (!arma::chol(proposal_cholesky, proposal_covariance, "lower")) {
        proposal_covariance.diag() += 1e-6;
        if (!arma::chol(proposal_cholesky, proposal_covariance, "lower")) {
          Rcpp::stop("Collapsed spatial proposal covariance is not positive definite.");
        }
      }
      arma::vec proposed_spatial = transformed_spatial +
        proposal_cholesky * arma::randn<arma::vec>(2);
      double rho_prop = std::exp(proposed_spatial(0));
      double sigma2_u_prop = std::exp(proposed_spatial(1));
      arma::vec F_u_prop;
      arma::vec values_u_prop;
      arma::sp_mat Tauinv_u_prop = build_spatial_precision(
        sigma2_u_prop, rho_prop, F_u_prop, values_u_prop
      );
      double log_target_curr = collapsed_spatial_log_target_from_state_cpp(
        transformed_spatial(0), transformed_spatial(1), w, phi, tau2,
        Tauinv_u, F_u, collapsed_spatial_sampler, n, T,
        a_rho, b_rho, a_sigma2_u, b_sigma2_u
      );
      double log_target_prop = collapsed_spatial_log_target_from_state_cpp(
        proposed_spatial(0), proposed_spatial(1), w, phi, tau2,
        Tauinv_u_prop, F_u_prop, collapsed_spatial_sampler, n, T,
        a_rho, b_rho, a_sigma2_u, b_sigma2_u
      );
      double log_alpha_spatial = log_target_prop - log_target_curr;
      double alpha_spatial = std::exp(std::min(0.0, log_alpha_spatial));
      if (std::log(arma::randu<double>()) < std::min(0.0, log_alpha_spatial)) {
        rho = rho_prop;
        sigma2_u = sigma2_u_prop;
        F_u = F_u_prop;
        values_u = values_u_prop;
        Tauinv_u = Tauinv_u_prop;
        acc_spatial += 1;
        acc_rho += 1;
        spatial_window_acc += 1;
        if ((iter + 1) <= burn) {
          acc_spatial_warmup += 1;
        } else {
          acc_spatial_sampling += 1;
        }
      }
      prop_spatial += 1;
      prop_rho += 1;
      spatial_window_prop += 1;
      if ((iter + 1) <= burn) {
        prop_spatial_warmup += 1;
      } else {
        prop_spatial_sampling += 1;
      }
      if ((iter + 1) <= burn) {
        spatial_adaptation.adapt_scale(alpha_spatial, 0.234, -10.0, 5.0);
      }
      timing_spatial_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - block_start
      ).count();
      block_start = std::chrono::steady_clock::now();
      u = update_u_cpp(u_sampler, w, phi, tau2, Tauinv_u, n, T);
      timing_u_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - block_start
      ).count();
    } else {
      double loglik_u_curr = compute_nngp_loglik_u_cpp(u, F_u, values_u, precomp, n, T);

      prop_rho += 1;
      double eta_rho_sd = std::sqrt(eta_rho);
      double log_rho_curr = std::log(rho);
      double log_rho_prop = R::rnorm(log_rho_curr, eta_rho_sd);
      double rho_prop = std::exp(log_rho_prop);

      arma::vec F_u_prop;
      arma::vec values_u_prop;
      arma::sp_mat Tauinv_u_prop = build_spatial_precision(
        sigma2_u, rho_prop, F_u_prop, values_u_prop
      );
      double loglik_u_prop = compute_nngp_loglik_u_cpp(u, F_u_prop, values_u_prop, precomp, n, T);

      double logprior_rho_curr =
        (a_rho - 1.0) * log_rho_curr -
        b_rho * rho +
        log_rho_curr;

      double logprior_rho_prop =
        (a_rho - 1.0) * log_rho_prop -
        b_rho * rho_prop +
        log_rho_prop;

      double log_alpha_rho = loglik_u_prop + logprior_rho_prop - loglik_u_curr - logprior_rho_curr;
      double alpha_rho = std::exp(std::min(0.0, log_alpha_rho));

      if (std::log(arma::randu()) < std::min(0.0, log_alpha_rho)) {
        rho = rho_prop;
        F_u = F_u_prop;
        values_u = values_u_prop;
        Tauinv_u = Tauinv_u_prop;
        acc_rho += 1;
      }

      if ((iter + 1) <= burn) {
        double rho_gamma = std::pow(static_cast<double>(prop_rho), -0.6);
        log_eta_rho += rho_gamma * (alpha_rho - rho_adapt_target);
        eta_rho = std::exp(log_eta_rho);
        eta_rho = std::max(eta_rho_min, std::min(eta_rho, eta_rho_max));
        log_eta_rho = std::log(eta_rho);
      }

      double sigma2_u_sse = compute_sigma2_u_sse_cpp(u, F_u, values_u, precomp, n, T, sigma2_u);
      double sigma2_u_shape = a_sigma2_u + 0.5 * static_cast<double>(n * std::max(T - 1, 0));
      double sigma2_u_scale = b_sigma2_u + 0.5 * sigma2_u_sse;
      sigma2_u = 1.0 / arma::randg<double>(arma::distr_param(sigma2_u_shape, 1.0 / sigma2_u_scale));
      Tauinv_u = build_spatial_precision(sigma2_u, rho, F_u, values_u);
      gibbs_sigma2_u_draws += 1;
      timing_spatial_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - block_start
      ).count();
    }

    block_start = std::chrono::steady_clock::now();
    if (debug_checks && (arma::any(phi <= 0.0) || arma::any(phi >= 1.0) ||
        arma::any(tau2 <= 0.0) || sigma2_u <= 0.0)) {
      Rcpp::stop("Invalid phi/tau2/sigma2_u values after parameter updates.");
    }

    if ((iter >= burn) && (((iter - burn) % thin) == 0)) {
      int save_idx = (iter - burn) / thin;
      betasample.row(save_idx) = beta.t();
      arma::vec current_w_projection;
      arma::vec current_u_projection;
      compute_latent_design_projections_cpp(
        X,
        design_crossproduct,
        w,
        u,
        phi,
        n,
        T,
        current_w_projection,
        current_u_projection
      );
      w_design_projection.row(save_idx) = current_w_projection.t();
      u_transition_design_projection.row(save_idx) =
        current_u_projection.t();
      latent_spatial_means(save_idx, 0) = arma::mean(w);
      latent_spatial_means(save_idx, 1) =
        T > 1 ? arma::mean(u.tail(n * (T - 1))) : 0.0;
      if (save_w) {
        wsample.row(save_idx) = w.t();
      }
      if (save_latent) {
        usample.row(save_idx) = u.t();
      }
      w_sum += w;
      u_sum += u;
      sigma2sample(save_idx, 0) = sigma2;
      xi2sample.row(save_idx) = xi2.t();
      phisample.row(save_idx) = phi.t();
      tau2sample.row(save_idx) = tau2.t();
      clustmemory.row(save_idx) = clu_current.t();
      phistarsample[save_idx] = phi_star;
      tau2starsample[save_idx] = tau2_star;
      nj[save_idx] = nj_current;
      k_clust(save_idx) = k_current;
      M_dp(save_idx) = Mdp;
      rhosample(save_idx) = rho;
      sigma2_usample(save_idx) = sigma2_u;
      loglik(save_idx) = observed_gaussian_loglik_cpp(
        y_observed, observed_indices, w, Xbeta, sigma2
      );
      if (n_missing > 0) {
        y_missing_sum += y.elem(missing_indices);
      }

      if (prediction_enabled && (save_idx % prediction_thin) == 0) {
        auto prediction_start = std::chrono::steady_clock::now();
        arma::mat u_observed(
          const_cast<double*>(u.memptr()), n, T, false, false
        );
        std::vector<SpatialClusterStatistics> prediction_cluster_statistics =
          build_spatial_cluster_statistics_cpp(
            clu_current,
            k_current,
            s1,
            s2
          );
        std::gamma_distribution<double> phi_numerator_distribution(a_phi, 1.0);
        std::gamma_distribution<double> phi_denominator_distribution(b_phi, 1.0);
        std::gamma_distribution<double> tau_precision_distribution(
          a_tau,
          1.0 / b_tau
        );

        for (int prediction = 0; prediction < n_prediction; ++prediction) {
          double prediction_s1 = prediction_similarity_matrix(prediction, 1);
          double prediction_s2 = prediction_similarity_matrix(prediction, 0);
          SpatialClusterStatistics singleton =
            empty_spatial_cluster_statistics_cpp();
          update_spatial_cluster_statistics_cpp(
            singleton,
            prediction_s1,
            prediction_s2,
            1
          );
          double singleton_log_cohesion = cached_cohesion_log_value_cpp(
            cohesion,
            singleton,
            mu0,
            k0,
            v0,
            L0
          );
          arma::vec log_weights(k_current + 1, arma::fill::zeros);
          for (int cluster_index = 0;
               cluster_index < k_current;
               ++cluster_index) {
            double old_log_cohesion = cached_cohesion_log_value_cpp(
              cohesion,
              prediction_cluster_statistics[cluster_index],
              mu0,
              k0,
              v0,
              L0
            );
            SpatialClusterStatistics proposed =
              prediction_cluster_statistics[cluster_index];
            update_spatial_cluster_statistics_cpp(
              proposed,
              prediction_s1,
              prediction_s2,
              1
            );
            double new_log_cohesion = cached_cohesion_log_value_cpp(
              cohesion,
              proposed,
              mu0,
              k0,
              v0,
              L0
            );
            log_weights(cluster_index) =
              std::log(static_cast<double>(nj_current(cluster_index))) +
              new_log_cohesion - old_log_cohesion;
          }
          log_weights(k_current) = std::log(Mdp) + singleton_log_cohesion;
          double maximum_weight = log_weights.max();
          arma::vec probabilities = arma::exp(log_weights - maximum_weight);
          probabilities /= arma::sum(probabilities);
          double draw_uniform = prediction_uniform(prediction_rng);
          double cumulative = 0.0;
          int selected_cluster = k_current;
          for (int candidate = 0; candidate < k_current + 1; ++candidate) {
            cumulative += probabilities(candidate);
            if (draw_uniform <= cumulative) {
              selected_cluster = candidate;
              break;
            }
          }

          double prediction_phi;
          double prediction_tau2;
          if (selected_cluster < k_current) {
            prediction_phi = phi_star(selected_cluster);
            prediction_tau2 = tau2_star(selected_cluster);
          } else {
            double first_gamma = phi_numerator_distribution(prediction_rng);
            double second_gamma = phi_denominator_distribution(prediction_rng);
            prediction_phi = first_gamma / (first_gamma + second_gamma);
            prediction_tau2 = 1.0 / tau_precision_distribution(prediction_rng);
          }
          prediction_phi_draws(prediction_draw_index, prediction) =
            static_cast<float>(prediction_phi);
          prediction_tau2_draws(prediction_draw_index, prediction) =
            static_cast<float>(prediction_tau2);

          arma::vec prediction_coefficients;
          double prediction_variance = 0.0;
          prediction_nngp_factor_cpp(
            prediction_geometry[prediction],
            sigma2_u,
            rho,
            nu,
            prediction_coefficients,
            prediction_variance
          );

          double predicted_u = 0.0;
          double predicted_w = std::sqrt(
            prediction_tau2 / (1.0 - prediction_phi * prediction_phi)
          ) * prediction_standard_normal(prediction_rng);
          for (int time = 0; time < T; ++time) {
            if (time > 0) {
              double conditional_mean = 0.0;
              for (arma::uword neighbor = 0;
                   neighbor < prediction_coefficients.n_elem;
                   ++neighbor) {
                conditional_mean += prediction_coefficients(neighbor) *
                  u_observed(
                    prediction_geometry[prediction].neighbors(neighbor),
                    time
                  );
              }
              predicted_u = conditional_mean + std::sqrt(prediction_variance) *
                prediction_standard_normal(prediction_rng);
              predicted_w = prediction_phi * predicted_w + predicted_u +
                std::sqrt(prediction_tau2) *
                  prediction_standard_normal(prediction_rng);
            }
            int target = time * n_prediction + prediction;
            double predicted_mean = arma::dot(
              prediction_design.row(target),
              beta.t()
            ) + predicted_w;
            double predicted_y = predicted_mean + std::sqrt(sigma2) *
              prediction_standard_normal(prediction_rng);
            prediction_u_draws(prediction_draw_index, target) =
              static_cast<float>(predicted_u);
            prediction_w_draws(prediction_draw_index, target) =
              static_cast<float>(predicted_w);
            prediction_y_draws(prediction_draw_index, target) =
              static_cast<float>(predicted_y);
          }
        }
        prediction_draw_index += 1;
        prediction_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - prediction_start
        ).count();
      }
    }

    if ((iter + 1) % 100 == 0 || iter == runs - 1) {
      std::vector<double> pairwise_scales;
      for (const auto& adaptation_entry : pairwise_variance_adaptation) {
        pairwise_scales.push_back(
          joint_variance_initial_scale *
            std::exp(adaptation_entry.second.log_scale)
        );
      }
      std::sort(pairwise_scales.begin(), pairwise_scales.end());
      double pairwise_scale_min = NA_REAL;
      double pairwise_scale_median = NA_REAL;
      double pairwise_scale_max = NA_REAL;
      if (!pairwise_scales.empty()) {
        pairwise_scale_min = pairwise_scales.front();
        pairwise_scale_max = pairwise_scales.back();
        std::size_t middle = pairwise_scales.size() / 2;
        pairwise_scale_median = pairwise_scales.size() % 2 == 0 ?
          0.5 * (pairwise_scales[middle - 1] + pairwise_scales[middle]) :
          pairwise_scales[middle];
      }

      double spatial_scale = NA_REAL;
      double spatial_sd_log_kappa = NA_REAL;
      double spatial_sd_log_sigma2_u = NA_REAL;
      if (spatial_adaptation.count > 0) {
        spatial_scale = spatial_initial_scale *
          std::exp(spatial_adaptation.log_scale);
        arma::mat spatial_covariance =
          spatial_adaptation.proposal_covariance(
            spatial_initial_scale,
            1e-10
          );
        spatial_sd_log_kappa = std::sqrt(spatial_covariance(0, 0));
        spatial_sd_log_sigma2_u = std::sqrt(spatial_covariance(1, 1));
      }

      adaptation_diagnostic_iteration.push_back(iter + 1);
      adaptation_diagnostic_pairwise_acceptance.push_back(
        pairwise_window_prop > 0 ?
          static_cast<double>(pairwise_window_acc) /
            static_cast<double>(pairwise_window_prop) : NA_REAL
      );
      adaptation_diagnostic_spatial_acceptance.push_back(
        spatial_window_prop > 0 ?
          static_cast<double>(spatial_window_acc) /
            static_cast<double>(spatial_window_prop) : NA_REAL
      );
      adaptation_diagnostic_pairwise_scale_min.push_back(pairwise_scale_min);
      adaptation_diagnostic_pairwise_scale_median.push_back(
        pairwise_scale_median
      );
      adaptation_diagnostic_pairwise_scale_max.push_back(pairwise_scale_max);
      adaptation_diagnostic_spatial_scale.push_back(spatial_scale);
      adaptation_diagnostic_spatial_sd_log_kappa.push_back(
        spatial_sd_log_kappa
      );
      adaptation_diagnostic_spatial_sd_log_sigma2_u.push_back(
        spatial_sd_log_sigma2_u
      );
      pairwise_window_prop = 0;
      pairwise_window_acc = 0;
      spatial_window_prop = 0;
      spatial_window_acc = 0;
    }

    recent_k_history.push_back(k_current);
    recent_k_sum += static_cast<double>(k_current);
    if (recent_k_history.size() > 100) {
      recent_k_sum -= static_cast<double>(recent_k_history.front());
      recent_k_history.pop_front();
    }

    double progress = 100.0 * static_cast<double>(iter + 1) / static_cast<double>(runs);
    int completed_tenth_percent = static_cast<int>(std::floor(10.0 * progress));
    if (iter == runs - 1) {
      completed_tenth_percent = 1000;
    }

    if (completed_tenth_percent >= next_progress_tenth_percent) {
      double elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - chain_start).count();
      double average_iter_seconds = elapsed_seconds / static_cast<double>(iter + 1);
      double remaining_seconds = average_iter_seconds * static_cast<double>(runs - iter - 1);
      double average_recent_k = recent_k_sum / static_cast<double>(recent_k_history.size());
      double rho_acceptance_rate =
        prop_rho > 0 ? static_cast<double>(acc_rho) / static_cast<double>(prop_rho) : 0.0;
      std::ostringstream progress_stream;
      progress_stream << std::fixed << std::setprecision(1)
                      << "Progress: " << progress
                      << "% | elapsed: " << elapsed_seconds
                      << "s | eta: " << remaining_seconds
                      << "s | avg_k_last100: " << average_recent_k
                      << " | acc_kappa: " << std::setprecision(3) << rho_acceptance_rate
                      << " | sigma2_u: "
                      << (use_collapsed_spatial ? "collapsed" : "gibbs");
      std::string progress_line = progress_stream.str();
      if (progress_line.size() < last_progress_width) {
        progress_line.append(last_progress_width - progress_line.size(), ' ');
      }
      last_progress_width = progress_line.size();

      if (!progress_file.empty()) {
        std::ofstream progress_output(progress_file.c_str(), std::ios::trunc);
        progress_output << std::fixed << std::setprecision(1)
                        << progress << std::endl;
      }
      if (show_progress) {
        Rprintf("\r%s", progress_line.c_str());
        R_FlushConsole();
      }
      next_progress_tenth_percent = completed_tenth_percent + 1;
    }
    timing_storage_diagnostics_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - block_start
    ).count();
    timing_iteration_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - iteration_start
    ).count();
  }

  if (show_progress) {
    Rprintf("\n");
    R_FlushConsole();
  }

  List out;
  out["beta"] = betasample;
  out["w_design_projection"] = w_design_projection;
  out["u_transition_design_projection"] =
    u_transition_design_projection;
  out["latent_spatial_means"] = latent_spatial_means;
  if (save_w) {
    out["w"] = wsample;
  } else {
    out["w"] = R_NilValue;
  }
  if (save_latent) {
    out["u"] = usample;
  } else {
    out["u"] = R_NilValue;
  }
  if (CL > 0) {
    out["w_mean"] = arma::reshape(w_sum / static_cast<double>(CL), n, T).t();
    out["u_mean"] = arma::reshape(u_sum / static_cast<double>(CL), n, T).t();
  } else {
    out["w_mean"] = arma::mat(T, n, arma::fill::value(NA_REAL));
    out["u_mean"] = arma::mat(T, n, arma::fill::value(NA_REAL));
  }
  out["save_latent"] = save_latent;
  out["save_u"] = save_latent;
  out["save_w"] = save_w;
  out["condition_on_fixed_u"] = condition_on_fixed_u;
  out["optimized_updates"] = CharacterVector::create(
    "sparse_parallel_timewise_u",
    "tridiagonal_sitewise_w",
    "parallel_tridiagonal_collapsed_global_beta",
    "cached_parallel_pairwise_variance_likelihood",
    "cached_cluster_membership_lists",
    "single_w_draw_pairwise_reordering",
    "cached_transition_statistics",
    "cached_spatial_cohesion_3_4"
  );
  const double measured_block_seconds =
    timing_u_seconds + timing_beta_seconds + timing_w_measurement_seconds +
    timing_allocation_seconds + timing_regime_seconds + timing_spatial_seconds +
    timing_storage_diagnostics_seconds;
  const double timing_other_seconds = std::max(
    0.0,
    timing_iteration_seconds - measured_block_seconds
  );
  CharacterVector timing_block_names = CharacterVector::create(
    "u", "beta", "w_measurement", "allocation", "regime", "spatial",
    "storage_diagnostics", "other"
  );
  NumericVector timing_block_seconds = NumericVector::create(
    timing_u_seconds,
    timing_beta_seconds,
    timing_w_measurement_seconds,
    timing_allocation_seconds,
    timing_regime_seconds,
    timing_spatial_seconds,
    timing_storage_diagnostics_seconds,
    timing_other_seconds
  );
  NumericVector timing_block_percent(timing_block_seconds.size());
  NumericVector timing_block_seconds_per_iteration(timing_block_seconds.size());
  for (R_xlen_t timing_index = 0;
       timing_index < timing_block_seconds.size();
       ++timing_index) {
    timing_block_percent[timing_index] = timing_iteration_seconds > 0.0 ?
      100.0 * timing_block_seconds[timing_index] / timing_iteration_seconds :
      NA_REAL;
    timing_block_seconds_per_iteration[timing_index] = runs > 0 ?
      timing_block_seconds[timing_index] / static_cast<double>(runs) :
      NA_REAL;
  }
  out["block_timing"] = DataFrame::create(
    Named("block") = timing_block_names,
    Named("seconds") = timing_block_seconds,
    Named("percent") = timing_block_percent,
    Named("seconds_per_iteration") = timing_block_seconds_per_iteration,
    Named("stringsAsFactors") = false
  );
  const double measured_regime_seconds =
    timing_regime_phi_seconds + timing_regime_variance_seconds +
    timing_regime_w_redraw_seconds + timing_regime_interweave_seconds;
  const double timing_regime_other_seconds = std::max(
    0.0,
    timing_regime_seconds - measured_regime_seconds
  );
  CharacterVector regime_timing_names = CharacterVector::create(
    "phi", "variance", "w_redraw", "interweave", "other"
  );
  NumericVector regime_timing_seconds = NumericVector::create(
    timing_regime_phi_seconds,
    timing_regime_variance_seconds,
    timing_regime_w_redraw_seconds,
    timing_regime_interweave_seconds,
    timing_regime_other_seconds
  );
  NumericVector regime_timing_percent(regime_timing_seconds.size());
  NumericVector regime_timing_total_percent(regime_timing_seconds.size());
  for (R_xlen_t timing_index = 0;
       timing_index < regime_timing_seconds.size();
       ++timing_index) {
    regime_timing_percent[timing_index] = timing_regime_seconds > 0.0 ?
      100.0 * regime_timing_seconds[timing_index] / timing_regime_seconds :
      NA_REAL;
    regime_timing_total_percent[timing_index] = timing_iteration_seconds > 0.0 ?
      100.0 * regime_timing_seconds[timing_index] / timing_iteration_seconds :
      NA_REAL;
  }
  out["regime_block_timing"] = DataFrame::create(
    Named("block") = regime_timing_names,
    Named("seconds") = regime_timing_seconds,
    Named("percent_of_regime") = regime_timing_percent,
    Named("percent_of_total") = regime_timing_total_percent,
    Named("stringsAsFactors") = false
  );
  out["timed_iterations"] = runs;
  out["timed_total_seconds"] = timing_iteration_seconds;
#ifdef _OPENMP
  out["openmp_threads"] = omp_get_max_threads();
#else
  out["openmp_threads"] = 1;
#endif
  out["sigma2"] = sigma2sample;
  out["sigma2_update"] = use_pairwise_interweaving ?
    "pairwise_interwoven_mh" : (use_pairwise_collapsed_variance ?
    "pairwise_collapsed_mh" : (use_joint_collapsed_variance ?
      "joint_collapsed_mh" : (use_collapsed_sigma2 ?
        "collapsed_mh" : "conditional_gibbs")));
  out["accsigma2"] = prop_sigma2 > 0 ?
    static_cast<double>(acc_sigma2) / static_cast<double>(prop_sigma2) : NA_REAL;
  out["n_sigma2_prop"] = prop_sigma2;
  out["n_sigma2_acc"] = acc_sigma2;
  out["eta_sigma2_final"] = eta_sigma2;
  out["acc_joint_variance"] = prop_joint_variance > 0 ?
    static_cast<double>(acc_joint_variance) /
      static_cast<double>(prop_joint_variance) : NA_REAL;
  out["n_joint_variance_prop"] = prop_joint_variance;
  out["n_joint_variance_acc"] = acc_joint_variance;
  out["acc_pairwise_variance"] = prop_pairwise_variance > 0 ?
    static_cast<double>(acc_pairwise_variance) /
      static_cast<double>(prop_pairwise_variance) : NA_REAL;
  out["n_pairwise_variance_prop"] = prop_pairwise_variance;
  out["n_pairwise_variance_acc"] = acc_pairwise_variance;
  out["acc_pairwise_variance_warmup"] =
    prop_pairwise_variance_warmup > 0 ?
      static_cast<double>(acc_pairwise_variance_warmup) /
        static_cast<double>(prop_pairwise_variance_warmup) : NA_REAL;
  out["acc_pairwise_variance_sampling"] =
    prop_pairwise_variance_sampling > 0 ?
      static_cast<double>(acc_pairwise_variance_sampling) /
        static_cast<double>(prop_pairwise_variance_sampling) : NA_REAL;
  out["n_pairwise_variance_prop_warmup"] =
    prop_pairwise_variance_warmup;
  out["n_pairwise_variance_prop_sampling"] =
    prop_pairwise_variance_sampling;

  IntegerVector pairwise_size_bin;
  IntegerVector pairwise_size_min;
  IntegerVector pairwise_size_max;
  IntegerVector pairwise_adaptation_count;
  IntegerVector pairwise_adaptation_proposals;
  NumericVector pairwise_adaptation_scale;
  NumericVector pairwise_adaptation_sd_log_sigma2;
  NumericVector pairwise_adaptation_sd_log_ratio;
  for (const auto& adaptation_entry : pairwise_variance_adaptation) {
    int size_bin = adaptation_entry.first;
    const RunningCovariance& adaptation = adaptation_entry.second;
    arma::mat proposal_covariance = adaptation.proposal_covariance(
      joint_variance_initial_scale
    );
    pairwise_size_bin.push_back(size_bin);
    pairwise_size_min.push_back(1 << size_bin);
    pairwise_size_max.push_back((1 << (size_bin + 1)) - 1);
    pairwise_adaptation_count.push_back(adaptation.count);
    pairwise_adaptation_proposals.push_back(adaptation.proposals);
    pairwise_adaptation_scale.push_back(
      joint_variance_initial_scale * std::exp(adaptation.log_scale)
    );
    pairwise_adaptation_sd_log_sigma2.push_back(
      std::sqrt(proposal_covariance(0, 0))
    );
    pairwise_adaptation_sd_log_ratio.push_back(
      std::sqrt(proposal_covariance(1, 1))
    );
  }
  out["pairwise_adaptation_by_cluster_size"] = DataFrame::create(
    Named("size_bin") = pairwise_size_bin,
    Named("cluster_size_min") = pairwise_size_min,
    Named("cluster_size_max") = pairwise_size_max,
    Named("covariance_updates") = pairwise_adaptation_count,
    Named("scale_updates") = pairwise_adaptation_proposals,
    Named("effective_scale") = pairwise_adaptation_scale,
    Named("proposal_sd_log_sigma2") = pairwise_adaptation_sd_log_sigma2,
    Named("proposal_sd_log_tau2_over_sigma2") =
      pairwise_adaptation_sd_log_ratio,
    Named("stringsAsFactors") = false
  );
  out["sigma2_u"] = sigma2_usample;
  out["sigma2beta"] = xi2sample;
  out["beta_parameterization"] = "global";
  out["sampler_variant"] = "real_data_2_confounding_diagnostics";
  out["latent_projection_definition"] =
    "w ~ X and u_t ~ X_t - diag(phi) X_{t-1}";
  out["design_order"] = "time-major: row t*n+i contains x_t(s_i)";
  out["phi"] = phisample;
  out["tau2"] = tau2sample;
  out["cluster"] = clustmemory;
  out["nj"] = nj;
  out["phistar"] = phistarsample;
  out["accphi"] = prop_phi > 0 ? acc_phi / static_cast<double>(prop_phi) : NA_REAL;
  out["n_phi_prop"] = prop_phi;
  out["n_phi_acc"] = acc_phi;
  out["eta_phi_final"] = eta_phi;
  out["tau2star"] = tau2starsample;
  out["acctau"] = prop_tau > 0 ?
    static_cast<double>(acc_tau) / static_cast<double>(prop_tau) : NA_REAL;
  out["n_tau_prop"] = prop_tau;
  out["n_tau_acc"] = acc_tau;
  out["tau2_update"] = use_pairwise_interweaving ?
    "pairwise_collapsed_plus_gibbs" : (use_pairwise_collapsed_variance ?
    "pairwise_collapsed_mh" : (use_joint_collapsed_variance ?
      "joint_collapsed_mh" : (use_collapsed_tau2 ?
        "collapsed_mh" : "conditional_gibbs")));
  out["n_tau_gibbs"] = gibbs_tau_draws;
  out["eta_tau2_final"] = eta_tau2;
  out["kappa"] = rhosample;
  out["acckappa"] = prop_rho > 0 ? acc_rho / static_cast<double>(prop_rho) : NA_REAL;
  out["n_kappa_prop"] = prop_rho;
  out["n_kappa_acc"] = acc_rho;
  out["spatial_update"] = use_collapsed_spatial ?
    "joint_collapsed_mh" : "conditional_kappa_plus_sigma2_u_gibbs";
  out["acc_spatial"] = prop_spatial > 0 ?
    static_cast<double>(acc_spatial) / static_cast<double>(prop_spatial) : NA_REAL;
  out["n_spatial_prop"] = prop_spatial;
  out["n_spatial_acc"] = acc_spatial;
  out["acc_spatial_warmup"] = prop_spatial_warmup > 0 ?
    static_cast<double>(acc_spatial_warmup) /
      static_cast<double>(prop_spatial_warmup) : NA_REAL;
  out["acc_spatial_sampling"] = prop_spatial_sampling > 0 ?
    static_cast<double>(acc_spatial_sampling) /
      static_cast<double>(prop_spatial_sampling) : NA_REAL;
  out["n_spatial_prop_warmup"] = prop_spatial_warmup;
  out["n_spatial_prop_sampling"] = prop_spatial_sampling;
  out["sigma2_u_update"] = use_collapsed_spatial ? "joint_collapsed_mh" : "gibbs";
  out["spatial_process_model"] = spatial_process_model;
  out["n_sigma2_u_gibbs"] = gibbs_sigma2_u_draws;
  out["eta_kappa_final"] = eta_rho;
  out["spatial_proposal_scale_final"] =
    spatial_initial_scale * std::exp(spatial_adaptation.log_scale);
  out["adaptation_diagnostics"] = DataFrame::create(
    Named("iteration") = wrap(adaptation_diagnostic_iteration),
    Named("pairwise_acceptance_last_100") =
      wrap(adaptation_diagnostic_pairwise_acceptance),
    Named("spatial_acceptance_last_100") =
      wrap(adaptation_diagnostic_spatial_acceptance),
    Named("pairwise_scale_min") =
      wrap(adaptation_diagnostic_pairwise_scale_min),
    Named("pairwise_scale_median") =
      wrap(adaptation_diagnostic_pairwise_scale_median),
    Named("pairwise_scale_max") =
      wrap(adaptation_diagnostic_pairwise_scale_max),
    Named("spatial_scale") = wrap(adaptation_diagnostic_spatial_scale),
    Named("spatial_sd_log_kappa") =
      wrap(adaptation_diagnostic_spatial_sd_log_kappa),
    Named("spatial_sd_log_sigma2_u") =
      wrap(adaptation_diagnostic_spatial_sd_log_sigma2_u),
    Named("stringsAsFactors") = false
  );
  out["k_clust"] = k_clust;
  out["Mdp"] = M_dp;
  out["loglik"] = loglik;
  out["n_missing"] = static_cast<double>(n_missing);
  if (n_missing > 0 && CL > 0) {
    arma::vec y_imputed_mean = y_observed;
    y_imputed_mean.elem(missing_indices) =
      y_missing_sum / static_cast<double>(CL);
    out["Y_imputed_mean"] = arma::reshape(y_imputed_mean, n, T).t();
  } else {
    out["Y_imputed_mean"] = Y;
  }
  out["prediction_enabled"] = prediction_enabled;
  if (prediction_enabled) {
    if (prediction_draw_index != prediction_draw_count) {
      Rcpp::stop("Internal prediction-draw count mismatch.");
    }
    List prediction_output = List::create(
      Named("u") = summarize_prediction_draws_cpp(
        prediction_u_draws,
        T,
        n_prediction,
        prediction_truth_u
      ),
      Named("w") = summarize_prediction_draws_cpp(
        prediction_w_draws,
        T,
        n_prediction,
        prediction_truth_w
      ),
      Named("y") = summarize_prediction_draws_cpp(
        prediction_y_draws,
        T,
        n_prediction,
        prediction_truth_y
      )
    );

    bool has_cluster_truth = !Rf_isNull(prediction_true_cluster) &&
      !Rf_isNull(prediction_cluster_phi) &&
      !Rf_isNull(prediction_cluster_tau2);
    arma::mat phi_truth_matrix;
    arma::mat tau2_truth_matrix;
    arma::ivec true_prediction_cluster;
    arma::vec true_cluster_phi;
    arma::vec true_cluster_tau2;
    if (has_cluster_truth) {
      true_prediction_cluster = Rcpp::as<arma::ivec>(prediction_true_cluster);
      true_cluster_phi = Rcpp::as<arma::vec>(prediction_cluster_phi);
      true_cluster_tau2 = Rcpp::as<arma::vec>(prediction_cluster_tau2);
      int true_clusters = static_cast<int>(true_cluster_phi.n_elem);
      if (true_prediction_cluster.n_elem !=
            static_cast<arma::uword>(n_prediction) ||
          true_cluster_tau2.n_elem != true_cluster_phi.n_elem ||
          true_clusters < 1 || arma::any(true_prediction_cluster < 1) ||
          arma::any(true_prediction_cluster > true_clusters) ||
          arma::any(true_cluster_phi <= 0.0) ||
          arma::any(true_cluster_phi >= 1.0) ||
          arma::any(true_cluster_tau2 <= 0.0)) {
        Rcpp::stop("Invalid prediction cluster truth.");
      }
      phi_truth_matrix.set_size(1, n_prediction);
      tau2_truth_matrix.set_size(1, n_prediction);
      for (int prediction = 0; prediction < n_prediction; ++prediction) {
        int label = true_prediction_cluster(prediction) - 1;
        phi_truth_matrix(0, prediction) = true_cluster_phi(label);
        tau2_truth_matrix(0, prediction) = true_cluster_tau2(label);
      }
    }
    List phi_summary = summarize_prediction_draws_cpp(
      prediction_phi_draws,
      1,
      n_prediction,
      has_cluster_truth ? Rcpp::wrap(phi_truth_matrix) : R_NilValue
    );
    List tau2_summary = summarize_prediction_draws_cpp(
      prediction_tau2_draws,
      1,
      n_prediction,
      has_cluster_truth ? Rcpp::wrap(tau2_truth_matrix) : R_NilValue
    );
    List regime_output = List::create(
      Named("phi") = phi_summary,
      Named("tau2") = tau2_summary
    );

    if (has_cluster_truth) {
      int true_clusters = static_cast<int>(true_cluster_phi.n_elem);
      arma::mat class_probability(
        n_prediction,
        true_clusters,
        arma::fill::zeros
      );
      double phi_scale = std::max(
        true_cluster_phi.max() - true_cluster_phi.min(),
        1e-8
      );
      arma::vec log_true_tau2 = arma::log(true_cluster_tau2);
      double tau2_scale = std::max(
        log_true_tau2.max() - log_true_tau2.min(),
        1e-8
      );
      for (int draw = 0; draw < prediction_draw_count; ++draw) {
        for (int prediction = 0; prediction < n_prediction; ++prediction) {
          double phi_value = prediction_phi_draws(draw, prediction);
          double log_tau2_value = std::log(
            static_cast<double>(prediction_tau2_draws(draw, prediction))
          );
          int closest = 0;
          double closest_distance = arma::datum::inf;
          for (int candidate = 0; candidate < true_clusters; ++candidate) {
            double phi_difference =
              (phi_value - true_cluster_phi(candidate)) / phi_scale;
            double tau2_difference =
              (log_tau2_value - log_true_tau2(candidate)) / tau2_scale;
            double distance = phi_difference * phi_difference +
              tau2_difference * tau2_difference;
            if (distance < closest_distance) {
              closest_distance = distance;
              closest = candidate;
            }
          }
          class_probability(prediction, closest) += 1.0;
        }
      }
      class_probability /= static_cast<double>(prediction_draw_count);
      arma::ivec modal_class(n_prediction);
      arma::vec true_class_probability(n_prediction);
      for (int prediction = 0; prediction < n_prediction; ++prediction) {
        modal_class(prediction) =
          class_probability.row(prediction).index_max() + 1;
        true_class_probability(prediction) = class_probability(
          prediction,
          true_prediction_cluster(prediction) - 1
        );
      }
      regime_output["class_probability"] = class_probability;
      regime_output["modal_class"] = modal_class;
      regime_output["true_class_probability"] = true_class_probability;
      regime_output["classification_accuracy"] = arma::mean(
        arma::conv_to<arma::vec>::from(
          modal_class == true_prediction_cluster
        )
      );
    }
    prediction_output["regime"] = regime_output;
    prediction_output["locations"] = prediction_location_matrix;
    prediction_output["retained_draws"] = prediction_draw_count;
    prediction_output["thin"] = prediction_thin;
    prediction_output["seconds"] = prediction_seconds;
    prediction_output["working_memory_mb"] =
      3.0 * prediction_draw_count * prediction_targets * sizeof(float) /
        (1024.0 * 1024.0) +
      2.0 * prediction_draw_count * n_prediction * sizeof(float) /
        (1024.0 * 1024.0);
    out["prediction"] = prediction_output;
  } else {
    out["prediction"] = R_NilValue;
  }
  return out;
}


