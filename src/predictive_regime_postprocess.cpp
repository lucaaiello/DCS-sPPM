// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

struct SpatialStatistics {
  int count = 0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum11 = 0.0;
  double sum12 = 0.0;
  double sum22 = 0.0;
};

inline void add_location(SpatialStatistics& statistics, double first, double second) {
  statistics.count += 1;
  statistics.sum1 += first;
  statistics.sum2 += second;
  statistics.sum11 += first * first;
  statistics.sum12 += first * second;
  statistics.sum22 += second * second;
}

inline double log_multivariate_gamma_2(double value) {
  return std::log(M_PI) + std::lgamma(value) + std::lgamma(value - 0.5);
}

inline double cohesion3_log(
    const SpatialStatistics& statistics,
    double mu01,
    double mu02,
    double k0,
    double v0,
    double l001,
    double l002,
    double l004) {
  const int count = statistics.count;
  if (count <= 0) return 0.0;
  const double n = static_cast<double>(count);
  const double mean1 = statistics.sum1 / n;
  const double mean2 = statistics.sum2 / n;
  const double scatter11 = std::max(0.0, statistics.sum11 - statistics.sum1 * statistics.sum1 / n);
  const double scatter12 = statistics.sum12 - statistics.sum1 * statistics.sum2 / n;
  const double scatter22 = std::max(0.0, statistics.sum22 - statistics.sum2 * statistics.sum2 / n);
  const double k_n = k0 + n;
  const double v_n = v0 + n;
  const double delta1 = mean1 - mu01;
  const double delta2 = mean2 - mu02;
  const double multiplier = k0 * n / k_n;
  const double ln11 = l001 + scatter11 + multiplier * delta1 * delta1;
  const double ln12 = l002 + scatter12 + multiplier * delta1 * delta2;
  const double ln22 = l004 + scatter22 + multiplier * delta2 * delta2;
  const double det_l0 = l001 * l004 - l002 * l002;
  const double det_ln = ln11 * ln22 - ln12 * ln12;
  if (!(det_l0 > 0.0) || !(det_ln > 0.0)) {
    return -std::numeric_limits<double>::infinity();
  }
  return -n * std::log(M_PI) +
    log_multivariate_gamma_2(0.5 * v_n) -
    log_multivariate_gamma_2(0.5 * v0) +
    0.5 * v0 * std::log(det_l0) -
    0.5 * v_n * std::log(det_ln) +
    std::log(k0) - std::log(k_n);
}

// [[Rcpp::export]]
Rcpp::List predictive_regime_summary_cpp(
    const Rcpp::IntegerMatrix& allocation,
    const Rcpp::List& phi_star,
    const Rcpp::List& tau2_star,
    const Rcpp::NumericVector& mass,
    const Rcpp::NumericMatrix& observed_similarity,
    const Rcpp::NumericMatrix& prediction_similarity,
    const Rcpp::IntegerVector& vi_partition,
    const Rcpp::IntegerVector& draw_indices,
    double k0 = 1.0,
    double v0 = 4.0,
    double a_phi = 1.0,
    double b_phi = 1.0,
    double a_tau = 2.0,
    double b_tau = 0.35,
    int threads = 1,
    std::string model_type = "dcs_sppm") {
  const bool uses_partition = model_type != "k1_common_shock";
  const bool uses_spatial_partition =
    model_type == "dcs_sppm" || model_type == "sppm_no_common_shock";
  if (model_type != "dcs_sppm" && model_type != "ppm_common_shock" &&
      model_type != "k1_common_shock" &&
      model_type != "sppm_no_common_shock") {
    Rcpp::stop("Unknown model_type.");
  }
  const int retained = allocation.nrow();
  const int sites = allocation.ncol();
  const int predictions = prediction_similarity.nrow();
  if (observed_similarity.nrow() != sites || observed_similarity.ncol() != 2 ||
      prediction_similarity.ncol() != 2 || vi_partition.size() != sites ||
      phi_star.size() != retained || tau2_star.size() != retained ||
      mass.size() != retained || draw_indices.size() < 1) {
    Rcpp::stop("Incompatible predictive-regime inputs.");
  }
  const int vi_groups = *std::max_element(vi_partition.begin(), vi_partition.end());
  if (vi_groups < 1 || *std::min_element(vi_partition.begin(), vi_partition.end()) < 1) {
    Rcpp::stop("VI labels must be consecutive positive integers.");
  }

  double mu01 = 0.0;
  double mu02 = 0.0;
  double singleton_log_reference = 0.0;
  if (uses_spatial_partition) {
    for (int site = 0; site < sites; ++site) {
      mu01 += observed_similarity(site, 1) / static_cast<double>(sites);
      mu02 += observed_similarity(site, 0) / static_cast<double>(sites);
    }
    SpatialStatistics singleton;
    add_location(singleton, prediction_similarity(0, 1), prediction_similarity(0, 0));
    singleton_log_reference = cohesion3_log(
      singleton, mu01, mu02, k0, v0, 1.0, 0.0, 1.0
    );
  }

  Rcpp::NumericMatrix regime_probability(predictions, vi_groups + 1);
  Rcpp::NumericVector phi_first(predictions), phi_second(predictions);
  Rcpp::NumericVector tau_first(predictions), tau_second(predictions);
  const double base_phi_mean = a_phi / (a_phi + b_phi);
  const double base_phi_second = a_phi * (a_phi + 1.0) /
    ((a_phi + b_phi) * (a_phi + b_phi + 1.0));
  if (a_tau <= 1.0 || a_tau <= 0.5) {
    Rcpp::stop("The inverse-gamma base measure requires a_tau > 1.");
  }
  const double base_tau_mean = std::sqrt(b_tau) *
    std::tgamma(a_tau - 0.5) / std::tgamma(a_tau);
  const double base_tau_second = b_tau / (a_tau - 1.0);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, threads));
#endif

  for (int selected = 0; selected < draw_indices.size(); ++selected) {
    const int draw = draw_indices[selected] - 1;
    if (draw < 0 || draw >= retained) Rcpp::stop("Invalid draw index.");
    int clusters = 0;
    for (int site = 0; site < sites; ++site) {
      clusters = std::max(clusters, allocation(draw, site));
    }
    if (!uses_partition && clusters != 1) {
      Rcpp::stop("The k1_common_shock model must contain one cluster.");
    }
    Rcpp::NumericVector draw_phi = phi_star[draw];
    Rcpp::NumericVector draw_tau2 = tau2_star[draw];
    if (draw_phi.size() < clusters || draw_tau2.size() < clusters) {
      Rcpp::stop("Representative parameters do not align with allocation labels.");
    }

    std::vector<SpatialStatistics> cluster_statistics(clusters);
    std::vector<int> cluster_sizes(clusters, 0);
    std::vector<std::vector<int>> overlap_counts(
      clusters, std::vector<int>(vi_groups, 0)
    );
    for (int site = 0; site < sites; ++site) {
      const int cluster = allocation(draw, site) - 1;
      const int vi_group = vi_partition[site] - 1;
      if (cluster < 0 || cluster >= clusters) Rcpp::stop("Invalid allocation label.");
      add_location(
        cluster_statistics[cluster],
        observed_similarity(site, 1),
        observed_similarity(site, 0)
      );
      cluster_sizes[cluster] += 1;
      overlap_counts[cluster][vi_group] += 1;
    }
    std::vector<double> old_log_cohesion(clusters);
    std::vector<std::vector<std::pair<int, double>>> overlap_proportions(clusters);
    for (int cluster = 0; cluster < clusters; ++cluster) {
      old_log_cohesion[cluster] = uses_spatial_partition ? cohesion3_log(
        cluster_statistics[cluster], mu01, mu02, k0, v0, 1.0, 0.0, 1.0
      ) : 0.0;
      for (int group = 0; group < vi_groups; ++group) {
        if (overlap_counts[cluster][group] > 0) {
          overlap_proportions[cluster].push_back(std::make_pair(
            group,
            static_cast<double>(overlap_counts[cluster][group]) /
              static_cast<double>(cluster_sizes[cluster])
          ));
        }
      }
    }

    const double draw_mass = mass[draw];
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int prediction = 0; prediction < predictions; ++prediction) {
      const double first = prediction_similarity(prediction, 1);
      const double second = prediction_similarity(prediction, 0);
      double singleton_log = 0.0;
      if (uses_spatial_partition) {
        SpatialStatistics prediction_singleton;
        add_location(prediction_singleton, first, second);
        singleton_log = prediction == 0 ? singleton_log_reference :
          cohesion3_log(
            prediction_singleton, mu01, mu02, k0, v0, 1.0, 0.0, 1.0
          );
      }
      std::vector<double> log_weights(clusters + 1);
      double maximum = uses_partition ?
        std::log(draw_mass) + singleton_log :
        -std::numeric_limits<double>::infinity();
      log_weights[clusters] = maximum;
      for (int cluster = 0; cluster < clusters; ++cluster) {
        double spatial_increment = 0.0;
        if (uses_spatial_partition) {
          SpatialStatistics proposed = cluster_statistics[cluster];
          add_location(proposed, first, second);
          spatial_increment = cohesion3_log(
            proposed, mu01, mu02, k0, v0, 1.0, 0.0, 1.0
          ) - old_log_cohesion[cluster];
        }
        log_weights[cluster] =
          std::log(static_cast<double>(cluster_sizes[cluster])) +
          spatial_increment;
        maximum = std::max(maximum, log_weights[cluster]);
      }
      if (!uses_partition) {
        log_weights[0] = 0.0;
        maximum = 0.0;
      }
      double denominator = 0.0;
      for (double value : log_weights) denominator += std::exp(value - maximum);
      const double new_probability = uses_partition ?
        std::exp(log_weights[clusters] - maximum) / denominator : 0.0;
      regime_probability(prediction, vi_groups) += new_probability;
      phi_first[prediction] += new_probability * base_phi_mean;
      phi_second[prediction] += new_probability * base_phi_second;
      tau_first[prediction] += new_probability * base_tau_mean;
      tau_second[prediction] += new_probability * base_tau_second;
      for (int cluster = 0; cluster < clusters; ++cluster) {
        const double probability = std::exp(log_weights[cluster] - maximum) / denominator;
        const double phi_value = draw_phi[cluster];
        const double tau_value = std::sqrt(draw_tau2[cluster]);
        phi_first[prediction] += probability * phi_value;
        phi_second[prediction] += probability * phi_value * phi_value;
        tau_first[prediction] += probability * tau_value;
        tau_second[prediction] += probability * draw_tau2[cluster];
        for (const auto& overlap : overlap_proportions[cluster]) {
          regime_probability(prediction, overlap.first) += probability * overlap.second;
        }
      }
    }
  }

  const double number_of_draws = static_cast<double>(draw_indices.size());
  Rcpp::IntegerVector dominant_regime(predictions);
  Rcpp::NumericVector maximum_probability(predictions), normalized_entropy(predictions);
  Rcpp::NumericVector phi_mean(predictions), phi_sd(predictions);
  Rcpp::NumericVector tau_mean(predictions), tau_sd(predictions), tau2_mean(predictions);
  for (int prediction = 0; prediction < predictions; ++prediction) {
    int dominant = 0;
    double entropy = 0.0;
    for (int group = 0; group < vi_groups + 1; ++group) {
      regime_probability(prediction, group) /= number_of_draws;
      const double probability = regime_probability(prediction, group);
      if (probability > regime_probability(prediction, dominant)) dominant = group;
      if (probability > 0.0) entropy -= probability * std::log(probability);
    }
    dominant_regime[prediction] = dominant + 1;
    maximum_probability[prediction] = regime_probability(prediction, dominant);
    normalized_entropy[prediction] = entropy / std::log(static_cast<double>(vi_groups + 1));
    phi_mean[prediction] = phi_first[prediction] / number_of_draws;
    const double phi_variance = phi_second[prediction] / number_of_draws -
      phi_mean[prediction] * phi_mean[prediction];
    phi_sd[prediction] = std::sqrt(std::max(0.0, phi_variance));
    tau_mean[prediction] = tau_first[prediction] / number_of_draws;
    tau2_mean[prediction] = tau_second[prediction] / number_of_draws;
    const double tau_variance = tau2_mean[prediction] -
      tau_mean[prediction] * tau_mean[prediction];
    tau_sd[prediction] = std::sqrt(std::max(0.0, tau_variance));
  }

  return Rcpp::List::create(
    Rcpp::Named("regime_probability") = regime_probability,
    Rcpp::Named("dominant_regime") = dominant_regime,
    Rcpp::Named("maximum_probability") = maximum_probability,
    Rcpp::Named("normalized_entropy") = normalized_entropy,
    Rcpp::Named("new_regime_probability") = regime_probability(Rcpp::_, vi_groups),
    Rcpp::Named("phi_mean") = phi_mean,
    Rcpp::Named("phi_sd") = phi_sd,
    Rcpp::Named("tau_mean") = tau_mean,
    Rcpp::Named("tau_sd") = tau_sd,
    Rcpp::Named("tau2_mean") = tau2_mean,
    Rcpp::Named("draw_indices") = draw_indices,
    Rcpp::Named("vi_groups") = vi_groups,
    Rcpp::Named("model_type") = model_type,
    Rcpp::Named("weight_rule") = uses_spatial_partition ?
      "spatial PPM cohesion ratio" :
      (uses_partition ? "ordinary PPM cluster-size/mass" : "fixed K=1"),
    Rcpp::Named("projection") =
      (uses_partition ?
        "Draw-cluster probability split across VI groups by observed membership proportions; singleton retained separately" :
        "Single global temporal regime; no label projection or singleton allocation")
  );
}
