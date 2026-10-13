clustering_recovery_geographies <- function() {
  data.frame(
    geography = c("regular_strips", "irregular_regions", "random_allocation"),
    description = c(
      "Three balanced vertical strips",
      "Three balanced contiguous regions with nonlinear boundaries",
      "Three balanced labels assigned independently of location"
    ),
    stringsAsFactors = FALSE
  )
}

clustering_recovery_separations <- function() {
  list(
    very_similar = list(
      phi = c(0.45, 0.55, 0.65),
      tau2 = c(0.30, 0.35, 0.40)
    ),
    intermediate = list(
      phi = c(0.25, 0.55, 0.80),
      tau2 = c(0.20, 0.35, 0.55)
    ),
    very_different = list(
      phi = c(0.10, 0.55, 0.90),
      tau2 = c(0.10, 0.35, 0.80)
    )
  )
}

clustering_recovery_locations <- function(n_side = 12L) {
  n_side <- as.integer(n_side)
  if (n_side < 3L || (n_side^2) %% 3L != 0L) {
    stop("`n_side^2` must be divisible by three and `n_side` must be at least three.")
  }
  locations <- expand.grid(
    x = seq(0, 5, length.out = n_side),
    y = seq(0, 5, length.out = n_side),
    KEEP.OUT.ATTRS = FALSE
  )
  locations <- locations[order(locations$x, locations$y), , drop = FALSE]
  rownames(locations) <- NULL
  as.matrix(locations)
}

balanced_labels_from_score <- function(score) {
  n <- length(score)
  if (n %% 3L != 0L || any(!is.finite(score))) {
    stop("Balanced three-region labels require a finite score and n divisible by three.")
  }
  ordered <- order(score, seq_along(score))
  labels <- integer(n)
  labels[ordered] <- rep(1:3, each = n / 3L)
  labels
}

clustering_recovery_truth <- function(
    geography,
    replicate,
    locations = clustering_recovery_locations(),
    seed_base = 20310000L
) {
  geography <- match.arg(
    geography,
    clustering_recovery_geographies()$geography
  )
  locations <- as.matrix(locations)
  n <- nrow(locations)
  if (n %% 3L != 0L) {
    stop("The number of locations must be divisible by three.")
  }

  x <- locations[, 1]
  y <- locations[, 2]
  region <- switch(
    geography,
    regular_strips = balanced_labels_from_score(x),
    irregular_regions = balanced_labels_from_score(
      x - 0.85 * sin(2 * pi * y / 5) - 0.20 * sin(4 * pi * y / 5)
    ),
    random_allocation = {
      set.seed(seed_base + 10000L + as.integer(replicate))
      sample(rep(1:3, each = n / 3L), size = n, replace = FALSE)
    }
  )

  geography_index <- match(
    geography,
    clustering_recovery_geographies()$geography
  )
  set.seed(seed_base + 1000L * geography_index + as.integer(replicate))
  region_to_parameter <- sample(1:3, size = 3L, replace = FALSE)
  true_cluster <- region_to_parameter[region]

  list(
    region = as.integer(region),
    true_cluster = as.integer(true_cluster),
    region_to_parameter = as.integer(region_to_parameter)
  )
}

clustering_recovery_scenarios <- function() {
  geographies <- clustering_recovery_geographies()
  separations <- clustering_recovery_separations()
  scenario <- expand.grid(
    geography = geographies$geography,
    separation = names(separations),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scenario$scenario <- paste(scenario$geography, scenario$separation, sep = "__")
  scenario$cluster_phi <- vapply(
    scenario$separation,
    function(level) paste(separations[[level]]$phi, collapse = ";"),
    character(1)
  )
  scenario$cluster_tau2 <- vapply(
    scenario$separation,
    function(level) paste(separations[[level]]$tau2, collapse = ";"),
    character(1)
  )
  scenario[, c(
    "scenario", "geography", "separation", "cluster_phi", "cluster_tau2"
  )]
}

clustering_recovery_tasks <- function(replications = 100L) {
  scenarios <- clustering_recovery_scenarios()
  tasks <- vector("list", nrow(scenarios) * replications)
  task_index <- 0L
  for (replicate in seq_len(replications)) {
    for (scenario_index in seq_len(nrow(scenarios))) {
      task_index <- task_index + 1L
      row <- scenarios[scenario_index, , drop = FALSE]
      tasks[[task_index]] <- list(
        task_index = task_index,
        scenario_index = scenario_index,
        scenario = row$scenario,
        geography = row$geography,
        separation = row$separation,
        replicate = replicate,
        display_name = sprintf("%s__r%03d", row$scenario, replicate)
      )
    }
  }
  tasks
}

simulate_clustering_recovery_data <- function(
    task,
    n_side = 12L,
    n_times = 365L,
    n_neighbors = 10L,
    seed_base = 20310000L
) {
  locations <- clustering_recovery_locations(n_side)
  truth <- clustering_recovery_truth(
    task$geography,
    task$replicate,
    locations,
    seed_base
  )
  representatives <- clustering_recovery_separations()[[task$separation]]
  geography_index <- match(
    task$geography,
    clustering_recovery_geographies()$geography
  )
  data_seed <- seed_base + 100000L * geography_index + task$replicate

  simulation <- simulate_dcs_sppm_global_beta(
    n_side = n_side,
    n_times = n_times,
    n_neighbors = n_neighbors,
    cluster_phi = representatives$phi,
    cluster_tau2 = representatives$tau2,
    spatial_dgp = "nngp",
    seed = data_seed,
    locations = locations,
    true_cluster = truth$true_cluster
  )
  simulation$true$region <- truth$region
  simulation$true$region_to_parameter <- truth$region_to_parameter
  simulation$true$geography <- task$geography
  simulation$true$separation <- task$separation
  simulation$data_seed <- data_seed
  simulation
}

clustering_recovery_partition_metrics <- function(fit, simulation, task) {
  similarity <- posterior_similarity_matrix(fit$cluster)
  representative <- dahl_partition(fit$cluster, similarity)$partition
  truth <- simulation$true$cluster
  pairwise <- pairwise_partition_accuracy(truth, representative)
  selected <- upper.tri(similarity)
  true_similarity <- outer(truth, truth, FUN = "==")
  ari <- adjusted_rand_index(truth, representative)

  metrics <- data.frame(
    scenario = task$scenario,
    geography = task$geography,
    separation = task$separation,
    replicate = task$replicate,
    adjusted_rand_index = ari,
    variation_of_information = variation_of_information(truth, representative),
    pairwise_sensitivity = unname(pairwise["sensitivity"]),
    pairwise_specificity = unname(pairwise["specificity"]),
    coclustering_brier = mean(
      (similarity[selected] - as.numeric(true_similarity[selected]))^2
    ),
    exact_partition_recovery = isTRUE(all.equal(ari, 1)),
    true_clusters = length(unique(truth)),
    posterior_mean_clusters = mean(fit$k_clust),
    posterior_mode_clusters = as.integer(names(which.max(table(fit$k_clust)))),
    spatial_acceptance = fit$acc_spatial,
    phi_acceptance = fit$accphi,
    elapsed_seconds = fit$elapsed_seconds,
    stringsAsFactors = FALSE
  )

  list(
    metrics = metrics,
    representative_partition = representative,
    posterior_similarity = similarity
  )
}

summarize_clustering_recovery_partitions <- function(metrics) {
  groups <- split(metrics, metrics$scenario)
  rows <- lapply(groups, function(group) {
    recovery <- mean(group$exact_partition_recovery)
    data.frame(
      scenario = group$scenario[1],
      geography = group$geography[1],
      separation = group$separation[1],
      completed_replications = nrow(group),
      mean_ari = mean(group$adjusted_rand_index),
      ari_mcse = stats::sd(group$adjusted_rand_index) / sqrt(nrow(group)),
      median_ari = stats::median(group$adjusted_rand_index),
      mean_variation_of_information = mean(group$variation_of_information),
      mean_pairwise_sensitivity = mean(group$pairwise_sensitivity),
      mean_pairwise_specificity = mean(group$pairwise_specificity),
      mean_coclustering_brier = mean(group$coclustering_brier),
      exact_recovery_rate = recovery,
      exact_recovery_mcse = sqrt(recovery * (1 - recovery) / nrow(group)),
      mean_posterior_clusters = mean(group$posterior_mean_clusters),
      correct_modal_cluster_rate = mean(
        group$posterior_mode_clusters == group$true_clusters
      ),
      mean_elapsed_seconds = mean(group$elapsed_seconds),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary[order(summary$geography, summary$separation), ]
}

write_clustering_worker_status <- function(
    file,
    worker,
    total,
    completed,
    state,
    current = "waiting",
    failures = 0L
) {
  atomic_save_rds(
    list(
      worker = as.integer(worker),
      total = as.integer(total),
      completed = as.integer(completed),
      state = state,
      current = current,
      failures = as.integer(failures),
      updated_at = Sys.time()
    ),
    file
  )
}

run_clustering_recovery_task <- function(
    task,
    results_dir,
    runs,
    burn,
    thin,
    n_side,
    n_times,
    n_neighbors,
    fixed_prior,
    deadline,
    reuse,
    sampler_source,
    chain_progress_file = ""
) {
  task_dir <- file.path(
    results_dir,
    task$scenario,
    sprintf("replicate_%03d", task$replicate)
  )
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  metrics_file <- file.path(task_dir, "metrics.rds")
  failure_file <- file.path(task_dir, "failure.rds")
  chain_seed <- 20360000L + 10000L * task$scenario_index + task$replicate
  expected_settings <- list(
    version = 1L,
    sampler_source = sampler_source,
    scenario = task$scenario,
    geography = task$geography,
    separation = task$separation,
    replicate = task$replicate,
    n_side = n_side,
    n = n_side^2,
    n_times = n_times,
    n_neighbors = n_neighbors,
    runs = runs,
    burn = burn,
    thin = thin,
    chain_seed = chain_seed
  )

  if (reuse && file.exists(metrics_file)) {
    cached <- tryCatch(readRDS(metrics_file), error = function(error) NULL)
    if (!is.null(cached) && identical(cached$settings, expected_settings)) {
      return(data.frame(
        scenario = task$scenario,
        replicate = task$replicate,
        status = "reused",
        stringsAsFactors = FALSE
      ))
    }
  }
  if (file.exists(metrics_file)) {
    unlink(metrics_file)
  }
  if (as.numeric(Sys.time()) >= deadline) {
    return(data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      status = "deadline",
      stringsAsFactors = FALSE
    ))
  }

  tryCatch({
    simulation <- simulate_clustering_recovery_data(
      task = task,
      n_side = n_side,
      n_times = n_times,
      n_neighbors = n_neighbors
    )
    fit <- run_global_beta_validation_chain(
      simulation = simulation,
      runs = runs,
      burn = burn,
      thin = thin,
      beta_update_mode = "collapsed",
      sigma2_update_mode = "collapsed",
      tau2_update_mode = "pairwise_interwoven",
      spatial_update_mode = "collapsed",
      spatial_process_model = "nngp",
      prior_override = fixed_prior,
      seed = chain_seed,
      overdispersed_initialization = TRUE,
      save_latent = FALSE,
      save_w = FALSE,
      progress_file = chain_progress_file,
      show_progress = FALSE
    )
    partition <- clustering_recovery_partition_metrics(fit, simulation, task)
    parameters <- full_parameter_metrics(
      fit,
      simulation,
      task$scenario,
      "clustering_recovery",
      task$replicate
    )
    compact <- list(
      settings = expected_settings,
      truth = list(
        cluster = simulation$true$cluster,
        region = simulation$true$region,
        region_to_parameter = simulation$true$region_to_parameter,
        cluster_phi = simulation$true$cluster_phi,
        cluster_tau2 = simulation$true$cluster_tau2
      ),
      partition = partition$metrics,
      parameters = parameters,
      representative_partition = partition$representative_partition,
      posterior_similarity = partition$posterior_similarity,
      mcmc = list(
        retained_draws = nrow(fit$cluster),
        k_clust_mean = mean(fit$k_clust),
        k_clust_ess = effective_sample_size(fit$k_clust),
        openmp_threads = fit$openmp_threads
      )
    )
    atomic_save_rds(compact, metrics_file)
    if (file.exists(failure_file)) {
      unlink(failure_file)
    }
    rm(fit, simulation)
    invisible(gc(FALSE))
    data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      status = "completed",
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    atomic_save_rds(
      list(settings = expected_settings, message = conditionMessage(error)),
      failure_file
    )
    data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      status = "failed",
      stringsAsFactors = FALSE
    )
  })
}

run_clustering_recovery_worker <- function(
    tasks,
    worker,
    status_file,
    results_dir,
    runs,
    burn,
    thin,
    n_side,
    n_times,
    n_neighbors,
    fixed_prior,
    deadline,
    reuse,
    sampler_source,
    chain_progress_file
) {
  total <- length(tasks)
  completed <- 0L
  failures <- 0L
  statuses <- vector("list", total)

  for (index in seq_along(tasks)) {
    task <- tasks[[index]]
    write_clustering_worker_status(
      status_file,
      worker,
      total,
      completed,
      "running",
      task$display_name,
      failures
    )
    if (nzchar(chain_progress_file)) {
      writeLines("0", chain_progress_file)
    }
    statuses[[index]] <- run_clustering_recovery_task(
      task = task,
      results_dir = results_dir,
      runs = runs,
      burn = burn,
      thin = thin,
      n_side = n_side,
      n_times = n_times,
      n_neighbors = n_neighbors,
      fixed_prior = fixed_prior,
      deadline = deadline,
      reuse = reuse,
      sampler_source = sampler_source,
      chain_progress_file = chain_progress_file
    )
    completed <- completed + 1L
    failures <- failures + as.integer(statuses[[index]]$status == "failed")
    if (nzchar(chain_progress_file)) {
      writeLines("0", chain_progress_file)
    }
    write_clustering_worker_status(
      status_file,
      worker,
      total,
      completed,
      if (completed == total) "done" else "running",
      task$display_name,
      failures
    )
  }
  do.call(rbind, statuses)
}
