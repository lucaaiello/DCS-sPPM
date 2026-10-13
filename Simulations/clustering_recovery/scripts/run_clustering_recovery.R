script_arguments <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
validation_root <- if (length(script_arguments)) {
  dirname(normalizePath(sub("^--file=", "", script_arguments[[1L]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
repeat {
  if (file.exists(file.path(validation_root, "R", "project_paths.R"))) break
  parent <- dirname(validation_root)
  if (identical(parent, validation_root)) stop("Could not locate the manuscript project.")
  validation_root <- parent
}
source(file.path(validation_root, "R", "project_paths.R"))
source(file.path(validation_root, "R", "simulate_dcs_sppm.R"))
source(file.path(validation_root, "R", "sampler_runner.R"))
source(file.path(validation_root, "R", "validation_utils.R"))
source(file.path(validation_root, "R", "pilot_simulation_utils.R"))
source(file.path(validation_root, "R", "full_simulation_utils.R"))
source(file.path(validation_root, "R", "clustering_recovery_utils.R"))

environment_flag <- function(name, default) {
  value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0")) {
    stop("`", name, "` must be true/false or 1/0.")
  }
  value %in% c("true", "1")
}

smoke_test <- environment_flag("DCS_CLUSTER_SMOKE", FALSE)
replications <- as.integer(Sys.getenv(
  "DCS_CLUSTER_REPLICATIONS",
  if (smoke_test) "1" else "100"
))
runs <- as.integer(Sys.getenv(
  "DCS_CLUSTER_RUNS",
  if (smoke_test) "200" else "8000"
))
burn <- as.integer(Sys.getenv(
  "DCS_CLUSTER_BURN",
  if (smoke_test) "100" else "4000"
))
thin <- as.integer(Sys.getenv("DCS_CLUSTER_THIN", "1"))
workers <- as.integer(Sys.getenv("DCS_CLUSTER_WORKERS", "6"))
omp_threads <- as.integer(Sys.getenv("DCS_CLUSTER_OMP_THREADS", "1"))
n_side <- as.integer(Sys.getenv("DCS_CLUSTER_N_SIDE", "12"))
n_times <- as.integer(Sys.getenv("DCS_CLUSTER_N_TIMES", "365"))
n_neighbors <- as.integer(Sys.getenv("DCS_CLUSTER_NEIGHBORS", "10"))
max_hours <- as.numeric(Sys.getenv("DCS_CLUSTER_MAX_HOURS", "18"))
reuse <- environment_flag("DCS_CLUSTER_REUSE", TRUE)
results_name <- Sys.getenv(
  "DCS_CLUSTER_RESULTS_NAME",
  if (smoke_test) "clustering_recovery_smoke" else "clustering_recovery_100rep"
)
sampler_source <- Sys.getenv(
  "DCS_CLUSTER_CPP_FILE",
  "sampler_global_beta_single_w.cpp"
)

integer_settings <- c(
  replications, runs, burn, thin, workers, omp_threads,
  n_side, n_times, n_neighbors
)
if (anyNA(integer_settings) || replications < 1L || runs <= burn || burn < 0L ||
    thin < 1L || workers < 1L || omp_threads < 1L || n_side < 3L ||
    n_times < 2L || n_neighbors < 1L || !is.finite(max_hours) ||
    max_hours <= 0) {
  stop("Invalid clustering-recovery run settings.")
}
if ((n_side^2) %% 3L != 0L) {
  stop("`DCS_CLUSTER_N_SIDE^2` must be divisible by three.")
}
if (!grepl("^[A-Za-z0-9_.-]+$", results_name) ||
    results_name %in% c(".", "..")) {
  stop("`DCS_CLUSTER_RESULTS_NAME` must be a safe directory name.")
}
if (!grepl("^[A-Za-z0-9_.-]+\\.cpp$", sampler_source)) {
  stop("`DCS_CLUSTER_CPP_FILE` must be a C++ filename in `src`.")
}

cpp_file <- file.path(validation_root, "src", sampler_source)
if (!file.exists(cpp_file)) {
  stop("Sampler source file not found: ", cpp_file)
}
results_dir <- file.path(simulation_clustering_results_dir, results_name)
status_dir <- file.path(results_dir, "worker_status")
cache_dir <- Sys.getenv(
  "DCS_CLUSTER_CACHE_DIR",
  file.path(dirname(validation_root), "tmp", "sourcecpp_clustering_recovery")
)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- normalizePath(cache_dir, winslash = "/", mustWork = TRUE)

scenarios <- clustering_recovery_scenarios()
tasks <- clustering_recovery_tasks(replications)
workers <- min(workers, length(tasks))
task_manifest <- do.call(rbind, lapply(tasks, function(task) {
  data.frame(
    task_index = task$task_index,
    scenario = task$scenario,
    geography = task$geography,
    separation = task$separation,
    replicate = task$replicate,
    display_name = task$display_name,
    stringsAsFactors = FALSE
  )
}))
task_manifest$runs <- runs
task_manifest$burn <- burn
task_manifest$thin <- thin
task_manifest$n <- n_side^2
task_manifest$n_times <- n_times
task_manifest$n_neighbors <- n_neighbors
task_manifest$sampler_source <- sampler_source
utils::write.csv(scenarios, file.path(results_dir, "scenario_manifest.csv"), row.names = FALSE)
utils::write.csv(task_manifest, file.path(results_dir, "task_manifest.csv"), row.names = FALSE)

run_settings <- list(
  version = 1L,
  smoke_test = smoke_test,
  replications = replications,
  runs = runs,
  burn = burn,
  thin = thin,
  workers = workers,
  omp_threads = omp_threads,
  n_side = n_side,
  n = n_side^2,
  n_times = n_times,
  n_neighbors = n_neighbors,
  max_hours = max_hours,
  reuse = reuse,
  sampler_source = sampler_source,
  started_at = Sys.time()
)
atomic_save_rds(run_settings, file.path(results_dir, "run_settings.rds"))

if (!requireNamespace("Rcpp", quietly = TRUE) ||
    !requireNamespace("Matrix", quietly = TRUE)) {
  stop("The `Rcpp` and `Matrix` packages are required.")
}
Sys.setenv(OMP_NUM_THREADS = as.character(omp_threads))
Rcpp::sourceCpp(cpp_file, cacheDir = cache_dir, rebuild = FALSE)

fixed_prior <- pilot_fixed_prior()
deadline <- as.numeric(Sys.time()) + max_hours * 60 * 60
worker_groups <- split(tasks, rep(seq_len(workers), length.out = length(tasks)))
status_files <- file.path(
  status_dir,
  sprintf("worker_%02d.rds", seq_len(workers))
)
chain_progress_files <- file.path(
  status_dir,
  sprintf("worker_%02d_chain_progress.txt", seq_len(workers))
)
for (worker in seq_len(workers)) {
  write_clustering_worker_status(
    status_files[worker],
    worker,
    length(worker_groups[[worker]]),
    0L,
    "starting"
  )
  writeLines("0", chain_progress_files[worker])
}

read_worker_status <- function(file, progress_file, worker, total) {
  status <- tryCatch(readRDS(file), error = function(error) NULL)
  if (is.null(status)) {
    return(list(
      worker = worker,
      total = total,
      completed = 0L,
      state = "starting",
      current = "waiting",
      failures = 0L,
      chain_progress = 0
    ))
  }
  chain_progress <- tryCatch(
    suppressWarnings(as.numeric(readLines(progress_file, n = 1L, warn = FALSE))),
    error = function(error) 0
  )
  if (length(chain_progress) != 1L || !is.finite(chain_progress)) {
    chain_progress <- 0
  }
  status$chain_progress <- min(100, max(0, chain_progress))
  status
}

worker_status_line <- function(status, width = 24L) {
  chain_fraction <- if (identical(status$state, "running")) {
    status$chain_progress / 100
  } else {
    0
  }
  fraction <- if (status$total > 0L) {
    (status$completed + chain_fraction) / status$total
  } else {
    1
  }
  filled <- min(width, floor(width * fraction))
  bar <- paste0(
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", width - filled), collapse = "")
  )
  current <- substr(status$current, 1L, 58L)
  sprintf(
    "Worker %d [%s] %3d/%3d | chain %5.1f%% | %-7s | failures %d | %s",
    status$worker,
    bar,
    status$completed,
    status$total,
    status$chain_progress,
    status$state,
    status$failures,
    current
  )
}

draw_worker_dashboard <- local({
  initialized <- FALSE
  function(statuses) {
    lines <- vapply(statuses, worker_status_line, character(1))
    if (isatty(stdout())) {
      if (initialized) {
        cat(sprintf("\033[%dA", length(lines)))
      }
      for (line in lines) {
        cat("\r\033[2K", line, "\n", sep = "")
      }
      initialized <<- TRUE
      flush.console()
    } else if (!initialized || all(vapply(statuses, `[[`, character(1), "state") == "done")) {
      cat(paste(lines, collapse = "\n"), "\n")
      initialized <<- TRUE
    }
  }
})

worker_job <- function(
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
  run_clustering_recovery_worker(
    tasks = tasks,
    worker = worker,
    status_file = status_file,
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
}

if (workers == 1L) {
  statuses <- list(worker_job(
    tasks = worker_groups[[1]],
    worker = 1L,
    status_file = status_files[1],
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
    chain_progress_file = chain_progress_files[1]
  ))
} else {
  worker_log <- file.path(results_dir, "worker_console.log")
  cluster <- parallel::makePSOCKcluster(workers, outfile = worker_log)
  statuses <- tryCatch({
    parallel::clusterCall(
      cluster,
      function(validation_root, cache_dir, cpp_file, omp_threads) {
        Sys.setenv(OMP_NUM_THREADS = as.character(omp_threads))
        source(file.path(validation_root, "R", "simulate_dcs_sppm.R"))
        source(file.path(validation_root, "R", "sampler_runner.R"))
        source(file.path(validation_root, "R", "validation_utils.R"))
        source(file.path(validation_root, "R", "pilot_simulation_utils.R"))
        source(file.path(validation_root, "R", "full_simulation_utils.R"))
        source(file.path(validation_root, "R", "clustering_recovery_utils.R"))
        Rcpp::sourceCpp(cpp_file, cacheDir = cache_dir, rebuild = FALSE)
        TRUE
      },
      validation_root = validation_root,
      cache_dir = cache_dir,
      cpp_file = cpp_file,
      omp_threads = omp_threads
    )

    for (worker in seq_len(workers)) {
      parallel:::sendCall(
        cluster[[worker]],
        worker_job,
        list(
          tasks = worker_groups[[worker]],
          worker = worker,
          status_file = status_files[worker],
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
          chain_progress_file = chain_progress_files[worker]
        ),
        tag = worker
      )
    }

    repeat {
      worker_statuses <- lapply(seq_len(workers), function(worker) {
        read_worker_status(
          status_files[worker],
          chain_progress_files[worker],
          worker,
          length(worker_groups[[worker]])
        )
      })
      draw_worker_dashboard(worker_statuses)
      if (all(vapply(worker_statuses, `[[`, character(1), "state") == "done")) {
        break
      }
      if (as.numeric(Sys.time()) > deadline + 15 * 60) {
        stop("Workers did not finish within 15 minutes of the requested deadline.")
      }
      Sys.sleep(2)
    }

    lapply(cluster, parallel:::recvResult)
  }, finally = {
    parallel::stopCluster(cluster)
  })
}

status_table <- do.call(rbind, statuses)
utils::write.csv(status_table, file.path(results_dir, "run_status.csv"), row.names = FALSE)
print(table(status_table$status))
