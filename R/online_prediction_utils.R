prediction_quantile <- function(summary, probability) {
  stopifnot(
    is.list(summary),
    is.matrix(summary$mean),
    is.matrix(summary$sd),
    length(probability) == 1L,
    is.finite(probability),
    probability > 0,
    probability < 1
  )
  summary$mean + stats::qnorm(probability) * summary$sd
}

prediction_interval <- function(summary, level = 0.95) {
  stopifnot(length(level) == 1L, is.finite(level), level > 0, level < 1)
  tail_probability <- (1 - level) / 2
  list(
    lower = prediction_quantile(summary, tail_probability),
    upper = prediction_quantile(summary, 1 - tail_probability),
    level = level,
    method = summary$quantile_method %||% "moment-matched Gaussian"
  )
}

prediction_monitor_diagnostics <- function(prediction) {
  if (is.null(prediction$monitor)) {
    stop("The prediction output does not contain representative traces.")
  }
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Install the `coda` package to compute prediction diagnostics.")
  }
  monitor <- prediction$monitor
  summarize_process <- function(draws, process) {
    draws <- as.matrix(draws)
    effective_size <- as.numeric(coda::effectiveSize(draws))
    posterior_sd <- apply(draws, 2L, stats::sd)
    data.frame(
      process = process,
      monitor = seq_len(ncol(draws)),
      time = monitor$time,
      location = monitor$location,
      mean = colMeans(draws),
      sd = posterior_sd,
      ess = effective_size,
      mcse = posterior_sd / sqrt(pmax(effective_size, 1)),
      mcse_over_sd = 1 / sqrt(pmax(effective_size, 1))
    )
  }
  rbind(
    summarize_process(monitor$u_draws, "u"),
    summarize_process(monitor$w_draws, "w"),
    summarize_process(monitor$y_draws, "y")
  )
}

`%||%` <- function(left, right) {
  if (is.null(left)) right else left
}

reconstruct_site_parameter <- function(fit, parameter = c("phi", "tau2")) {
  parameter <- match.arg(parameter)
  representatives <- if (parameter == "phi") fit$phistar else fit$tau2star
  labels <- as.matrix(fit$cluster)
  if (is.null(representatives) || length(representatives) != nrow(labels)) {
    stop("The representative chain is not aligned with the allocation chain.")
  }
  reconstructed <- matrix(NA_real_, nrow(labels), ncol(labels))
  for (draw in seq_len(nrow(labels))) {
    draw_labels <- labels[draw, ]
    draw_representatives <- as.numeric(representatives[[draw]])
    if (anyNA(draw_labels) || any(draw_labels < 1L) ||
        any(draw_labels > length(draw_representatives))) {
      stop("Invalid allocation/representative alignment at draw ", draw, ".")
    }
    reconstructed[draw, ] <- draw_representatives[draw_labels]
  }
  reconstructed
}
