prepare_europe_model_comparison_data <- function(
    validation_root,
    validation_results_dir,
    holdout_result_name = "spatial_holdout_100_coordinate_groups",
    preparation_name = "production_10000_single_chain",
    n_neighbors = 10L
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required.")
  }
  preparation_file <- file.path(validation_root, "Data", paste0(preparation_name, "_preparation.rds"))
  if (!file.exists(preparation_file)) stop("Missing frozen preparation input: ", preparation_file)
  preparation <- readRDS(preparation_file)
  if (!identical(preparation$schema_version, 1L)) stop("Unsupported preparation input schema.")
  station_information <- data.table::as.data.table(
    preparation$station_information
  )
  analysis_dates <- as.Date(preparation$analysis_dates)
  scaling <- preparation$data_scaling
  unique_data_file <- Sys.getenv("DCS_REAL_FULL_DATA_FILE", "")
  unique_sample <- nzchar(unique_data_file)
  if (unique_sample) {
    unique_data <- readRDS(unique_data_file)
    selected <- data.table::as.data.table(unique_data$station_information)
    if (anyDuplicated(selected$Samplingpoint) ||
        anyDuplicated(selected[, .(Longitude, Latitude)])) {
      stop("The selected sample must contain unique identifiers and coordinates.")
    }
    original_inventory <- station_information
    weather_match <- vapply(seq_len(nrow(selected)), function(index) {
      same_station <- which(original_inventory$Samplingpoint == selected$Samplingpoint[index])
      if (length(same_station) == 1L) return(same_station)
      which.min((original_inventory$Longitude - selected$Longitude[index])^2 +
        (original_inventory$Latitude - selected$Latitude[index])^2)
    }, integer(1L))
    selected[, c("weather_longitude", "weather_latitude", "weather_cell") :=
      original_inventory[weather_match, .(weather_longitude, weather_latitude, weather_cell)]]
    station_information <- selected
    analysis_dates <- as.Date(unique_data$analysis_dates)
    if (!identical(as.character(colnames(unique_data$Y_imputed)), selected$Samplingpoint) ||
        !identical(as.character(rownames(unique_data$Y_imputed)), as.character(analysis_dates))) {
      stop("Selected response matrix does not align with stations and dates.")
    }
  }
  if (is.null(holdout_result_name)) {
    original_holdout <- NULL
    holdout <- integer()
  } else {
    holdout_file <- file.path(
      validation_results_dir, "Real_data_analysis", "holdout", "results", holdout_result_name,
      "spatial_holdout_result.rds"
    )
    if (!file.exists(holdout_file)) stop("Missing holdout result: ", holdout_file)
    original_holdout <- readRDS(holdout_file)
    holdout_ids <- original_holdout$held_out_station_information$Samplingpoint
    holdout <- match(holdout_ids, station_information$Samplingpoint)
    if (anyNA(holdout) || anyDuplicated(holdout)) {
      stop("Saved holdout stations do not match the production station inventory.")
    }
    holdout <- sort(holdout)
  }
  training <- setdiff(seq_len(nrow(station_information)), holdout)

  if (unique_sample) {
    Y_raw <- unique_data$Y_imputed
  } else {
  pm10 <- data.table::as.data.table(readRDS(
    file.path(validation_root, "Data", "pm10_combined_2024.rds")
  ))
  pm10[, Date := as.Date(as.numeric(unclass(Date)), origin = "1970-01-01")]
  pm10 <- pm10[
    Samplingpoint %in% station_information$Samplingpoint &
      Date %in% analysis_dates
  ]
  Y_raw <- matrix(
    NA_real_, length(analysis_dates), nrow(station_information),
    dimnames = list(
      as.character(analysis_dates), station_information$Samplingpoint
    )
  )
  row_index <- match(pm10$Date, analysis_dates)
  column_index <- match(pm10$Samplingpoint, station_information$Samplingpoint)
  if (anyDuplicated(cbind(row_index, column_index))) {
    stop("Duplicate station-date PM10 values.")
  }
  Y_raw[cbind(row_index, column_index)] <- pm10$Value
  impute_series <- function(values) {
    observed <- which(is.finite(values))
    if (length(observed) < 2L) stop("A station has fewer than two observations.")
    stats::approx(observed, values[observed], seq_along(values), rule = 2)$y
  }
  Y_raw <- apply(Y_raw, 2L, impute_series)
  }
  Y <- (Y_raw - scaling$response_center) / scaling$response_scale

  weather <- data.table::as.data.table(readRDS(
    file.path(validation_root, "Data", "era5_combined_2024.rds")
  ))
  weather[, Date := as.Date("2024-01-01") + Date]
  weather[, `:=`(
    Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)
  )]
  weather_cells <- unique(station_information[, .(
    weather_longitude, weather_latitude, weather_cell
  )])
  weather <- weather[
    weather_cells,
    on = c(
      "Longitude" = "weather_longitude", "Latitude" = "weather_latitude"
    ),
    nomatch = 0L,
    allow.cartesian = TRUE
  ][Date %in% analysis_dates]
  weather_matrix <- function(variable) {
    result <- matrix(
      NA_real_, length(analysis_dates), max(weather_cells$weather_cell)
    )
    result[cbind(
      match(weather$Date, analysis_dates), weather$weather_cell
    )] <- weather[[variable]]
    result[, station_information$weather_cell, drop = FALSE]
  }
  temperature <- weather_matrix("Temperature")
  precipitation <- weather_matrix("Precipitation")
  wind_speed <- weather_matrix("Wind")
  if (!all(is.finite(c(temperature, precipitation, wind_speed)))) {
    stop("Weather matching produced missing covariates.")
  }
  time_major <- function(values) as.vector(t(values))
  X <- cbind(
    intercept = 1,
    temperature = time_major(
      (temperature - scaling$temperature_center) / scaling$temperature_scale
    ),
    precipitation = time_major(
      (precipitation - scaling$precipitation_center) /
        scaling$precipitation_scale
    ),
    wind_speed = time_major(
      (wind_speed - scaling$wind_speed_center) / scaling$wind_speed_scale
    )
  )

  latitude_reference <- mean(station_information$Latitude) * pi / 180
  locations <- cbind(
    x_1000km = station_information$Longitude * cos(latitude_reference) *
      111.32 / 1000,
    y_1000km = station_information$Latitude * 110.57 / 1000
  )
  rows_for_sites <- function(sites) {
    unlist(lapply(seq_len(nrow(Y)), function(time) {
      (time - 1L) * ncol(Y) + sites
    }), use.names = FALSE)
  }
  locations_training <- locations[training, , drop = FALSE]
  locations_holdout <- locations[holdout, , drop = FALSE]
  similarity_center <- colMeans(locations_training)
  similarity_scale <- apply(locations_training, 2L, stats::sd)
  similarity_training <- sweep(
    sweep(locations_training, 2L, similarity_center),
    2L, similarity_scale, "/"
  )
  similarity_holdout <- sweep(
    sweep(locations_holdout, 2L, similarity_center),
    2L, similarity_scale, "/"
  )

  list(
    Y = Y,
    Y_raw = Y_raw,
    X = X,
    locations = locations,
    station_information = station_information,
    analysis_dates = analysis_dates,
    scaling = scaling,
    training = training,
    holdout = holdout,
    Y_training = Y[, training, drop = FALSE],
    Y_holdout = Y[, holdout, drop = FALSE],
    X_training = X[rows_for_sites(training), , drop = FALSE],
    X_holdout = X[rows_for_sites(holdout), , drop = FALSE],
    locations_training = locations_training,
    locations_holdout = locations_holdout,
    similarity_training = similarity_training,
    similarity_holdout = similarity_holdout,
    NN_training = make_previous_neighbor_matrix(
      locations_training, n_neighbors
    ),
    prediction_neighbors = if (length(holdout)) {
      nearest_observed_neighbors(
        locations_training, locations_holdout, n_neighbors
      )
    } else {
      matrix(integer(), 0L, n_neighbors)
    },
    preparation = preparation,
    unique_sample = unique_sample,
    data_file = if (unique_sample) normalizePath(unique_data_file, winslash = "/") else NULL,
    data_md5 = if (unique_sample) unname(tools::md5sum(unique_data_file)) else NULL,
    original_holdout = original_holdout,
    holdout_result_name = holdout_result_name,
    preparation_name = preparation_name,
    n_neighbors = n_neighbors
  )
}

europe_initial_values <- function(prepared, model) {
  model <- validate_dcs_model_type(model)
  training <- prepared$training
  if (isTRUE(prepared$unique_sample)) {
    return(list(
      beta = rep(0, ncol(prepared$X)), xi2 = rep(1, ncol(prepared$X)),
      sigma2 = 0.2, phi = rep(0.5, length(training)),
      tau2 = rep(0.35, length(training)), M = 1, kappa = 1.5,
      sigma2_u = 0.3
    ))
  }
  initial <- prepared$preparation$initial_values
  if (is.null(initial)) stop("This preparation input has no historical warm start; use the unique-location dataset.")
  phi <- initial$phi[training]
  tau2 <- initial$tau2[training]
  if (model == "k1_common_shock") {
    phi[] <- mean(phi)
    tau2[] <- mean(tau2)
  }
  list(
    beta = initial$beta,
    xi2 = initial$xi2,
    sigma2 = initial$sigma2,
    phi = phi,
    tau2 = tau2,
    M = initial$M,
    kappa = initial$kappa,
    sigma2_u = initial$sigma2_u
  )
}

prepare_europe_full_prediction_inputs <- function(prepared, n_neighbors = 10L) {
  production_prediction <- prepared$preparation$prediction
  if (is.null(production_prediction) || is.null(production_prediction$grid)) {
    stop("The preparation input does not contain the manuscript prediction grid.")
  }
  grid <- as.data.frame(production_prediction$grid)
  if (!all(c("Longitude", "Latitude") %in% names(grid)) || !nrow(grid)) {
    stop("The stored manuscript prediction grid is invalid.")
  }

  station_information <- prepared$station_information
  latitude_reference <- mean(station_information$Latitude) * pi / 180
  prediction_locations <- cbind(
    x_1000km = grid$Longitude * cos(latitude_reference) * 111.32 / 1000,
    y_1000km = grid$Latitude * 110.57 / 1000
  )
  similarity_center <- colMeans(prepared$locations)
  similarity_scale <- apply(prepared$locations, 2L, stats::sd)
  prediction_similarity <- sweep(
    sweep(prediction_locations, 2L, similarity_center, "-"),
    2L, similarity_scale, "/"
  )
  prediction_neighbors <- nearest_observed_neighbors(
    prepared$locations, prediction_locations, n_neighbors
  )
  nearest_station <- prediction_neighbors[, 1L]
  stored_nearest <- production_prediction$nearest_observed_station
  if (!isTRUE(prepared$unique_sample) && !is.null(stored_nearest) &&
      !identical(as.integer(stored_nearest), as.integer(nearest_station))) {
    stop("Reconstructed nearest stations differ from the production grid.")
  }

  n_sites <- ncol(prepared$Y)
  prediction_rows <- unlist(lapply(seq_len(nrow(prepared$Y)), function(time) {
    (time - 1L) * n_sites + nearest_station
  }), use.names = FALSE)
  prediction_X <- prepared$X[prediction_rows, , drop = FALSE]

  list(
    grid = grid,
    locations = prediction_locations,
    similarity = prediction_similarity,
    neighbors = prediction_neighbors,
    X = prediction_X,
    nearest_observed_station = nearest_station
  )
}

normal_crps <- function(observation, mean, sd) {
  sd <- pmax(sd, sqrt(.Machine$double.eps))
  z <- (observation - mean) / sd
  sd * (z * (2 * stats::pnorm(z) - 1) +
    2 * stats::dnorm(z) - 1 / sqrt(pi))
}

europe_holdout_metrics <- function(prediction, prepared, model) {
  scaling <- prepared$scaling
  truth <- prepared$Y_raw[, prepared$holdout, drop = FALSE]
  mean_original <- prediction$y$mean * scaling$response_scale +
    scaling$response_center
  sd_original <- prediction$y$sd * scaling$response_scale
  lower <- mean_original + stats::qnorm(0.025) * sd_original
  upper <- mean_original + stats::qnorm(0.975) * sd_original
  error <- mean_original - truth
  absolute_error <- abs(error)
  crps <- normal_crps(truth, mean_original, sd_original)
  overall <- data.frame(
    model = model,
    model_label = dcs_model_label(model),
    RMSE = sqrt(mean(error^2)),
    MAE = mean(absolute_error),
    median_absolute_error = stats::median(absolute_error),
    bias = mean(error),
    correlation = stats::cor(as.vector(mean_original), as.vector(truth)),
    CRPS = mean(crps),
    coverage_95 = mean(truth >= lower & truth <= upper),
    mean_width_95 = mean(upper - lower),
    stringsAsFactors = FALSE
  )
  station <- do.call(rbind, lapply(seq_along(prepared$holdout), function(index) {
    data.frame(
      model = model,
      Samplingpoint = prepared$station_information$Samplingpoint[
        prepared$holdout[index]
      ],
      Longitude = prepared$station_information$Longitude[prepared$holdout[index]],
      Latitude = prepared$station_information$Latitude[prepared$holdout[index]],
      RMSE = sqrt(mean(error[, index]^2)),
      MAE = mean(absolute_error[, index]),
      median_absolute_error = stats::median(absolute_error[, index]),
      bias = mean(error[, index]),
      correlation = stats::cor(mean_original[, index], truth[, index]),
      CRPS = mean(crps[, index]),
      coverage_95 = mean(
        truth[, index] >= lower[, index] & truth[, index] <= upper[, index]
      ),
      mean_width_95 = mean(upper[, index] - lower[, index]),
      stringsAsFactors = FALSE
    )
  }))
  coordinate_group <- interaction(
    station$Longitude, station$Latitude, drop = TRUE
  )
  group <- do.call(rbind, lapply(split(station, coordinate_group), function(values) {
    data.frame(
      model = model,
      Longitude = values$Longitude[1L],
      Latitude = values$Latitude[1L],
      n_series = nrow(values),
      RMSE = sqrt(mean(values$RMSE^2)),
      MAE = mean(values$MAE),
      CRPS = mean(values$CRPS),
      coverage_95 = mean(values$coverage_95),
      mean_width_95 = mean(values$mean_width_95)
    )
  }))
  list(overall = overall, station = station, group = group)
}
