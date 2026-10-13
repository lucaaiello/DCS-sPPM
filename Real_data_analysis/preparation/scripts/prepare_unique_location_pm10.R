required_packages <- "data.table"
if (!requireNamespace(required_packages, quietly = TRUE)) {
  stop("Install data.table before running this script.")
}
library(data.table)

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_argument)) {
  stop("Run this standalone script with Rscript.")
}
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]), winslash = "/")
project_dir <- dirname(script_path)
while (!file.exists(file.path(project_dir, "R", "project_paths.R"))) {
  parent <- dirname(project_dir)
  if (identical(parent, project_dir)) stop("Could not locate the manuscript project.")
  project_dir <- parent
}
repository_dir <- dirname(project_dir)
source(file.path(project_dir, "R", "unique_location_selection_helpers.R"))
preparation_path <- file.path(project_dir, "Data", "production_30000_online_grid_preparation.rds")
data_path <- file.path(project_dir, "Data", "pm10_combined_2024.rds")
metadata_path <- file.path(repository_dir, "Data", "eea_pm10", "metadata", "DataExtract.csv")
source_paths <- list.files(file.path(repository_dir, "Data", "eea_pm10", "2024", "daily_chunks"),
  pattern = "[.]rds$", full.names = TRUE)
output_dir <- file.path(project_dir, "Real_data_analysis", "preparation", "results", "unique_location_preparation")
output_path <- file.path(project_dir, "Data", "pm10_unique_locations_2024.rds")
override_path <- file.path(output_dir, "series_selection_overrides.csv")
minimum_coverage <- 0.90
input_paths <- c(preparation_path, data_path, metadata_path, source_paths)
overrides <- if (file.exists(override_path)) fread(override_path) else data.table()
if (nrow(overrides)) input_paths <- c(input_paths, override_path)
if (!length(source_paths)) stop("Existing country-level daily caches were not found.")
if (any(!file.exists(input_paths))) {
  stop("Missing input(s): ", paste(input_paths[!file.exists(input_paths)], collapse = ", "))
}

normalize_sampling_point <- function(values) {
  toupper(gsub(":", "_", sub("^[A-Z]{2}/", "", trimws(as.character(values))), fixed = TRUE))
}
collapse_values <- function(values) {
  paste(sort(unique(values[!is.na(values) & nzchar(values)])), collapse = "; ")
}
longest_missing_run <- function(observed_indices, number_of_days) {
  observed <- rep(FALSE, number_of_days)
  observed[observed_indices] <- TRUE
  runs <- rle(!observed)
  if (any(runs$values)) max(runs$lengths[runs$values]) else 0L
}

message("Reading the final analysis cohort and existing daily observations.")
preparation <- readRDS(preparation_path)
stations <- as.data.table(preparation$station_information)[, .(Samplingpoint, Longitude, Latitude)]
analysis_dates <- as.Date(preparation$analysis_dates)
study_bbox <- preparation$study_bbox
rm(preparation)
invisible(gc())
if (anyDuplicated(stations$Samplingpoint) || any(!is.finite(stations$Longitude) | !is.finite(stations$Latitude))) {
  stop("The final analysis cohort contains invalid identifiers or coordinates.")
}
stations[, normalized_id := normalize_sampling_point(Samplingpoint)]
stations[, country_code := sub("/.*$", "", Samplingpoint)]
stations[, `:=`(original_longitude = Longitude, original_latitude = Latitude,
  metadata_lookup_id = normalized_id, decision_applied = FALSE, resolution_evidence = "",
  station_identity_verified = FALSE)]
approved_overrides <- if (nrow(overrides)) overrides[decision == "retain"] else data.table()
if (nrow(approved_overrides)) {
  if (anyDuplicated(approved_overrides$station_key) || any(!approved_overrides$Samplingpoint %in% stations$Samplingpoint)) {
    stop("Overrides must identify one existing sampling point per station key.")
  }
  for (override_index in seq_len(nrow(approved_overrides))) {
    override <- approved_overrides[override_index]
    candidate_index <- match(override$Samplingpoint, stations$Samplingpoint)
    if (!is.na(override$metadata_sampling_point) && nzchar(override$metadata_sampling_point)) {
      stations$metadata_lookup_id[candidate_index] <- normalize_sampling_point(override$metadata_sampling_point)
    }
    if (xor(is.finite(override$Longitude), is.finite(override$Latitude))) stop("Incomplete override coordinates.")
    if (is.finite(override$Longitude) && is.finite(override$Latitude)) {
      stations$Longitude[candidate_index] <- override$Longitude
      stations$Latitude[candidate_index] <- override$Latitude
    }
    stations$decision_applied[candidate_index] <- TRUE
    stations$station_identity_verified[candidate_index] <- TRUE
    stations$resolution_evidence[candidate_index] <- paste(override$evidence_source, override$rationale, sep = "; ")
  }
}
observations <- as.data.table(readRDS(data_path))
observations[, Date := as.Date(as.numeric(unclass(Date)), origin = "1970-01-01")]
observations <- observations[Samplingpoint %in% stations$Samplingpoint & Date %in% analysis_dates]
if (anyDuplicated(observations[, .(Samplingpoint, Date)])) {
  stop("Duplicate sampling-point/date records in the existing daily dataset.")
}
completeness <- observations[, .(
  observed_days = sum(is.finite(Value)),
  longest_missing_run = longest_missing_run(match(Date[is.finite(Value)], analysis_dates), length(analysis_dates))
), by = Samplingpoint]
stations <- merge(stations, completeness, by = "Samplingpoint", all.x = TRUE)
stations[is.na(observed_days), `:=`(observed_days = 0L, longest_missing_run = length(analysis_dates))]
stations[, coverage := observed_days / length(analysis_dates)]

metadata_all <- pm10_activity_metadata(fread(metadata_path, encoding = "UTF-8", check.names = FALSE),
  min(analysis_dates), max(analysis_dates))
metadata_all <- metadata_all[normalized_id %in% stations$metadata_lookup_id]
metadata <- pm10_rankable_metadata(metadata_all)
metadata_summary <- metadata[, .(
  station_code = collapse_values(`Air Quality Station EoI Code`),
  station_name = collapse_values(`Air Quality Station Name`),
  measurement_type = collapse_values(`Measurement Type`),
  measurement_method = collapse_values(`Measurement Method`),
  analytical_technique = collapse_values(`Analytical Technique`),
  process_ids = collapse_values(`Process Id`),
  equipment = collapse_values(`Measurement Equipment`),
  equivalence = collapse_values(`Equivalence Demonstrated`),
  metadata_year = max(metadata_year, na.rm = TRUE),
  method_rank = max(method_rank),
  method_rank_conflict = uniqueN(method_rank) > 1L,
  activity_dates_unparsed = any(activity_dates_unparsed),
  activity_start_documented = any(process_start_documented),
  station_area = collapse_values(`Air Quality Station Area`),
  station_type = collapse_values(`Air Quality Station Type`)
), by = normalized_id]
setnames(metadata_summary, "normalized_id", "metadata_lookup_id")
stations <- merge(stations, metadata_summary, by = "metadata_lookup_id", all.x = TRUE)
stations[, metadata_matched := !is.na(station_code) & nzchar(station_code)]
stations[, station_code_inferred := FALSE]
unmatched_indices <- which(!stations$metadata_matched)
for (station_index in unmatched_indices) {
  identifier <- stations$normalized_id[station_index]
  code_match <- regmatches(identifier, regexpr("[A-Z]{2}[0-9]{4}A", identifier))
  if (length(code_match) && nzchar(code_match)) {
    stations$station_code[station_index] <- code_match
    stations$station_code_inferred[station_index] <- TRUE
  } else {
    stations$station_code[station_index] <- paste0("UNRESOLVED:", identifier)
  }
}
stations[is.na(method_rank), `:=`(method_rank = 3L, method_rank_conflict = FALSE)]
stations[, station_key := paste(country_code, station_code, sep = "/")]

message("Checking native-daily versus hourly-derived provenance against the stored observations.")
sources <- rbindlist(lapply(source_paths, function(path) as.data.table(readRDS(path))), fill = TRUE)
sources[, normalized_id := normalize_sampling_point(Samplingpoint)]
sources[, Date := as.Date(date)]
sources <- sources[normalized_id %in% stations$normalized_id & Date %in% analysis_dates]
source_candidates <- unique(sources[, .(normalized_id, Date, source, source_value = value)])
comparison <- merge(
  observations[, .(normalized_id = normalize_sampling_point(Samplingpoint), Date, Value)],
  source_candidates, by = c("normalized_id", "Date"), allow.cartesian = TRUE
)
comparison[, matches := is.finite(Value) & is.finite(source_value) &
  abs(Value - source_value) <= 1e-7 * pmax(1, abs(Value))]
source_summary <- comparison[, .(
  native_daily_matching_days = uniqueN(Date[matches & source == "day"]),
  hourly_derived_matching_days = uniqueN(Date[matches & source == "hour"]),
  provenance_matching_days = uniqueN(Date[matches]),
  available_sources = collapse_values(source)
), by = normalized_id]
stations <- merge(stations, source_summary, by = "normalized_id", all.x = TRUE)
count_columns <- c("native_daily_matching_days", "hourly_derived_matching_days", "provenance_matching_days")
for (column in count_columns) set(stations, which(is.na(stations[[column]])), column, 0L)
stations[, source_rank := fifelse(
  native_daily_matching_days >= 0.90 * observed_days & hourly_derived_matching_days == 0L & observed_days > 0L, 1L,
  fifelse(hourly_derived_matching_days >= 0.90 * observed_days & native_daily_matching_days == 0L & observed_days > 0L, 2L, 3L)
)]
stations[, source_rank_verified := FALSE]
if (nrow(approved_overrides)) {
  for (override_index in seq_len(nrow(approved_overrides))) {
    override <- approved_overrides[override_index]
    if (!is.na(override$source_rank_verified)) {
      if (!override$source_rank_verified %in% c(1L, 2L)) stop("Invalid source-rank override.")
      stations[Samplingpoint == override$Samplingpoint,
        `:=`(source_rank = override$source_rank_verified, source_rank_verified = TRUE)]
    }
  }
}
stations[, eligible := coverage >= minimum_coverage]
stations[, selected_station := FALSE]
stations[, selected_location := FALSE]
stations[, station_representative := NA_character_]
stations[, location_representative := NA_character_]
stations[, selection_reason := "Below 90% finite-observation coverage"]

stations[, override_rank := fifelse(decision_applied, 0L, 1L)]
ranking_columns <- c("override_rank", "method_rank", "source_rank", "observed_days", "longest_missing_run", "Samplingpoint")
ranking_order <- c(1L, 1L, 1L, -1L, 1L, 1L)
eligible <- copy(stations[eligible == TRUE])
setorderv(eligible, c("station_key", ranking_columns), c(1L, ranking_order))
station_winners <- eligible[, .(station_representative = Samplingpoint[1L]), by = station_key]
stations[station_winners, on = "station_key", station_representative := i.station_representative]
stations[, selected_station := eligible & Samplingpoint == station_representative]
station_candidates <- copy(stations[selected_station == TRUE])
setorderv(station_candidates, c("Longitude", "Latitude", ranking_columns), c(1L, 1L, ranking_order))
location_winners <- station_candidates[, .(location_representative = Samplingpoint[1L]), by = .(Longitude, Latitude)]
stations[location_winners, on = c("Longitude", "Latitude"), location_representative := i.location_representative]
stations[, selected_location := selected_station & Samplingpoint == location_representative]
stations[eligible == TRUE & selected_station == FALSE,
  selection_reason := paste0("Lower-ranked series for EEA station; retained ", station_representative)]
stations[selected_station == TRUE & selected_location == FALSE,
  selection_reason := paste0("Different station code at identical coordinates; retained ", location_representative)]
stations[selected_location == TRUE, selection_reason := "Highest-ranked eligible series at this unique location"]
stations[selected_location == TRUE & decision_applied == TRUE,
  selection_reason := paste0("Documented exception decision: ", resolution_evidence)]
stations[, station_coordinate_conflict := uniqueN(paste(Longitude, Latitude)) > 1L, by = station_key]
stations[, review_required := !metadata_matched | is.na(metadata_year) | metadata_year != 2024L |
  method_rank_conflict | station_coordinate_conflict |
  provenance_matching_days < observed_days | grepl(";", station_code, fixed = TRUE)]
stations[, review_reason := paste0(
  fifelse(!metadata_matched, "No direct metadata match; ", ""),
  fifelse(is.na(metadata_year) | metadata_year != 2024L, "Metadata year differs from 2024; ", ""),
  fifelse(method_rank_conflict, "Multiple procedure ranks; ", ""),
  fifelse(station_coordinate_conflict, "Station code reported at multiple coordinates; ", ""),
  fifelse(provenance_matching_days < observed_days, "Some values unmatched in existing source caches; ", "")
)]
stations[, original_coordinate_series := .N, by = .(Longitude, Latitude)]
stations[, station_series := .N, by = station_key]
if (nrow(approved_overrides) && any(!approved_overrides$Samplingpoint %in% stations[selected_location == TRUE]$Samplingpoint)) {
  stop("An approved override conflicts with coverage or another exact-coordinate selection; inspect decisions.")
}

selected_stations <- copy(stations[selected_location == TRUE])
setorder(selected_stations, Samplingpoint)
if (anyDuplicated(selected_stations[, .(Longitude, Latitude)]) || anyDuplicated(selected_stations$station_key)) {
  stop("Selection did not produce unique locations and unique station keys.")
}
selected_observations <- observations[Samplingpoint %in% selected_stations$Samplingpoint]
setorder(selected_observations, Samplingpoint, Date)
Y_observed <- matrix(NA_real_, length(analysis_dates), nrow(selected_stations),
  dimnames = list(as.character(analysis_dates), selected_stations$Samplingpoint))
Y_observed[cbind(match(selected_observations$Date, analysis_dates),
  match(selected_observations$Samplingpoint, selected_stations$Samplingpoint))] <- selected_observations$Value
Y_observed[!is.finite(Y_observed)] <- NA_real_
Y_imputed <- vapply(seq_len(ncol(Y_observed)), function(station_index) {
  values <- Y_observed[, station_index]
  observed <- which(is.finite(values))
  if (length(observed) < 2L) stop("Insufficient values for interpolation.")
  approx(observed, values[observed], xout = seq_along(values), rule = 2)$y
}, numeric(nrow(Y_observed)))
dimnames(Y_imputed) <- dimnames(Y_observed)
stopifnot(all(is.finite(Y_imputed)), all(Y_imputed[is.finite(Y_observed)] == Y_observed[is.finite(Y_observed)]))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setorder(stations, country_code, station_key, method_rank, source_rank, Samplingpoint)
fwrite(stations, file.path(output_dir, "sampling_point_selection_audit.csv"))
fwrite(selected_stations, file.path(output_dir, "selected_unique_locations.csv"))
fwrite(stations[review_required == TRUE], file.path(output_dir, "metadata_and_provenance_review.csv"))
fwrite(stations[original_coordinate_series > 1L | station_series > 1L],
  file.path(output_dir, "multiple_series_selection_audit.csv"))
fwrite(stations[station_coordinate_conflict == TRUE],
  file.path(output_dir, "station_coordinate_conflicts.csv"))
fwrite(metadata_all[activity_compatible_with_2024 == FALSE],
  file.path(output_dir, "excluded_historical_procedures.csv"))
baseline_path <- file.path(output_dir, "exception_resolution", "selection_before_resolution.csv")
if (file.exists(baseline_path)) {
  baseline <- fread(baseline_path)[, .(station_key, previous_series = Samplingpoint,
    previous_longitude = Longitude, previous_latitude = Latitude)]
  updated <- selected_stations[, .(station_key, updated_series = Samplingpoint,
    updated_longitude = Longitude, updated_latitude = Latitude, selection_reason)]
  changes <- merge(baseline, updated, by = "station_key", all = TRUE)
  changes <- changes[is.na(previous_series) | is.na(updated_series) | previous_series != updated_series |
    previous_longitude != updated_longitude | previous_latitude != updated_latitude]
  fwrite(changes, file.path(output_dir, "selection_changes_after_resolution.csv"))
}
country_summary <- stations[, .(
  original_series = .N,
  selected_locations = sum(selected_location),
  discarded_series = sum(!selected_location),
  retained_requiring_review = sum(selected_location & review_required)
), by = country_code][order(country_code)]
fwrite(country_summary, file.path(output_dir, "selection_by_country.csv"))
settings <- list(
  minimum_coverage = minimum_coverage,
  input_paths = input_paths,
  input_md5 = tools::md5sum(input_paths),
  study_bbox = study_bbox,
  ranking = c("reference or documented gravimetric", "equivalent", "unknown", "no equivalence or in progress"),
  tie_breaking = c("native daily", "hourly-derived daily", "unknown or mixed provenance", "more finite days", "shorter longest missing run", "identifier"),
  provenance_policy = "Source class requires at least 90% of observed values matching one source class and no matches to the other; unmatched days remain flagged",
  metadata_policy = "Filter procedures by activity overlap with 2024, then prefer 2024 metadata; fallback year flagged; ignore blank equivalence duplicates within a documented process, conservatively rank conflicting documented procedures",
  cohort_policy = "Candidates restricted to the exact final production analysis cohort",
  imputation = "Linear interpolation after selection, with nearest-endpoint extension",
  units = "PM10 concentrations retained in the original dataset units; no standardization"
)
saveRDS(list(
  observations = selected_observations,
  station_information = selected_stations,
  analysis_dates = analysis_dates,
  Y_observed = Y_observed,
  Y_imputed = Y_imputed,
  imputed_mask = is.na(Y_observed),
  selection_audit = stations,
  settings = settings
), output_path)
report <- c(
  "Unique-location PM10 preparation audit",
  paste("Original final-analysis series:", nrow(stations)),
  paste("Original distinct coordinate locations:", uniqueN(stations[, .(Longitude, Latitude)])),
  paste("Selected unique locations:", nrow(selected_stations)),
  paste("Discarded series:", sum(!stations$selected_location)),
  paste("Retained series requiring metadata/provenance review:", sum(selected_stations$review_required)),
  paste("Station codes reported at multiple coordinates:", uniqueN(stations[station_coordinate_conflict == TRUE]$station_key)),
  paste("Selected series using 2025 metadata fallback:", sum(selected_stations$metadata_year == 2025L, na.rm = TRUE)),
  paste("Selected series without a direct metadata match:", sum(!selected_stations$metadata_matched)),
  paste("Applied documented selection decisions:", sum(selected_stations$decision_applied)),
  paste("Excluded historical/future metadata records:", sum(!metadata_all$activity_compatible_with_2024)),
  paste("Missing daily cells before interpolation:", sum(is.na(Y_observed))),
  paste("Output dataset:", output_path),
  "",
  "The saved dataset is a named list; downstream scripts must use Y_imputed and station_information together.",
  "Finite observed values are preserved exactly. Imputation happens only after series selection.",
  "The retained cohort inherits the final production geographic support and recovered coordinates.",
  "Reference status comes from Equivalence Demonstrated=ref or explicit gravimetric metadata.",
  "Active measurement alone is not treated as proof of a reference method.",
  "Unmatched metadata, conflicting procedure ranks, and unmatched provenance are flagged for review.",
  "Metadata fallback and imperfect cache matching are audit flags, not automatic evidence of invalid observations.",
  "Selection describes method metadata, not a judgment that the discarded series are invalid.",
  "No sampler or fitted analysis is run by this script."
)
writeLines(report, file.path(output_dir, "preparation_summary.txt"))
cat(paste(report, collapse = "\n"), "\n")
print(country_summary)
