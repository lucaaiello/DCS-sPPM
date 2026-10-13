required <- c("data.table", "arrow")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Required packages: ", paste(missing, collapse = ", "))
library(data.table)
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_argument)) stop("Run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[1L]), winslash = "/")
project_dir <- dirname(script_path)
while (!file.exists(file.path(project_dir, "R", "project_paths.R"))) {
  parent <- dirname(project_dir)
  if (identical(parent, project_dir)) stop("Could not locate the manuscript project.")
  project_dir <- parent
}
repository_dir <- dirname(project_dir)
source(file.path(project_dir, "R", "unique_location_selection_helpers.R"))
preparation_dir <- file.path(project_dir, "Real_data_analysis", "preparation", "results", "unique_location_preparation")
resolution_dir <- file.path(preparation_dir, "exception_resolution")
dir.create(resolution_dir, recursive = TRUE, showWarnings = FALSE)
raw_dir <- file.path(resolution_dir, "source_parquet")
dir.create(raw_dir, showWarnings = FALSE)
stations <- fread(file.path(preparation_dir, "sampling_point_selection_audit.csv"))
previous_selection <- stations[selected_location == TRUE]
baseline_path <- file.path(resolution_dir, "selection_before_resolution.csv")
if (!file.exists(baseline_path)) fwrite(previous_selection, baseline_path)
observations <- as.data.table(readRDS(file.path(project_dir, "Data", "pm10_combined_2024.rds")))
observations[, Date := as.Date(as.numeric(unclass(Date)), origin = "1970-01-01")]
observations <- observations[Date >= as.Date("2024-01-01") & Date <= as.Date("2024-12-31") &
  Samplingpoint %in% stations$Samplingpoint, .(Samplingpoint, Date, Value)]
metadata <- pm10_activity_metadata(
  fread(file.path(repository_dir, "Data", "eea_pm10", "metadata", "DataExtract.csv"), encoding = "UTF-8"),
  as.Date("2024-01-01"), as.Date("2024-12-31")
)
active_metadata <- pm10_rankable_metadata(metadata)
catalogue <- fread(file.path(repository_dir, "Data", "eea_pm10", "eea_pm10_parquet_urls_day_hour.csv"))
catalogue[, normalized_id := normalize_pm10_id(sub("[.]parquet$", "", basename(ParquetFileUrl)))]
target_ids <- c(paste0("SPO_PL", c("0836A", "0839A", "0844A", "0846A"), "_5_002"),
  sort(unique(stations[country_code == "MK", sub("_[0-9]+$", "_100", normalized_id)])),
  "SPO-MK0043A_00005_101")
evidence <- list()
comparisons <- list()
verification <- list()
for (target_index in seq_along(target_ids)) {
  identifier <- target_ids[target_index]
  message("[", target_index, "/", length(target_ids), "] Checking original EEA record ", identifier)
  country <- if (startsWith(identifier, "SPO_PL")) "PL" else "MK"
  raw_path <- file.path(raw_dir, paste0(identifier, ".parquet"))
  source_reference <- NA_character_
  if (country == "MK") {
    zip_path <- file.path(repository_dir, "Data", "eea_pm10", "2024", "extracts", "MK_hour.zip")
    entries <- unzip(zip_path, list = TRUE)$Name
    member <- entries[normalize_pm10_id(sub("[.]parquet$", "", basename(entries))) == identifier]
    if (length(member) == 1L) {
      unzip(zip_path, files = member, exdir = raw_dir, junkpaths = TRUE, overwrite = TRUE)
      source_reference <- paste0(zip_path, "::", member)
    }
  }
  if (!file.exists(raw_path)) {
    urls <- unique(catalogue[normalized_id == identifier, ParquetFileUrl])
    if (!length(urls) && identifier == "SPO-MK0043A_00005_101") {
      urls <- "https://eeadmz1batchservice02.blob.core.windows.net/airquality-p-e1a/MK/SPO-MK0043A_00005_101.parquet"
    }
    if (length(urls) != 1L) {
      evidence[[length(evidence) + 1L]] <- data.table(normalized_id = identifier,
        status = "unresolved", evidence = "No unique local or catalogued source file")
      next
    }
    download.file(urls, raw_path, mode = "wb", quiet = TRUE, method = "libcurl")
    source_reference <- urls
  } else if (is.na(source_reference)) {
    urls <- unique(catalogue[normalized_id == identifier, ParquetFileUrl])
    source_reference <- if (length(urls) == 1L) urls else if (identifier == "SPO-MK0043A_00005_101") {
      "https://eeadmz1batchservice02.blob.core.windows.net/airquality-p-e1a/MK/SPO-MK0043A_00005_101.parquet"
    } else raw_path
  }
  raw <- as.data.table(arrow::read_parquet(raw_path))
  raw[, normalized_id := normalize_pm10_id(Samplingpoint)]
  raw[, Date := as.Date(Start, tz = "UTC")]
  raw <- raw[normalized_id == identifier & Date >= as.Date("2024-01-01") & Date <= as.Date("2024-12-31")]
  if ("Pollutant" %in% names(raw)) {
    raw <- raw[grepl("PM10|(^|/)5$", as.character(Pollutant), ignore.case = TRUE)]
  }
  if (any(!grepl("ug.m-3|µg/m|μg/m", as.character(raw$Unit)))) stop("Unexpected concentration units for ", identifier)
  raw[, valid_verified := as.integer(Validity) == 1L & as.integer(Verification) == 1L]
  raw[, finite_value := is.finite(as.numeric(Value))]
  native_daily <- raw[tolower(AggType) == "day" & finite_value == TRUE,
    .(raw_value = mean(as.numeric(Value)), raw_records = .N, all_valid_verified = all(valid_verified)), by = Date]
  hourly_daily <- raw[tolower(AggType) == "hour" & finite_value == TRUE,
    .(raw_value = mean(as.numeric(Value)), raw_records = .N, all_valid_verified = all(valid_verified)), by = Date]
  daily <- if (nrow(native_daily)) native_daily else hourly_daily
  raw_source_rank <- if (nrow(native_daily)) 1L else if (nrow(hourly_daily)) 2L else 3L
  valid_daily <- raw[valid_verified == TRUE & finite_value == TRUE,
    .(valid_value = mean(as.numeric(Value))), by = .(Date, AggType)]
  valid_daily <- valid_daily[tolower(AggType) == if (raw_source_rank == 1L) "day" else "hour"]
  daily <- merge(daily, valid_daily[, .(Date, valid_value)], by = "Date", all.x = TRUE)
  target_station_code <- if (country == "PL") sub("^SPO_(PL[0-9]+A).*$", "\\1", identifier) else {
    sub("^SPO-(MK[0-9]+A).*$", "\\1", identifier)
  }
  candidate_ids <- stations[station_code == target_station_code & country_code == country, Samplingpoint]
  for (candidate_id in candidate_ids) {
    stored <- observations[Samplingpoint == candidate_id & is.finite(Value)]
    paired <- merge(stored, daily, by = "Date", all.x = TRUE)
    paired[, matches_all := is.finite(raw_value) & abs(Value - raw_value) <= 1e-7 * pmax(1, abs(Value))]
    paired[, matches_valid := is.finite(valid_value) & abs(Value - valid_value) <= 1e-7 * pmax(1, abs(Value))]
    comparisons[[length(comparisons) + 1L]] <- paired[, .(
      source_identifier = identifier, stored_identifier = candidate_id, station_code = target_station_code,
      stored_days = .N, raw_overlapping_days = sum(is.finite(raw_value)),
      matching_days_all_finite = sum(matches_all), matching_days_valid_verified = sum(matches_valid),
      daily_values_including_unverified_records = sum(matches_all & !matches_valid),
      max_absolute_difference = if (any(is.finite(raw_value))) max(abs(Value - raw_value), na.rm = TRUE) else NA_real_,
      raw_source_rank = raw_source_rank, source_reference = source_reference
    )]
    verification[[length(verification) + 1L]] <- paired[, .(
      Samplingpoint = candidate_id, Date, Value, source_identifier = identifier,
      raw_value, valid_value, raw_records, all_valid_verified, matches_all, matches_valid
    )]
  }
  evidence[[length(evidence) + 1L]] <- data.table(normalized_id = identifier,
    status = "source checked", evidence = source_reference, md5 = unname(tools::md5sum(raw_path)))
}
pair_checks <- rbindlist(comparisons, fill = TRUE)
daily_checks <- rbindlist(verification, fill = TRUE)
fwrite(pair_checks, file.path(resolution_dir, "original_source_comparison.csv"))
fwrite(daily_checks, file.path(resolution_dir, "original_source_daily_checks.csv"))
fwrite(rbindlist(evidence, fill = TRUE), file.path(resolution_dir, "original_source_manifest.csv"))
decisions <- list()
for (target_station_key in unique(stations[country_code == "MK", station_key])) {
  candidates <- stations[station_key == target_station_key]
  exact <- candidates[normalized_id %in% active_metadata$normalized_id]
  if (!nrow(exact)) {
    candidate <- candidates[order(-observed_days, Samplingpoint)][1L]
    original_check <- pair_checks[stored_identifier == candidate$Samplingpoint & source_identifier == candidate$normalized_id]
    station_records <- active_metadata[`Air Quality Station EoI Code` == candidate$station_code]
    station_coordinates <- unique(station_records[, .(Longitude, Latitude)])
    coordinates_agree <- nrow(station_coordinates) == 1L &&
      abs(station_coordinates$Longitude - candidate$Longitude) < 1e-5 &&
      abs(station_coordinates$Latitude - candidate$Latitude) < 1e-5
    station_supported <- nrow(original_check) == 1L &&
      original_check$matching_days_valid_verified >= 0.99 * candidate$observed_days && coordinates_agree
    decisions[[length(decisions) + 1L]] <- data.table(station_key = target_station_key,
      Samplingpoint = candidate$Samplingpoint,
      metadata_sampling_point = NA_character_, Longitude = NA_real_, Latitude = NA_real_,
      source_rank_verified = if (station_supported) original_check$raw_source_rank else NA_integer_,
      decision = if (station_supported) "retain" else "unresolved",
      metadata_scope = "station_only",
      evidence_source = if (nrow(original_check)) original_check$source_reference else "Local EEA metadata",
      rationale = if (station_supported) {
        paste0("Original sampling-point file reproduces ", original_check$matching_days_valid_verified,
          "/", candidate$observed_days, " stored days under validity filtering; station code and coordinates agree with EEA metadata. Instrument method remains unknown; unmatched day is recorded, not altered.")
      } else "No exact metadata match and insufficient original-source/coordinate evidence")
    next
  }
  chosen <- exact[order(-observed_days, Samplingpoint)][1L]
  paired_checks <- pair_checks[stored_identifier == chosen$Samplingpoint & source_identifier == chosen$normalized_id]
  support <- nrow(paired_checks) && max(paired_checks$matching_days_all_finite[1L],
    paired_checks$matching_days_valid_verified[1L]) >= 0.90 * chosen$observed_days
  decisions[[length(decisions) + 1L]] <- data.table(station_key = target_station_key,
    Samplingpoint = chosen$Samplingpoint, metadata_sampling_point = chosen$Samplingpoint,
    Longitude = NA_real_, Latitude = NA_real_,
    source_rank_verified = if (support) paired_checks$raw_source_rank[1L] else NA_integer_,
    decision = "retain", evidence_source = "Exact identifier in EEA metadata and original MK_hour extract",
    rationale = if (support) "Prefer the metadata-backed _100 series; raw-source agreement checked; identical duplicate partners are redundant" else {
      "Prefer exact metadata match; raw-source coverage differences remain documented"
    })
}
for (identifier in target_ids[startsWith(target_ids, "SPO_PL")]) {
  candidate <- stations[normalized_id == identifier]
  if (!nrow(candidate)) next
  checked <- pair_checks[stored_identifier == candidate$Samplingpoint & source_identifier == identifier]
  confirmed <- nrow(checked) && checked$matching_days_all_finite == candidate$observed_days && checked$raw_source_rank == 1L
  decisions[[length(decisions) + 1L]] <- data.table(station_key = candidate$station_key,
    Samplingpoint = candidate$Samplingpoint, metadata_sampling_point = candidate$Samplingpoint,
    Longitude = NA_real_, Latitude = NA_real_, source_rank_verified = if (confirmed) 1L else NA_integer_,
    decision = if (confirmed) "retain" else "unresolved",
    evidence_source = if (nrow(checked)) checked$source_reference else "Original EEA source unavailable",
    rationale = if (confirmed) "Native daily original file reproduces every stored finite value; 2024-active reference procedure documented" else {
      "Stored values do not all match the original native daily file; inspect original_source_daily_checks.csv"
    })
}
for (code in c("SI0038A", "LU0105A", "PL0052A")) {
  candidates <- stations[station_code == code]
  record_candidates <- active_metadata[`Air Quality Station EoI Code` == code &
    normalized_id %in% candidates$normalized_id]
  if (!nrow(record_candidates)) next
  reference <- record_candidates[method_rank == min(method_rank)]
  reference_ids <- unique(reference$normalized_id)
  ranked <- candidates[normalized_id %in% reference_ids][order(-observed_days, Samplingpoint)]
  chosen <- ranked[1L]
  coordinates <- unique(reference[normalized_id == chosen$normalized_id, .(Longitude, Latitude)])
  coordinates <- coordinates[is.finite(Longitude) & is.finite(Latitude)]
  decisions[[length(decisions) + 1L]] <- data.table(station_key = chosen$station_key,
    Samplingpoint = chosen$Samplingpoint, metadata_sampling_point = chosen$Samplingpoint,
    Longitude = if (nrow(coordinates) == 1L) coordinates$Longitude else NA_real_,
    Latitude = if (nrow(coordinates) == 1L) coordinates$Latitude else NA_real_,
    source_rank_verified = NA_integer_, decision = if (nrow(coordinates) == 1L) "retain" else "unresolved",
    evidence_source = "Data/eea_pm10/metadata/DataExtract.csv: activity-compatible reference procedure",
    rationale = "One physical EEA station; choose its reference series and the unique coordinates attached to the procedure active in 2024"
  )
}
decision_log <- rbindlist(decisions, fill = TRUE)
fwrite(decision_log, file.path(preparation_dir, "series_selection_overrides.csv"))
decision_codes <- sub("^[A-Z]{2}/", "", decision_log$station_key)
fwrite(metadata[`Air Quality Station EoI Code` %in% decision_codes],
  file.path(resolution_dir, "exception_metadata_records.csv"))
writeLines(c(
  "Unique-location exception resolution",
  paste("Documented decisions:", nrow(decision_log)),
  paste("Unresolved decisions:", sum(decision_log$decision == "unresolved")),
  "Original source comparison does not certify validity: both unfiltered and valid/verified comparisons are reported.",
  "Coordinates are taken from the reference procedure whose documented activity overlaps 2024.",
  "One series per physical EEA station is the selection policy, including where instrument coordinates differ.",
  "Raw data and fitted analysis results are preserved. Run scripts 54 and 55 after this script."
), file.path(resolution_dir, "resolution_summary.txt"))
print(decision_log[, .(station_key, Samplingpoint, Longitude, Latitude, decision)])
