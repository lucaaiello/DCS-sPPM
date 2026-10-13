normalize_pm10_id <- function(values) {
  toupper(gsub(":", "_", sub("^[A-Z]{2}/", "", trimws(as.character(values))), fixed = TRUE))
}

parse_pm10_activity_date <- function(values) {
  values <- trimws(as.character(values))
  result <- as.Date(substr(values, 1L, 10L), format = "%d/%m/%Y")
  missing <- is.na(result)
  result[missing] <- as.Date(substr(values[missing], 1L, 10L), format = "%Y-%m-%d")
  result
}

pm10_activity_metadata <- function(metadata, start_date, end_date) {
  metadata <- data.table::copy(data.table::as.data.table(metadata))
  metadata[, normalized_id := normalize_pm10_id(`Sampling Point Id`)]
  metadata[, metadata_year := suppressWarnings(as.integer(Year))]
  metadata[, process_begin := parse_pm10_activity_date(`Process Activity Begin`)]
  metadata[, process_end := parse_pm10_activity_date(`Process Activity End`)]
  metadata[, station_begin := parse_pm10_activity_date(`Operational Activity Begin`)]
  metadata[, station_end := parse_pm10_activity_date(`Operational Activity End`)]
  metadata[, activity_compatible_with_2024 :=
    (is.na(process_begin) | process_begin <= end_date) &
    (is.na(process_end) | process_end >= start_date) &
    (is.na(station_begin) | station_begin <= end_date) &
    (is.na(station_end) | station_end >= start_date)]
  metadata[, process_start_documented := !is.na(process_begin)]
  metadata[, activity_dates_unparsed :=
    (!is.na(`Process Activity Begin`) & nzchar(trimws(`Process Activity Begin`)) & is.na(process_begin)) |
    (!is.na(`Process Activity End`) & nzchar(trimws(`Process Activity End`)) & is.na(process_end)) |
    (!is.na(`Operational Activity Begin`) & nzchar(trimws(`Operational Activity Begin`)) & is.na(station_begin)) |
    (!is.na(`Operational Activity End`) & nzchar(trimws(`Operational Activity End`)) & is.na(station_end))]
  metadata[, method_description := tolower(paste(
    `Measurement Method`, `Other Measurement Method`, `Analytical Technique`,
    `Other Analytical Technique`, `Sampling Method`, `Other Sampling Method`, `Process Id`
  ))]
  metadata[, activity_method_rank := data.table::fifelse(
    tolower(trimws(`Equivalence Demonstrated`)) == "ref" | grepl("gravi", method_description), 1L,
    data.table::fifelse(tolower(trimws(`Equivalence Demonstrated`)) == "yes", 2L,
      data.table::fifelse(tolower(trimws(`Equivalence Demonstrated`)) %in% c("no", "noref", "inprogress"), 4L, 3L))
  )]
  metadata
}

pm10_rankable_metadata <- function(metadata, analysis_year = 2024L) {
  active <- data.table::copy(metadata[activity_compatible_with_2024 == TRUE])
  active <- active[, {
    available_years <- metadata_year[!is.na(metadata_year)]
    chosen_year <- if (analysis_year %in% available_years) analysis_year else if (any(available_years <= analysis_year)) {
      max(available_years[available_years <= analysis_year])
    } else if (length(available_years)) min(available_years) else NA_integer_
    if (is.na(chosen_year)) .SD else .SD[metadata_year == chosen_year]
  }, by = normalized_id]
  active[, process_group := paste(normalized_id, `Process Id`, sep = "::")]
  active[, method_rank := {
    documented_ranks <- activity_method_rank[activity_method_rank != 3L]
    if (length(documented_ranks)) max(documented_ranks) else 3L
  }, by = process_group]
  active
}
