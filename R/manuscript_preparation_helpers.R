preparation_project <- function() {
  argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(argument)) stop("Run with Rscript.")
  current <- dirname(normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/"))
  while (!file.exists(file.path(current, "R", "project_paths.R"))) {
    parent <- dirname(current)
    if (identical(parent, current)) stop("Cannot locate project.")
    current <- parent
  }
  current
}

preparation_arguments <- function(defaults) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) %% 2L) stop("Arguments must be --name value pairs.")
  for (index in if (length(arguments)) seq.int(1L, length(arguments), by = 2L) else integer()) {
    name <- sub("^--", "", arguments[index])
    if (!startsWith(arguments[index], "--") || !name %in% names(defaults)) stop("Unknown argument: ", arguments[index])
    defaults[[name]] <- arguments[index + 1L]
  }
  defaults
}

preparation_fingerprint <- function(values) {
  temporary <- tempfile()
  on.exit(unlink(temporary))
  lines <- ifelse(is.na(values), "NA", sprintf("%.12g", as.numeric(values)))
  writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), "\n")), temporary)
  unname(tools::md5sum(temporary))
}

preparation_write_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path, compress = "xz")
}

preparation_day <- function(values) {
  if (inherits(values, "Date")) return(as.Date(as.numeric(values), origin = "1970-01-01"))
  if (inherits(values, "POSIXt")) return(as.Date(values, tz = "UTC"))
  if (is.numeric(values)) return(as.Date("2024-01-01") + values)
  as.Date(substr(as.character(values), 1L, 10L))
}
