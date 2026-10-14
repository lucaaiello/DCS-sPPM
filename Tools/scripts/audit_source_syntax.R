script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("Run this read-only audit with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
project <- dirname(script_path)
while (!file.exists(file.path(project, "R", "project_paths.R"))) {
  parent <- dirname(project)
  if (identical(parent, project)) stop("Could not locate the manuscript project.")
  project <- parent
}
source_roots <- c("R", "Simulations", "Real_data_analysis", "Manuscript_figures", "Tools")
files <- unlist(lapply(source_roots, function(directory) {
  list.files(file.path(project, directory), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
}), use.names = FALSE)
files <- files[!grepl("/results/", files, fixed = TRUE)]
for (path in files) parse(file = path, encoding = "UTF-8")
cat("R source files parsed without executing them:", length(files), "\n")
