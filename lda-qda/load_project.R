# Source this file to load the project-local packages and analysis functions.
# No analysis or output is started by loading it.
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  locale_result <- Sys.setlocale("LC_CTYPE", "English_United States.utf8")
  if (!nzchar(locale_result)) stop("Cannot enable UTF-8 for this R session.")
}
analysis_project_dir <- local({
  paths <- lapply(sys.frames(), function(fr) fr$ofile)
  paths <- Filter(function(z) is.character(z) && length(z) == 1L, paths)
  if (!length(paths)) stop("Load this file using source().")
  normalizePath(dirname(tail(paths, 1L)[[1L]]), winslash = "/", mustWork = TRUE)
})
project_library <- file.path(analysis_project_dir, ".r-library")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))
source(file.path(analysis_project_dir, "scripts", "run_analysis.R"),
       encoding = "UTF-8")

