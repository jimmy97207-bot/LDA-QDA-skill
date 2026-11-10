# Run with Rscript tests/validate.R or source("tests/validate.R").
if (.Platform$OS.type == "windows")
  invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
entry <- local({
  paths <- lapply(sys.frames(), function(fr) fr$ofile)
  paths <- Filter(function(z) is.character(z) && length(z) == 1L, paths)
  if (length(paths)) tail(paths, 1L)[[1L]] else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(arg)) stop("Cannot locate the validation script.")
    sub("^--file=", "", arg[[1L]])
  }
})
project <- normalizePath(file.path(dirname(entry), ".."),
                         winslash = "/", mustWork = TRUE)
source(file.path(project, "load_project.R"), encoding = "UTF-8")
has_mvn <- requireNamespace("MVN", quietly = TRUE)
has_biotools <- requireNamespace("biotools", quietly = TRUE)
parent <- file.path(project, "validation")
dir.create(parent, showWarnings = FALSE)
run_dir <- tempfile(pattern = paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-"),
                    tmpdir = parent)
stopifnot(dir.create(run_dir))
count <- 0L
check <- function(...) {
  stopifnot(...)
  count <<- count + 1L
}
expect_error <- function(expr) {
  outcome <- try(force(expr), silent = TRUE)
  if (!inherits(outcome, "try-error"))
    stop("Expected an error: ", paste(deparse(substitute(expr)), collapse = " "))
  check(TRUE)
}
x <- names(iris)[1:4]
set.seed(742)
rng <- .Random.seed
result <- run_pipeline(iris, "Species", x)
check(identical(rng, .Random.seed))
check(result$status == "completed",
      result$model_status == "completed",
      result$diagnostic_status == "completed",
      all(result$diagnostics$status == "computed"))
# Repeated CV: 10 repeats x 150 rows x 2 models.
check(result$settings$repeats == 10L, result$settings$folds == 5L,
      nrow(result$cv_predictions) == 3000L,
      nrow(result$cv_probabilities) == 9000L,
      all(result$metrics$coverage == 1),
      nrow(result$repeat_metrics) == 20L,
      nrow(result$fold_assignments) == 1500L)
# Repeat 1 equals the single split of the previous version (147 / 146 correct).
r1 <- subset(result$repeat_metrics, repeat_id == 1)
check(identical(r1$correct, c(147L, 146L)),
      isTRUE(all.equal(r1$macro_recall, c(0.98, 146/150))))
# Reported means are the means of the per-repeat values.
for (m in c("lda", "qda")) {
  z <- subset(result$repeat_metrics, model == m)
  mm <- subset(result$metrics, model == m)
  check(isTRUE(all.equal(mm$macro_recall, mean(z$macro_recall))),
        isTRUE(all.equal(mm$correct_mean, mean(z$correct))),
        mm$correct_min == min(z$correct), mm$correct_max == max(z$correct))
}
check(isTRUE(all.equal(result$metrics$correct_mean, c(146.9, 146.0))))
# Every repeat is a different stratified partition with 10 rows per class per fold.
fa <- result$fold_assignments
check(length(unique(lapply(split(fa$fold, fa$repeat_id), identity))) == 10L,
      all(table(fa$repeat_id, fa$fold, rep(iris$Species, 10)) == 10L))
# Pooled confusion matrix: counts sum to 10 x 150; row percentages sum to 1.
cm <- result$confusion_matrix
check(all(tapply(cm$Freq, cm$model, sum) == 1500L),
      all(abs(tapply(cm$row_pct, list(cm$model, cm$truth), sum) - 1) < 1e-12))
again <- run_pipeline(iris, "Species", x)
check(identical(result$fold_assignments, again$fold_assignments),
      identical(result$metrics, again$metrics))
psum <- aggregate(probability ~ model + repeat_id + row_id,
                  result$cv_probabilities, sum)
check(all(abs(psum$probability - 1) < 1e-8))
# Internal Mardia / Box's M against the values MVN 6.3 and biotools 4.3
# produced for iris on 2026-11-05 (R 4.6.0, Windows). When the packages are
# installed, they are also compared directly.
reference_diag <- data.frame(
  test = c(rep(c("Mardia Skewness", "Mardia Kurtosis"), 3), "Box M"),
  group = c(rep(c("setosa", "versicolor", "virginica"), each = 2), "ALL"),
  statistic = c(25.6643445196298, 1.29499223711605, 25.1850115362466,
                -0.57186635893429, 26.2705981752915, 0.152614173978342,
                140.943049923498),
  p_value = c(0.177185884467652, 0.195322907441935, 0.194444483140265,
              0.567412516528727, 0.157059707690356, 0.878702546726567,
              3.35203417831723e-20))
check(identical(result$diagnostics$test, reference_diag$test),
      identical(result$diagnostics$group, reference_diag$group),
      isTRUE(all.equal(result$diagnostics$statistic, reference_diag$statistic,
                       tolerance = 1e-10)),
      isTRUE(all.equal(result$diagnostics$p_value, reference_diag$p_value,
                       tolerance = 1e-10)))
if (has_mvn) for (cl in levels(iris$Species)) {
  direct <- MVN::mardia(iris[iris$Species == cl, x])
  got <- subset(result$diagnostics, group == cl)
  check(isTRUE(all.equal(got$statistic, direct$Statistic)),
        isTRUE(all.equal(got$p_value, direct$p.value)))
}
if (has_biotools) {
  direct_box <- biotools::boxM(iris[x], iris$Species)
  got_box <- subset(result$diagnostics, test == "Box M")
  check(isTRUE(all.equal(got_box$statistic, unname(direct_box$statistic))),
        isTRUE(all.equal(got_box$p_value, direct_box$p.value)))
}
# Fold-level metrics and macro-F1 are internally consistent.
check(nrow(result$fold_metrics) == 100L,
      all(result$fold_metrics$n == 30L),
      all(result$metrics$macro_f1 > 0.97),
      nrow(result$data_checks) == 7L,
      result$data_checks$value[result$data_checks$check == "duplicate_rows"] == 1)

out <- file.path(run_dir, "iris")
save_analysis(result, out)
manifest <- readRDS(file.path(out, "export_status.rds"))
check(manifest$status == "completed",
      all(file.exists(file.path(out, manifest$files))),
      all(file.info(file.path(out, manifest$files))$size > 0))
restored <- readRDS(file.path(out, "analysis.rds"))
check(identical(result, restored))
check(identical(result$settings, readRDS(file.path(out, "settings.rds"))))
model <- readRDS(file.path(out, "final_model.rds"))
check(identical(predict(model, iris[x])$class,
                predict(result$final_model, iris[x])$class))
for (name in sub("\\.csv$", "", grep("\\.csv$", manifest$files, value = TRUE))) {
  loaded <- read.csv(file.path(out, paste0(name, ".csv")),
                     check.names = FALSE, fileEncoding = "UTF-8")
  expected <- result[[name]]
  check(identical(names(loaded), names(expected)),
        nrow(loaded) == nrow(expected))
  # Empty CSV tables have headers but cannot preserve column storage types.
  if (!nrow(expected)) next
  for (column in names(expected)) {
    if (is.numeric(expected[[column]])) {
      check(isTRUE(all.equal(loaded[[column]], expected[[column]],
                             tolerance = 1e-12, check.attributes = FALSE)))
    } else {
      # CSV blank fields denote absent values; RDS preserves exact types/strings.
      text_value <- function(z) {
        z <- as.character(z)
        z[!is.na(z) & z == ""] <- NA_character_
        z
      }
      check(identical(text_value(loaded[[column]]),
                      text_value(expected[[column]])))
    }
  }
}
text <- readLines(file.path(out, "report.md"), encoding = "UTF-8")
check(identical(text, result$report_lines),
      any(grepl("\u6574\u9ad4\u72c0\u614b\uff1a\u5b8c\u6210", text, fixed = TRUE)))
hashes <- tools::md5sum(file.path(out, manifest$files))
expect_error(save_analysis(result, out))
check(identical(hashes, tools::md5sum(file.path(out, manifest$files))))

# A real CSV failure must leave a failed export status rather than success.
bad_export <- result
as.character.validation_unsupported <- function(x, ...)
  stop("Intentional unsupported CSV column serialization.")
bad_export$cv_predictions$unsupported_column <- structure(
  rep("bad", nrow(bad_export$cv_predictions)), class = "validation_unsupported")
failed_out <- file.path(run_dir, "intentional-export-failure")
expect_error(save_analysis(bad_export, failed_out))
rm(as.character.validation_unsupported)
failed_manifest <- readRDS(file.path(failed_out, "export_status.rds"))
check(failed_manifest$status == "failed",
      nzchar(failed_manifest$error),
      !file.exists(file.path(failed_out, "report.md")))

# Verify Unicode paths, column names, class labels and UTF-8 report round-trips.
unicode_dir <- file.path(run_dir, "\u4e2d\u6587\u8cc7\u6599\u8207\u8f38\u51fa")
dir.create(unicode_dir)
unicode <- iris
names(unicode)[c(1, 5)] <- c("\u82b1\u843c\u9577\u5ea6", "\u54c1\u7a2e")
unicode[["\u54c1\u7a2e"]] <- factor(
  c("\u5c71\u9cf6\u5c3e", "\u8b8a\u8272\u9cf6\u5c3e", "\u7dad\u5409\u5c3c\u4e9e\u9cf6\u5c3e")[as.integer(iris$Species)])
csv <- file.path(unicode_dir, "\u9cf6\u5c3e\u82b1.csv")
write.csv(unicode, csv, row.names = FALSE, fileEncoding = "UTF-8")
u <- run_pipeline(csv, "\u54c1\u7a2e", names(unicode)[1:4])
u_out <- file.path(unicode_dir, "\u5206\u6790\u7d50\u679c")
save_analysis(u, u_out)
check(u$status == "completed",
      identical(readRDS(file.path(u_out, "analysis.rds")), u))
u_diag <- read.csv(file.path(u_out, "diagnostics.csv"),
                   fileEncoding = "UTF-8", check.names = FALSE)
check(setequal(u_diag$group, c(levels(unicode[["\u54c1\u7a2e"]]), "ALL")))
u_report <- readLines(file.path(u_out, "report.md"), encoding = "UTF-8")
check(any(grepl("\u82b1\u843c\u9577\u5ea6", u_report, fixed = TRUE)),
      any(grepl("\u7dad\u5409\u5c3c\u4e9e\u9cf6\u5c3e", u_report, fixed = TRUE)))

# Input and singularity cases.
missing <- iris
missing[1, 1] <- NA_real_
expect_error(run_pipeline(missing, "Species", x))
expect_error(run_pipeline(iris, "unknown", x))
singular <- iris
singular$Petal.Width <- 2 * singular$Petal.Length
blocked <- run_pipeline(singular, "Species", x)
check(blocked$status == "blocked", is.null(blocked$final_model))
set.seed(99)
small <- data.frame(x1 = rnorm(9), x2 = rnorm(9),
                    y = rep(c("A", "B", "C"), each = 3))
partial <- run_pipeline(small, "y", c("x1", "x2"), folds = 3)
check(partial$model_status == "partial",
      partial$metrics$status[partial$metrics$model == "lda"] == "complete",
      partial$cv_leader == "lda", partial$selection == "only_feasible",
      !is.null(partial$final_model),
      partial$metrics$status[partial$metrics$model == "qda"] == "cv_failed",
      is.na(partial$metrics$macro_recall[partial$metrics$model == "qda"]))
# A model that fails in a single fold of a single repeat is not compared
# ("validation incomplete"); the other model is evaluated alone.
one_fail <- local({
  orig <- MASS::qda
  calls <- 0L
  fake <- function(...) {
    calls <<- calls + 1L
    if (calls == 37L) stop("forced failure in one fold") else orig(...)
  }
  ns <- asNamespace("MASS")
  unlockBinding("qda", ns); assign("qda", fake, envir = ns)
  on.exit({ assign("qda", orig, envir = ns); lockBinding("qda", ns) })
  run_analysis(iris, "Species", x)
})
check(one_fail$metrics$status[one_fail$metrics$model == "qda"] == "cv_failed",
      one_fail$status == "partial", one_fail$cv_leader == "lda",
      one_fail$selection == "only_feasible", inherits(one_fail$final_model, "lda"),
      any(grepl("\u53ea\u6709 LDA", selection_text(c(one_fail,
        list(class_counts = data.frame(n = 150)))))),
      nrow(subset(one_fail$issues, model == "qda" & type == "error")) == 1L,
      nrow(subset(one_fail$repeat_metrics, model == "qda")) == 0L)
# Deliberately emulate an unavailable diagnostic in the report only.
missing_diag <- result
missing_diag$diagnostics$status <- "missing_dependency"
missing_diag$diagnostics$statistic <- NA_real_
missing_diag$diagnostics$p_value <- NA_real_
missing_diag$diagnostics$note <- "dependency unavailable (report test)"
missing_diag$diagnostic_status <- "partial"
missing_diag$status <- "partial"
missing_text <- make_report(missing_diag)
check(any(grepl("\u672a\u57f7\u884c\uff1a\u7f3a\u5c11\u5957\u4ef6", missing_text, fixed = TRUE)),
      any(grepl("\u6574\u9ad4\u72c0\u614b\uff1a\u90e8\u5206\u5b8c\u6210", missing_text, fixed = TRUE)))

# Actual EDA images for visual review, produced separately from save_analysis().
fig_dir <- file.path(out, "figures")
dir.create(fig_dir)
draw <- function() {
  grDevices::png(file.path(fig_dir, "eda-%02d.png"),
                 width = 1600, height = 1000, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_eda(iris, "Species", x)
}
draw()
saved <- save_figures(result, iris, file.path(out, "figures-saved"))
check(length(saved) == 6L,
      all(file.info(file.path(out, "figures-saved", saved))$size > 0))
figures <- list.files(fig_dir, pattern = "\\.png$", full.names = TRUE)
check(length(figures) > 0L, all(file.info(figures)$size > 0))

# Command-line entry point on the bundled demo data.
cli_out <- file.path(run_dir, "cli-iris")
rscript <- file.path(R.home("bin"), "Rscript")
status <- system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")),
  shQuote(file.path(project, "assets", "iris.csv")),
  "--target", "Species", "--out", shQuote(cli_out)),
  stdout = FALSE, stderr = FALSE)
check(status == 0L,
      all(file.exists(file.path(cli_out, c("report.md", "report.txt",
        "summary.json", "metrics.csv", "variable_info.csv", "report.html",
        "report.qmd", "input_data.csv", "run_analysis.R")))),
      length(list.files(file.path(cli_out, "figures"))) == 6L)
cli_miss <- file.path(run_dir, "cli-iris-missing")
status_m <- system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")),
  shQuote(file.path(project, "assets", "iris_missing.csv")),
  "--target", "Species", "--impute", "median", "--out", shQuote(cli_miss)),
  stdout = FALSE, stderr = FALSE)
check(status_m == 0L, file.exists(file.path(cli_miss, "imputation.csv")),
      any(grepl("impute = \"median\"", readLines(file.path(cli_miss, "report.qmd"),
                                                  encoding = "UTF-8"), fixed = TRUE)))
# Fixed 5-fold: a class with < 5 rows stops the run; --folds is rejected;
# repeated IDs / dates warn; constant columns are excluded from X.
small_csv <- file.path(run_dir, "small.csv")
write.csv(iris[c(1:4, 51:150), ], small_csv, row.names = FALSE)
small_out <- file.path(run_dir, "cli-small")
st_small <- suppressWarnings(system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")), shQuote(small_csv),
  "--target", "Species", "--out", shQuote(small_out)), stdout = FALSE, stderr = FALSE))
st_folds <- suppressWarnings(system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")),
  shQuote(file.path(project, "assets", "iris.csv")), "--target", "Species",
  "--folds", "3", "--out", shQuote(file.path(run_dir, "cli-folds"))),
  stdout = FALSE, stderr = FALSE))
check(st_small != 0L, !file.exists(small_out), st_folds != 0L)
set.seed(3)
rep_csv <- file.path(run_dir, "repeated.csv")
write.csv(data.frame(patient_id = rep(1:20, each = 3), const = 1,
                     x1 = rnorm(60), x2 = rnorm(60), y = rep(c("a", "b"), 30)),
          rep_csv, row.names = FALSE)
rep_out <- file.path(run_dir, "cli-repeated")
st_rep <- system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")), shQuote(rep_csv),
  "--target", "y", "--out", shQuote(rep_out)), stdout = FALSE, stderr = FALSE)
rep_log <- readLines(file.path(rep_out, "report.txt"), encoding = "UTF-8")
rep_vars <- read.csv(file.path(rep_out, "variable_info.csv"), encoding = "UTF-8")
check(st_rep == 0L, any(grepl("possible repeated subjects", rep_log, fixed = TRUE)),
      !("const" %in% rep_vars$column[grepl("^X", rep_vars$role)]),
      grepl("\u5e38\u6578", rep_vars$role[rep_vars$column == "const"]))
# A column that only LOOKS like an ID (all-distinct whole numbers) stays in X
# with a warning; a column NAMED like an ID is excluded. --plan analyses nothing.
set.seed(4)
look_csv <- file.path(run_dir, "lookalike.csv")
look <- iris
look$leaf_count <- sample(1000:5000, 150)
look$plant_id <- rep(1:50, each = 3)
write.csv(look, look_csv, row.names = FALSE)
plan_out <- suppressWarnings(system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")), shQuote(look_csv),
  "--plan"), stdout = TRUE, stderr = TRUE))
look_out <- file.path(run_dir, "cli-lookalike")
st_look <- system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")), shQuote(look_csv),
  "--target", "Species", "--out", shQuote(look_out)), stdout = FALSE, stderr = FALSE)
look_vars <- read.csv(file.path(look_out, "variable_info.csv"), encoding = "UTF-8")
check(any(grepl("Y candidates: Species", plan_out, fixed = TRUE)),
      any(grepl("CONFIRM: leaf_count", plan_out, fixed = TRUE)),
      any(grepl("CONFIRM: plant_id", plan_out, fixed = TRUE)),
      !file.exists(file.path(dirname(look_csv), "lda_qda_lookalike")),
      st_look == 0L,
      grepl("^X", look_vars$role[look_vars$column == "leaf_count"]),
      grepl("^\u6392\u9664", look_vars$role[look_vars$column == "plant_id"]),
      any(grepl("look like an ID were kept in X: leaf_count",
                readLines(file.path(look_out, "report.txt"), encoding = "UTF-8"), fixed = TRUE)),
      grepl('"suspect_id_columns": ["leaf_count"]',
            paste(readLines(file.path(look_out, "summary.json")), collapse = ""), fixed = TRUE))
# Blank Y stops with the row numbers (never imputed, never treated as a class).
ymiss <- iris
ymiss$Species <- as.character(ymiss$Species)
ymiss$Species[c(3, 77)] <- ""
ymiss_csv <- file.path(run_dir, "ymiss.csv")
write.csv(ymiss, ymiss_csv, row.names = FALSE)
ymiss_out <- suppressWarnings(system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")), shQuote(ymiss_csv),
  "--target", "Species", "--out", shQuote(file.path(run_dir, "cli-ymiss"))),
  stdout = TRUE, stderr = TRUE))
check(!is.null(attr(ymiss_out, "status")), !file.exists(file.path(run_dir, "cli-ymiss")),
      any(grepl("2 row(s) (3, 77)", ymiss_out, fixed = TRUE)),
      inherits(tryCatch(run_analysis(ymiss, "Species", x), error = identity), "error"))
cli_metrics <- read.csv(file.path(cli_out, "metrics.csv"))
check(isTRUE(all.equal(cli_metrics$accuracy, result$metrics$accuracy)),
      file.exists(file.path(cli_out, "repeat_metrics.csv")))
cli_json <- paste(readLines(file.path(cli_out, "summary.json")), collapse = "\n")
check(grepl('"repeats": 10', cli_json, fixed = TRUE),
      grepl('"wins": {"lda": 6, "qda": 0, "tie": 4}', cli_json, fixed = TRUE))
st_rep_opt <- suppressWarnings(system2(rscript, c(
  shQuote(file.path(project, "scripts", "lda_qda_pipeline.R")),
  shQuote(file.path(project, "assets", "iris.csv")), "--target", "Species",
  "--repeats", "3", "--out", shQuote(file.path(run_dir, "cli-repeats"))),
  stdout = FALSE, stderr = FALSE))
check(st_rep_opt != 0L, !file.exists(file.path(run_dir, "cli-repeats")))
html <- paste(readLines(file.path(cli_out, "report.html"), encoding = "UTF-8"),
              collapse = "\n")
check(grepl("data:image/png;base64,", html, fixed = TRUE),
      lengths(regmatches(html, gregexpr("<figure>", html, fixed = TRUE))) == 6L,
      grepl("\u9078\u64c7\u898f\u5247", html, fixed = TRUE))
qmd <- readLines(file.path(cli_out, "report.qmd"), encoding = "UTF-8")
check(any(grepl("^```\\{r\\}", qmd)),
      any(grepl("run_pipeline(data, y = \"Species\"", qmd, fixed = TRUE)),
      any(grepl("output-file: report-quarto.html", qmd, fixed = TRUE)))

# Model selection: higher mean macro recall over the repeats wins; tie -> LDA.
check(result$cv_leader == "lda", result$selection == "lda_better",
      isTRUE(all.equal(result$selection_gap,
        mean(subset(result$repeat_metrics, model == "lda")$macro_recall) -
        mean(subset(result$repeat_metrics, model == "qda")$macro_recall))),
      identical(unname(result$selection_wins), c(6L, 0L, 4L)),
      is.null(result$selection_se))
check(any(grepl("\u5168\u90e8 150 \u7b46", selection_text(result), fixed = TRUE)),
      any(grepl("LDA \u8f03\u9ad8 6 \u6b21", selection_text(result), fixed = TRUE)))
# Tie: identical per-repeat scores keep the simpler LDA.
tie <- result
tie$metrics$macro_recall[2] <- tie$metrics$macro_recall[1]
check(grepl("\u9078\u8f03\u7c21\u55ae\u7684 LDA",
            selection_text(within(tie, selection <- "tie_prefer_lda"))[1]))
# Joint reading of Mardia and Box's M: iris is "Mardia not significant,
# Box's M significant"; a skewed class triggers the automatic Mardia note.
dn <- diagnostic_notes(result$diagnostics)
check(nrow(dn$table) == 4L, which(nzchar(dn$table[[4]])) == 2L,
      !any(grepl("\u6ce8\u610f", dn$notes, fixed = TRUE)))
set.seed(11)
skew <- data.frame(x1 = c(rexp(60), rnorm(60)), x2 = c(rexp(60), rnorm(60)),
                   y = rep(c("skewed", "normal"), each = 60))
skew_case <- run_pipeline(skew, "y", c("x1", "x2"))
dn_skew <- diagnostic_notes(skew_case$diagnostics)
check(any(grepl("\u7f85\u5409\u65af\u8ff4\u6b78", dn_skew$notes, fixed = TRUE)),
      any(grepl("\u6ce8\u610f\uff1askewed", dn_skew$notes, fixed = TRUE)),
      !any(grepl("normal", dn_skew$notes[1], fixed = TRUE)),
      which(nzchar(dn_skew$table[[4]])) %in% 3:4,
      any(grepl("\u6ce8\u610f\uff1askewed", make_report(skew_case), fixed = TRUE)))
# Same means, very different spreads: only QDA can separate the classes.
set.seed(7)
spread <- data.frame(x1 = c(rnorm(100, 0, 1), rnorm(100, 0, 4)),
                     x2 = c(rnorm(100, 0, 1), rnorm(100, 0, 4)),
                     y = rep(c("narrow", "wide"), each = 100))
qda_case <- run_pipeline(spread, "y", c("x1", "x2"))
gap <- diff(qda_case$metrics$macro_recall)
check(gap > 0, qda_case$cv_leader == "qda",
      qda_case$selection == "qda_better")

# Missing values: stop without consent; with impute = "median" the medians
# come from the training rows of each fold only (no leakage).
miss <- read.csv(file.path(project, "assets", "iris_missing.csv"))
check(sum(is.na(miss[x])) == 30L,
      inherits(try(run_analysis(miss, "Species", x), silent = TRUE), "try-error"))
imp <- run_analysis(miss, "Species", x, impute = "median")
check(imp$status == "completed", all(imp$metrics$coverage == 1),
      nrow(imp$imputation) == 200L,
      imp$selection == "qda_better", imp$cv_leader == "qda",
      imp$metrics$macro_recall[imp$metrics$model == "qda"] >
        imp$metrics$macro_recall[imp$metrics$model == "lda"])
fa1 <- subset(imp$fold_assignments, repeat_id == 1)
f1_train <- fa1$fold != 1
expected <- apply(miss[f1_train, x], 2, median, na.rm = TRUE)
in_f1 <- imp$imputation$repeat_id == 1 & imp$imputation$fold == 1
check(isTRUE(all.equal(unname(expected), imp$imputation$training_median[in_f1])))
fa3 <- subset(imp$fold_assignments, repeat_id == 3)
expected3 <- apply(miss[fa3$fold != 2, x], 2, median, na.rm = TRUE)
in_r3f2 <- imp$imputation$repeat_id == 3 & imp$imputation$fold == 2
check(isTRUE(all.equal(unname(expected3), imp$imputation$training_median[in_r3f2])))
leak <- miss
v1 <- which(fa1$fold == 1)
leak[v1, x] <- leak[v1, x] + 100   # wildly change validation rows of fold 1
leak_imp <- run_analysis(leak, "Species", x, impute = "median")
check(identical(leak_imp$imputation$training_median[in_f1],
                imp$imputation$training_median[in_f1]))
piped <- run_pipeline(miss, "Species", x, impute = "median")
check(piped$data_checks$value[piped$data_checks$check == "missing_cells"] == 30,
      sum(piped$class_counts$n) == 150L,
      piped$diagnostic_status == "completed")
summary <- c(
  "# LDA\uff0fQDA \u539f\u578b\u9a57\u6536", "",
  paste0("- \u901a\u904e\u6aa2\u67e5\u6578\uff1a", count),
  "- \u8cc7\u6599\uff1aR \u5167\u5efa iris\uff0c150 \u7b46\u30013 \u985e\u5225\u30014 \u500b\u9810\u6e2c\u8b8a\u6578\u3002",
  "- Mardia \u8207 Box's M \u7531 base R \u8a08\u7b97\uff0c\u8207 MVN 6.3\uff0fbiotools 4.3 \u7684\u5be6\u6e2c\u503c\u4e00\u81f4\u3002",
  "- \u5beb\u51fa\u4e26\u8b80\u56de\uff1aCSV\u3001RDS\u3001\u6700\u7d42\u6a21\u578b\u3001UTF-8 \u7e41\u9ad4\u4e2d\u6587\u5831\u544a\u3002",
  "- \u9a57\u8b49\uff1a\u76f8\u540c\u8a2d\u5b9a\u53ef\u91cd\u73fe\u3001\u6a5f\u7387\u548c\u70ba\u4e00\u3001\u7d50\u679c\u8207\u5957\u4ef6\u76f4\u63a5\u547c\u53eb\u4e00\u81f4\u3002",
  "- \u9a57\u8b49\uff1a\u4e2d\u6587\u8def\u5f91\u3001\u6b04\u540d\u8207\u985e\u5225\u540d\u7a31\u53ef\u8b80\u56de\uff1b\u62d2\u7d55\u8986\u84cb\u65e2\u6709\u8f38\u51fa\u3002",
  "- CSV \u7a7a\u767d\u6b04\u4f4d\u6309\u672a\u63d0\u4f9b\u8655\u7406\uff1b\u5b8c\u6574\u5132\u5b58\u578b\u5225\u8207\u7a7a\u5b57\u4e32\u5dee\u7570\u4ee5 RDS \u4fdd\u5b58\u3002",
  "- \u6545\u610f\u5beb\u6a94\u5931\u6557\u6848\u4f8b\u4fdd\u7559\u65bc intentional-export-failure\uff0c\u5c6c\u9810\u671f\u6e2c\u8a66\u3002",
  "- EDA PNG \u5df2\u8f38\u51fa\uff0c\u7e6a\u5716\u54c1\u8cea\u9700\u53e6\u5916\u76ee\u8996\u6aa2\u67e5\u3002",
  "",
  paste0("R\uff1a", R.version.string, "\uff1bMASS\uff1a",
         as.character(packageVersion("MASS")),
         "\uff1bMVN \u76f4\u63a5\u6bd4\u5c0d\uff1a", has_mvn,
         "\uff1bbiotools \u76f4\u63a5\u6bd4\u5c0d\uff1a", has_biotools),
  "", "## Iris \u4ea4\u53c9\u9a57\u8b49\u7d50\u679c", "",
  markdown_table(result$metrics), "", "## \u8a3a\u65b7\u7d50\u679c", "",
  markdown_table(result$diagnostics), "",
  "## \u9a57\u6536\u7bc4\u570d", "",
  "\u672c\u9a57\u6536\u6db5\u84cb R \u51fd\u5f0f\u8207\u547d\u4ee4\u5217\u5165\u53e3\uff1b\u4e0d\u6db5\u84cb AI \u4f9d SKILL.md \u64cd\u4f5c\u7684\u5c0d\u8a71\u6d41\u7a0b\u3002",
  "\u5c1a\u672a\u5be6\u4f5c\u7fa4\u7d44\uff0f\u6642\u9593\u5207\u5206\u3001\u81ea\u8a02\u5148\u9a57\u3002"
)
writeLines(enc2utf8(summary), file.path(run_dir, "validation-summary.md"),
           useBytes = TRUE)
cat("ALL CHECKS PASSED:", count, "\n")
cat("VALIDATION_DIR:", normalizePath(run_dir, winslash = "/"), "\n")
cat("EDA_IMAGES:", length(figures), "\n")

