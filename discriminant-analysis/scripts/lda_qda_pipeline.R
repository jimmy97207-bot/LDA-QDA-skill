# LDA/QDA command-line pipeline.
#
# Usage (from the skill folder):
#   Rscript scripts/lda_qda_pipeline.R data.csv --target <Y>
#           [--x a,b,c] [--drop id,other] [--dict data_dictionary.csv]
#           [--seed 2026] [--impute none|median] [--dedup] [--out lda_qda_<name>]
#   Rscript scripts/lda_qda_pipeline.R data.csv --plan [--target <Y>]
#           Prints candidate Y columns and the proposed X list (with columns
#           to confirm, and rows identical in Y and X), then exits without
#           analysing. Used for checkpoints A/B/C.
#   --dedup   Remove rows identical to an earlier row in Y and X. Only when the
#             user has confirmed they are double entries; recorded in the report.
#
# Validation is fixed: stratified 5-fold CV repeated 10 times (new random
# folds each repeat). Every class needs >= 5 rows, otherwise the run stops and
# explains why (classes are never merged/dropped).
#
# Steps: read CSV -> variable introduction -> choose X -> data checks ->
# diagnostics (Mardia, Box's M) -> LDA/QDA on the same repeated folds ->
# report.html, report.qmd, report.md, CSV tables, RDS models, figures/*.png,
# report.txt, summary.json.
# Nothing in the data is deleted or transformed; X is imputed only with
# --impute median (fold-wise training medians), and duplicate rows are removed
# only with --dedup.

if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]]))
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))

script_dir <- local({
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
  else "scripts"
})
source(file.path(script_dir, "run_analysis.R"), encoding = "UTF-8")

usage <- paste(
  "Usage: Rscript scripts/lda_qda_pipeline.R data.csv --target Y",
  "[--x a,b,c] [--drop id] [--dict dict.csv] [--seed 2026]",
  "[--impute none|median] [--dedup] [--out dir]",
  "Rscript scripts/lda_qda_pipeline.R data.csv --plan [--target Y]",
  sep = "\n       ")

parse_args <- function(args) {
  opts <- list(seed = "2026", impute = "none")
  positional <- character()
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) { cat(usage, "\n"); quit(status = 0) }
    if (a == "--plan") { opts$plan <- TRUE; i <- i + 1L; next }
    if (a == "--dedup") { opts$dedup <- TRUE; i <- i + 1L; next }
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (grepl("=", key)) {
        opts[[sub("=.*", "", key)]] <- sub("^[^=]*=", "", key)
      } else {
        if (i == length(args)) stop("Missing value for ", a)
        opts[[key]] <- args[[i + 1L]]
        i <- i + 1L
      }
    } else positional <- c(positional, a)
    i <- i + 1L
  }
  if (length(positional) != 1L) stop("Give exactly one CSV file.\n", usage)
  if (is.null(opts$target) && !isTRUE(opts$plan))
    stop("--target is required.\n", usage)
  if (!is.null(opts$folds))
    stop("This version always uses stratified 5-fold CV; --folds is not available.")
  if (!is.null(opts$repeats))
    stop("This version always repeats the 5-fold CV 10 times; --repeats is not available.")
  known <- c("target", "x", "drop", "dict", "seed", "out", "impute", "plan", "dedup")
  bad <- setdiff(names(opts), known)
  if (length(bad)) stop("Unknown option(s): ", paste(bad, collapse = ", "))
  opts$csv <- positional
  opts
}
split_list <- function(s) if (is.null(s)) character() else
  trimws(strsplit(s, ",", fixed = TRUE)[[1L]])

opts <- parse_args(commandArgs(trailingOnly = TRUE))
csv <- opts$csv
if (!file.exists(csv)) stop("File not found: ", csv)
stem <- tools::file_path_sans_ext(basename(csv))
out <- if (is.null(opts$out)) paste0("lda_qda_", stem) else opts$out
plan_mode <- isTRUE(opts$plan)
if (!plan_mode && file.exists(out))
  stop("Output folder already exists: ", out,
       "\nChoose a new --out so earlier results are kept.")

log_lines <- character()
say <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

data <- read.csv(csv, check.names = FALSE, fileEncoding = "UTF-8-BOM")
# Columns whose NAME marks them as an ID are excluded automatically; columns
# that only LOOK like an ID (all distinct whole numbers) stay in X and are
# flagged so the AI asks the user (checkpoint B) instead of dropping them.
id_name <- function(name) grepl("(^id$|_id$|^id_|編號)", name, ignore.case = TRUE)
looks_like_id <- function(v) is.numeric(v) && length(v) > 20 &&
  all(v == round(v), na.rm = TRUE) && !anyDuplicated(v[!is.na(v)])
y_candidates <- Filter(function(nm) {
  v <- data[[nm]]; k <- length(unique(v[!is.na(v)]))
  !id_name(nm) && k >= 2 && ((!is.numeric(v) && k <= 20) ||
    (is.numeric(v) && all(v == round(v), na.rm = TRUE) && k <= 10))
}, names(data))
y <- opts$target
if (is.null(y)) y <- if (length(y_candidates)) y_candidates[[1]] else names(data)[1]
if (!y %in% names(data)) stop("Target column not found: ", y)
# Y is the class label: missing or blank values are never imputed, guessed or
# dropped (checkpoint D). Report the rows and stop before anything else.
y_raw <- data[[y]]
bad_y <- which(is.na(y_raw) | !nzchar(trimws(as.character(y_raw))))
if (length(bad_y) && !isTRUE(opts$plan))
  stop("Missing or blank values in Y (", y, "): ", length(bad_y), " row(s) (",
       paste(utils::head(bad_y, 20), collapse = ", "),
       if (length(bad_y) > 20) ", ..." else "",
       "). Y is the class label and is never imputed, guessed or dropped; ",
       "stop and ask the user (SKILL.md checkpoint D).")

# ---- dictionary (optional): columns 欄位,說明,單位,角色 ----
dict_path <- opts$dict
if (is.null(dict_path)) {
  guess <- file.path(dirname(csv), paste0(stem, "_dictionary.csv"))
  if (file.exists(guess)) dict_path <- guess
}
dict <- NULL
dataset_description <- NULL
if (!is.null(dict_path)) {
  dict <- read.csv(dict_path, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                   stringsAsFactors = FALSE)
  names(dict)[1:min(4, ncol(dict))] <-
    c("column", "description", "unit", "role")[1:min(4, ncol(dict))]
  ds <- dict$column == "__dataset__"
  if (any(ds)) dataset_description <- dict$description[ds][1]
  dict <- dict[!ds, , drop = FALSE]
}

# ---- choose X ----
drop <- split_list(opts$drop)
is_constant <- function(v) is.numeric(v) &&
  length(unique(v[!is.na(v)])) <= 1L
roles <- setNames(rep("", ncol(data)), names(data))
roles[y] <- "Y（分類變數）"
if (length(opts$x)) {
  x <- split_list(opts$x)
  if (!all(x %in% names(data)))
    stop("Unknown --x column(s): ", paste(setdiff(x, names(data)), collapse = ", "))
  roles[x] <- "X（預測變數）"
  rest <- setdiff(names(data), c(y, x))
  roles[rest] <- "未使用（未列於 --x）"
} else {
  x <- character()
  for (nm in setdiff(names(data), y)) {
    v <- data[[nm]]
    if (nm %in% drop) roles[nm] <- "排除（--drop）"
    else if (id_name(nm)) roles[nm] <- "排除（名稱為編號）"
    else if (!is.numeric(v)) roles[nm] <- "排除（非數值，本版不支援）"
    else if (is_constant(v)) roles[nm] <- "排除（常數欄位）"
    else if (looks_like_id(v)) {
      roles[nm] <- "X（預測變數；疑似編號，請確認）"; x <- c(x, nm)
    } else { roles[nm] <- "X（預測變數）"; x <- c(x, nm) }
  }
}
suspect_ids <- names(roles)[grepl("疑似編號", roles)]

# Rows identical to an earlier row in Y and every X. Kept by default; above
# 5% of the rows the AI pauses and asks (checkpoint C) before analysing.
dup_share_limit <- 0.05
dup_rows <- if (length(x)) which(duplicated(data[c(y, x)])) else integer()
dup_text <- if (length(dup_rows)) {
  key <- do.call(paste, c(lapply(data[c(y, x)], as.character), sep = "\r"))
  first <- match(key[dup_rows], key)
  pairs <- paste0(dup_rows, "=", first)
  paste0(length(dup_rows), " row(s) identical in Y and X to an earlier row (row=earlier row: ",
         paste(utils::head(pairs, 20), collapse = ", "), if (length(pairs) > 20) ", ..." else "",
         "; rows counted without the header), ",
         sprintf("%.1f%%", 100 * length(dup_rows) / nrow(data)), " of the data")
} else ""
dup_many <- length(dup_rows) > dup_share_limit * nrow(data)

if (plan_mode) {
  cat("== Plan (nothing analysed) ==\n")
  cat("Data:", csv, "(", nrow(data), "rows,", ncol(data), "columns )\n")
  cat("Y candidates:", if (length(y_candidates)) paste(vapply(y_candidates, function(nm) {
    tb <- table(data[[nm]])
    paste0(nm, " (", length(tb), " classes: ",
           paste(names(tb), tb, sep = "=", collapse = ", "), ")")
  }, ""), collapse = " | ") else "none found; ask the user", "\n")
  if (is.null(opts$target))
    cat("X plan below assumes Y =", y, "; rerun with --plan --target <Y> if the user picks another.\n")
  else cat("Y:", y, "\n")
  cat("X (", length(x), "): ", paste(x, collapse = ", "), "\n", sep = "")
  ex <- roles[grepl("^排除|^未使用", roles)]
  if (length(ex)) cat("Excluded: ", paste0(names(ex), " [", ex, "]", collapse = "; "), "\n", sep = "")
  for (nm in suspect_ids)
    cat("CONFIRM: ", nm, " has all-distinct whole numbers and looks like an ID; ",
        "kept in X until the user confirms (checkpoint B).\n", sep = "")
  rep_ids <- Filter(function(nm) id_name(nm) &&
    anyDuplicated(data[[nm]][!is.na(data[[nm]])]) > 0, names(data))
  if (length(rep_ids))
    cat("CONFIRM: ", paste(rep_ids, collapse = ", "),
        " has repeated values; one subject may have several rows (checkpoint C).\n", sep = "")
  if (length(dup_rows))
    cat(if (dup_many) "CONFIRM: " else "NOTE: ", dup_text,
        if (dup_many) ". More than 5%: pause and ask whether these are double entries (checkpoint C) before analysing.\n"
        else ". Kept by default; mention it in question 2 and use --dedup only if the user says they are double entries.\n",
        sep = "")
  if (length(bad_y))
    cat("STOP: Y (", y, ") has ", length(bad_y), " missing or blank row(s): ",
        paste(utils::head(bad_y, 20), collapse = ", "),
        "; Y is never imputed (checkpoint D).\n", sep = "")
  quit(status = 0)
}
if (!length(x)) stop("No numeric predictor columns were selected.")

variable_info <- data.frame(
  column = names(data),
  description = "", unit = "",
  role = unname(roles),
  type = vapply(data, function(v) class(v)[1], character(1)),
  missing = vapply(data, function(v) sum(is.na(v)), integer(1)),
  stringsAsFactors = FALSE, row.names = NULL)
if (!is.null(dict)) {
  m <- match(variable_info$column, dict$column)
  has <- !is.na(m)
  if ("description" %in% names(dict))
    variable_info$description[has] <- dict$description[m[has]]
  if ("unit" %in% names(dict))
    variable_info$unit[has] <- dict$unit[m[has]]
}

say("== LDA/QDA pipeline ==")
say("Data: ", csv, " (", nrow(data), " rows, ", ncol(data), " columns)")
if (!is.null(dataset_description)) say("Description: ", dataset_description)
if (!is.null(dict_path)) say("Dictionary: ", dict_path)
say("Y: ", y, " | classes: ",
    paste(names(table(data[[y]])), table(data[[y]]), sep = "=", collapse = ", "))
say("X (", length(x), "): ", paste(x, collapse = ", "))
excluded <- variable_info[grepl("^排除|^未使用", variable_info$role), ]
if (nrow(excluded))
  say("Excluded: ", paste0(excluded$column, " [", excluded$role, "]",
                           collapse = "; "))
if (length(suspect_ids))
  say("WARNING column(s) that look like an ID were kept in X: ",
      paste(suspect_ids, collapse = ", "),
      ". Confirm they are real measurements (SKILL.md checkpoint B); otherwise rerun with --drop.")

folds <- 5L
class_n <- table(droplevels(factor(data[[y]])))
if (any(class_n < folds))
  stop("Class(es) with fewer than 5 rows: ",
       paste0(names(class_n)[class_n < folds], "=", class_n[class_n < folds], collapse = ", "),
       ". Stratified 5-fold CV needs >= 5 rows per class. ",
       "This version stops here; it does not merge or drop classes.")

# Signs that one subject may have several rows (checkpoint C): an ID-like
# column with repeated values, or a date/time column. The script only warns;
# the AI must confirm the unit of analysis with the user before trusting CV.
dup_ids <- names(data)[vapply(names(data), function(nm)
  grepl("(^id$|_id$|^id_|編號)", nm, ignore.case = TRUE) &&
    anyDuplicated(data[[nm]][!is.na(data[[nm]])]) > 0, logical(1))]
time_cols <- names(data)[grepl("(date|time|日期|時間|年月)", names(data), ignore.case = TRUE) |
  vapply(data, function(v) inherits(v, c("Date", "POSIXt")), logical(1))]
if (length(dup_ids))
  say("WARNING possible repeated subjects: ID-like column(s) with repeated values: ",
      paste(dup_ids, collapse = ", "),
      ". Confirm the unit of analysis (SKILL.md checkpoint C) before using these results.")
if (length(time_cols))
  say("WARNING date/time column(s) present: ", paste(time_cols, collapse = ", "),
      ". Confirm rows are independent (SKILL.md checkpoint C).")
removed <- integer()
if (length(dup_rows)) {
  if (isTRUE(opts$dedup)) {
    removed <- dup_rows
    data <- data[-dup_rows, , drop = FALSE]
    rownames(data) <- NULL
    say("Removed (--dedup, user confirmed double entries): ", dup_text, ".")
  } else {
    say(if (dup_many) "WARNING " else "NOTE ", dup_text, "; kept",
        if (dup_many) ". More than 5% of the data: confirm with the user (SKILL.md checkpoint C)." else ".")
  }
}
seed <- as.integer(opts$seed)
if (!opts$impute %in% c("none", "median")) stop("--impute must be none or median.")
result <- run_pipeline(data, y, x, folds = folds, seed = seed,
                       variable_info = variable_info,
                       dataset_description = dataset_description,
                       impute = opts$impute, removed_duplicates = removed)

say("")
say("-- Data checks --")
for (i in seq_len(nrow(result$data_checks))) {
  r <- result$data_checks[i, ]
  say(sprintf("  %-28s %8s  %s", r$check, r$value, r$note))
}
say("")
say("-- Diagnostics --")
for (i in seq_len(nrow(result$diagnostics))) {
  r <- result$diagnostics[i, ]
  say(sprintf("  %-16s %-12s stat=%10.4f  p=%-11s %s", r$test, r$group,
              r$statistic, format.pval(r$p_value, digits = 4, eps = 1e-300), r$status))
}
say("")
say(sprintf("-- Repeated stratified %d-fold CV (%d repeats, seed %d), out-of-fold --",
            folds, result$settings$repeats, seed))
for (i in seq_len(nrow(result$metrics))) {
  r <- result$metrics[i, ]
  say(sprintf("  %-4s %-11s coverage=%.3f  mean macro_recall=%.4f (sd %.4f)  mean correct=%.1f [%s-%s]",
              toupper(r$model), r$status, r$coverage, r$macro_recall,
              r$macro_recall_sd, r$correct_mean, r$correct_min, r$correct_max))
}
rm_tab <- result$repeat_metrics
if (nrow(rm_tab)) for (m in unique(rm_tab$model)) {
  z <- rm_tab[rm_tab$model == m, ]
  say(sprintf("  %-4s correct per repeat: %s", toupper(m),
              paste(z$correct, collapse = " ")))
}
if (!anyNA(result$selection_wins))
  say(sprintf("  Wins per repeat: LDA %d, QDA %d, tie %d",
              result$selection_wins[["lda"]], result$selection_wins[["qda"]],
              result$selection_wins[["tie"]]))
if (identical(opts$impute, "median") && nrow(result$imputation))
  say("Imputation: fold-wise training medians (", result$settings$missing_cells,
      " missing cells); see imputation.csv")
say("Chosen model: ",
    if (is.na(result$cv_leader)) "none" else toupper(result$cv_leader),
    " (", result$selection, "; mean macro recall gap = ",
    sprintf("%.4f", result$selection_gap), "; tie -> LDA)")
for (line in diagnostic_notes(result$diagnostics)$notes) say("Note: ", line)
if (!is.null(result$final_model))
  say("Final model: ", toupper(result$cv_leader), " refitted on all ", nrow(data),
      " rows -> final_model.rds")
say("Status: model=", result$model_status,
    ", diagnostics=", result$diagnostic_status, ", overall=", result$status)
if (nrow(result$issues)) {
  say("Issues:")
  for (i in seq_len(nrow(result$issues)))
    say("  ", paste(result$issues[i, ], collapse = " | "))
}

fig_tmp <- tempfile("figures-")
result$figures <- save_figures(result, data, fig_tmp)
result$report_lines <- make_report(result)
save_analysis(result, out)
dir.create(file.path(out, "figures"))
invisible(file.copy(file.path(fig_tmp, result$figures), file.path(out, "figures")))

# Reports for people: report.html (direct, no Quarto) and report.qmd (live R
# chunks; Render recomputes from the copied data and analysis code).
write_report_html(result, out)
# input_data.csv is the data actually analysed (after --dedup, if used).
if (length(removed)) {
  utils::write.csv(data, file.path(out, "input_data.csv"), row.names = FALSE,
                   fileEncoding = "UTF-8")
} else invisible(file.copy(csv, file.path(out, "input_data.csv")))
invisible(file.copy(file.path(script_dir, "run_analysis.R"),
                    file.path(out, "run_analysis.R")))
write_report_qmd(result, out, title = paste0("LDA／QDA 分析報告：", stem))

json_value <- function(v) {
  if (length(v) != 1L || is.na(v)) return("null")
  if (is.numeric(v)) return(format(signif(v, 10), trim = TRUE))
  paste0('"', gsub('(["\\\\])', "\\\\\\1", as.character(v)), '"')
}
mt <- result$metrics
summary_json <- c(
  "{",
  paste0('  "data": ', json_value(basename(csv)), ","),
  paste0('  "n": ', nrow(data), ","),
  paste0('  "target": ', json_value(y), ","),
  paste0('  "predictors": [', paste(vapply(x, json_value, ""), collapse = ", "), "],"),
  paste0('  "folds": ', folds, ', "repeats": ', result$settings$repeats,
         ', "seed": ', seed, ","),
  paste0('  "status": ', json_value(result$status), ","),
  paste0('  "cv_leader": ', json_value(result$cv_leader), ","),
  paste0('  "selection": ', json_value(result$selection), ","),
  paste0('  "selection_gap": ', json_value(result$selection_gap), ","),
  paste0('  "wins": {"lda": ', json_value(result$selection_wins[["lda"]]),
         ', "qda": ', json_value(result$selection_wins[["qda"]]),
         ', "tie": ', json_value(result$selection_wins[["tie"]]), "},"),
  paste0('  "suspect_id_columns": [', paste(vapply(suspect_ids, json_value, ""), collapse = ", "), "],"),
  paste0('  "duplicate_rows": ', length(dup_rows), ","),
  paste0('  "removed_duplicate_rows": [', paste(removed, collapse = ", "), "],"),
  paste0('  "possible_repeated_ids": [', paste(vapply(dup_ids, json_value, ""), collapse = ", "), "],"),
  paste0('  "impute": ', json_value(opts$impute), ","),
  '  "metrics": {',
  paste0(vapply(seq_len(nrow(mt)), function(i) paste0(
    '    "', mt$model[i], '": {"status": ', json_value(mt$status[i]),
    ', "accuracy": ', json_value(mt$accuracy[i]),
    ', "macro_recall": ', json_value(mt$macro_recall[i]),
    ', "macro_f1": ', json_value(mt$macro_f1[i]),
    ', "correct_mean": ', json_value(mt$correct_mean[i]), "}"), ""),
    c(rep(",", nrow(mt) - 1L), "")),
  "  },",
  paste0('  "box_m_p_value": ',
         json_value(result$diagnostics$p_value[result$diagnostics$test == "Box M"])),
  "}")
writeLines(enc2utf8(summary_json), file.path(out, "summary.json"), useBytes = TRUE)

say("")
say("Output: ", out)
say("  report.html, report.qmd (Render -> report-quarto.html), report.md, summary.json,")
say("  report.txt, *.csv, *.rds, figures/ (", length(result$figures), " PNG)")
writeLines(enc2utf8(log_lines), file.path(out, "report.txt"), useBytes = TRUE)
