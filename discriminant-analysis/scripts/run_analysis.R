# LDA/QDA analysis: core models, diagnostics, figures and reporting.
# Loading this file defines functions only; it does not run analysis or write files.
# Only base R + MASS are required. Mardia and Box's M are computed here directly
# (Mardia, 1970; Box, 1949) and were checked against MVN 6.3 / biotools 4.3.

# ---------------------------------------------------------------------------
# 1. Models and cross-validation
# ---------------------------------------------------------------------------

stratified_folds <- function(Y, folds) {
  # Each class is spread over the folds; leftover observations rotate across
  # classes so the same fold does not receive every class's remainder.
  fold_id <- integer(length(Y))
  offset <- 0L
  for (cl in levels(Y)) {
    i <- which(Y == cl)
    base <- rep(seq_len(folds), times = length(i) %/% folds)
    r <- length(i) %% folds
    extra <- if (r) ((offset + seq_len(r) - 1L) %% folds) + 1L else integer()
    offset <- offset + r
    fold_id[i] <- sample(c(base, extra))
  }
  fold_id
}

class_metrics <- function(truth, predicted, lev) {
  cm <- table(truth = factor(truth, levels = lev),
              predicted = factor(predicted, levels = lev))
  recall <- diag(cm) / rowSums(cm)
  precision <- diag(cm) / colSums(cm)
  f1 <- ifelse(is.finite(precision) & (precision + recall) > 0,
               2 * precision * recall / (precision + recall), 0)
  list(cm = cm, recall = as.numeric(recall),
       precision = as.numeric(precision), f1 = as.numeric(f1),
       accuracy = sum(diag(cm)) / sum(cm))
}

# impute = "none": any missing X stops the analysis.
# impute = "median": missing X values are filled with medians computed from the
# training rows only (inside every fold, then from all rows for the final
# model), so the validation fold never influences the filled values.
run_analysis <- function(data, y, x, folds = 5L, repeats = 10L, seed = 2026L,
                         impute = c("none", "median")) {
  impute <- match.arg(impute)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS is required.")
  if (is.character(data) && length(data) == 1L)
    data <- read.csv(data, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  if (!is.data.frame(data)) stop("data must be a data.frame or CSV path.")
  if (anyDuplicated(names(data))) stop("Duplicate column names.")
  if (!is.character(y) || length(y) != 1L || is.na(y) ||
      !is.character(x) || !length(x) || anyNA(x) ||
      anyDuplicated(x) || y %in% x) stop("Invalid y or x.")
  if (!all(c(y, x) %in% names(data))) stop("Unknown column names.")
  whole_number <- function(z, low) is.numeric(z) && length(z) == 1L &&
    is.finite(z) && z == floor(z) && z >= low && z <= .Machine$integer.max
  if (!whole_number(folds, 2L) || !whole_number(repeats, 1L) ||
      !whole_number(seed, 0L))
    stop("folds must be >= 2, repeats >= 1, seed a nonnegative integer.")
  if (!all(vapply(data[x], is.numeric, logical(1)))) stop("X must be numeric.")
  X <- as.matrix(data[x])
  raw_y <- data[[y]]
  if (!is.atomic(raw_y) || !is.null(dim(raw_y))) stop("Y must be a vector.")
  bad_y <- which(is.na(raw_y) | !nzchar(trimws(as.character(raw_y))) |
                   (is.numeric(raw_y) & !is.finite(suppressWarnings(as.numeric(raw_y)))))
  if (length(bad_y))
    stop("Missing, blank, or non-finite values in Y: ", length(bad_y), " row(s) (",
         paste(utils::head(bad_y, 20), collapse = ", "),
         if (length(bad_y) > 20) ", ..." else "", "). ",
         "Y is the class label and is never imputed, guessed or dropped; ",
         "stop and ask the user (SKILL.md checkpoint D).")
  if (any(is.infinite(X)) || any(is.nan(X))) stop("Non-finite values in X.")
  if (anyNA(X) && impute == "none")
    stop("Missing values in X: ", sum(is.na(X)), " cells. ",
         "Ask the user how to handle them (see SKILL.md, checkpoint D); ",
         "impute = \"median\" fills them inside each training fold.")
  Y <- droplevels(factor(raw_y))
  lev <- levels(Y)
  n <- nrow(X)
  if (length(lev) < 2L) stop("Y needs at least two classes.")
  if (any(table(Y) < folds)) stop("Each class needs at least folds observations.")

  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  # One seed drives all repeats; repeat 1 equals the single split of earlier
  # versions. Every repeat is a new stratified 5-fold partition.
  set.seed(seed)
  fold_ids <- vapply(seq_len(repeats), function(r) stratified_folds(Y, folds),
                     integer(n))
  if (!is.matrix(fold_ids)) fold_ids <- matrix(fold_ids, ncol = repeats)

  issues <- data.frame(model = character(), stage = character(),
                       repeat_id = integer(), fold = integer(),
                       type = character(), message = character())
  attempt <- function(expr, model, stage, r = NA_integer_, fold = NA_integer_) {
    record <- function(type, msg) {
      issues <<- rbind(issues, data.frame(
        model = model, stage = stage, repeat_id = r, fold = fold,
        type = type, message = msg))
    }
    tryCatch(withCallingHandlers(expr, warning = function(w) {
      record("warning", conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      record("error", conditionMessage(e))
      NULL
    })
  }
  train_medians <- function(i) {
    med <- apply(X[i, , drop = FALSE], 2, stats::median, na.rm = TRUE)
    if (any(!is.finite(med)))
      stop("A predictor is entirely missing in the training rows.")
    med
  }
  fill <- function(a, med) {
    if (!anyNA(a)) return(a)
    for (j in seq_len(ncol(a))) a[is.na(a[, j]), j] <- med[[j]]
    a
  }
  fit <- function(model, i) {
    a <- X[i, , drop = FALSE]
    med <- if (impute == "median") train_medians(i) else NULL
    if (!is.null(med)) a <- fill(a, med)
    g <- factor(Y[i], levels = lev)
    if (any(table(g) == 0L)) stop("A training class is missing.")
    centered <- a
    for (cl in lev) {
      j <- which(g == cl)
      centered[j, ] <- sweep(a[j, , drop = FALSE], 2,
                             colMeans(a[j, , drop = FALSE]), "-")
      if (model == "qda" &&
          qr(centered[j, , drop = FALSE])$rank < ncol(a))
        stop(paste("Singular within-class covariance:", cl))
    }
    if (model == "lda" && qr(centered)$rank < ncol(a))
      stop("Singular pooled within-class covariance.")
    fn <- if (model == "lda") MASS::lda else MASS::qda
    out <- fn(a, grouping = g, prior = as.numeric(prop.table(table(g))),
              method = "moment")
    attr(out, "impute_medians") <- med
    out
  }
  imputation_log <- list()
  if (impute == "median" && anyNA(X)) {
    for (r in seq_len(repeats)) for (f in seq_len(folds)) {
      train <- which(fold_ids[, r] != f); test <- which(fold_ids[, r] == f)
      med <- tryCatch(train_medians(train), error = function(e)
        rep(NA_real_, ncol(X)))
      imputation_log[[length(imputation_log) + 1L]] <- data.frame(
        repeat_id = r, fold = f, variable = colnames(X),
        training_median = unname(med),
        filled_in_training = colSums(is.na(X[train, , drop = FALSE])),
        filled_in_validation = colSums(is.na(X[test, , drop = FALSE])),
        row.names = NULL)
    }
  }
  predictions <- probabilities <- metrics <- recalls <- confusions <- list()
  fold_metrics <- repeat_metrics <- list()
  full_models <- list()
  for (model in c("lda", "qda")) {
    full <- attempt(fit(model, seq_len(n)), model, "full_fit")
    full_models[model] <- list(full)
    pooled_truth <- pooled_pred <- character()
    all_complete <- !is.null(full)
    n_predicted <- 0L
    for (r in seq_len(repeats)) {
      predicted <- rep(NA_character_, n)
      if (!is.null(full)) {
        for (f in seq_len(folds)) {
          train <- which(fold_ids[, r] != f)
          test <- which(fold_ids[, r] == f)
          p <- attempt({
            fitted <- fit(model, train)
            newx <- X[test, , drop = FALSE]
            med <- attr(fitted, "impute_medians")
            if (!is.null(med)) newx <- fill(newx, med)
            z <- predict(fitted, newx)
            posterior <- z$posterior[, lev, drop = FALSE]
            if (anyNA(z$class) || any(!as.character(z$class) %in% lev) ||
                any(!is.finite(posterior)) || any(posterior < 0) ||
                any(posterior > 1) ||
                any(abs(rowSums(posterior) - 1) > 1e-8))
              stop("Invalid predictions or posterior probabilities.")
            z
          }, model, "cross_validation", r, f)
          if (!is.null(p)) {
            predicted[test] <- as.character(p$class)
            probabilities[[length(probabilities) + 1L]] <- data.frame(
              model = model, repeat_id = r, row_id = rep(test, times = length(lev)),
              fold = f, class = rep(lev, each = length(test)),
              probability = as.vector(p$posterior[, lev, drop = FALSE]))
            fm <- class_metrics(Y[test], predicted[test], lev)
            fold_metrics[[length(fold_metrics) + 1L]] <- data.frame(
              model = model, repeat_id = r, fold = f, n = length(test),
              accuracy = fm$accuracy, macro_recall = mean(fm$recall, na.rm = TRUE))
          } else {
            fold_metrics[[length(fold_metrics) + 1L]] <- data.frame(
              model = model, repeat_id = r, fold = f, n = length(test),
              accuracy = NA_real_, macro_recall = NA_real_)
          }
        }
      }
      n_predicted <- n_predicted + sum(!is.na(predicted))
      complete_r <- !is.null(full) && all(!is.na(predicted))
      all_complete <- all_complete && complete_r
      if (complete_r) {
        cm_r <- class_metrics(Y, predicted, lev)
        repeat_metrics[[length(repeat_metrics) + 1L]] <- data.frame(
          model = model, repeat_id = r, correct = sum(predicted == as.character(Y)),
          accuracy = cm_r$accuracy, macro_recall = mean(cm_r$recall),
          macro_f1 = mean(cm_r$f1))
        pooled_truth <- c(pooled_truth, as.character(Y))
        pooled_pred <- c(pooled_pred, predicted)
      }
      predictions[[length(predictions) + 1L]] <- data.frame(
        model = model, repeat_id = r, row_id = seq_len(n), fold = fold_ids[, r],
        truth = as.character(Y), predicted = predicted)
    }
    state <- if (is.null(full)) "infeasible" else if (all_complete) "complete" else "cv_failed"
    m_row <- data.frame(model = model, status = state,
      coverage = n_predicted / (n * repeats), accuracy = NA_real_,
      macro_recall = NA_real_, macro_recall_sd = NA_real_, macro_f1 = NA_real_,
      correct_mean = NA_real_, correct_min = NA_real_, correct_max = NA_real_)
    if (all_complete) {
      rm_df <- do.call(rbind, Filter(function(z) z$model[1] == model, repeat_metrics))
      m_row$accuracy <- mean(rm_df$accuracy)
      m_row$macro_recall <- mean(rm_df$macro_recall)
      m_row$macro_recall_sd <- if (repeats > 1L) stats::sd(rm_df$macro_recall) else NA_real_
      m_row$macro_f1 <- mean(rm_df$macro_f1)
      m_row$correct_mean <- mean(rm_df$correct)
      m_row$correct_min <- min(rm_df$correct)
      m_row$correct_max <- max(rm_df$correct)
      # Pooled over all repeats: every row predicted once per repeat.
      pc <- class_metrics(pooled_truth, pooled_pred, lev)
      recalls[[model]] <- data.frame(model = model, class = lev,
        recall = pc$recall, precision = pc$precision, f1 = pc$f1)
      cmd <- as.data.frame(pc$cm)
      totals <- tapply(cmd$Freq, cmd$truth, sum)
      cmd$row_pct <- cmd$Freq / as.numeric(totals[as.character(cmd$truth)])
      confusions[[model]] <- cbind(model = model, cmd)
    }
    # A model that misses any fold of any repeat is "validation incomplete":
    # none of its scores are reported or compared.
    if (!all_complete)
      repeat_metrics <- Filter(function(z) z$model[1] != model, repeat_metrics)
    metrics[[model]] <- m_row
  }
  bind <- function(z) if (length(z)) do.call(rbind, z) else data.frame()
  metrics <- bind(metrics)
  repeat_metrics <- bind(repeat_metrics)
  ok <- metrics$status == "complete"
  status <- if (any(metrics$status == "cv_failed")) "partial" else
    if (any(ok)) "completed" else "blocked"
  # Selection rule: the higher mean macro recall over the repeated stratified
  # 5-fold CV wins; an exact tie keeps the simpler LDA. Win counts per repeat
  # are reported as a stability indicator, not used for the decision.
  # "partial" is a warning only: if one model completed every repeat and the
  # other is "validation incomplete" (cv_failed) or infeasible, the complete
  # model is selected alone.
  leader <- NA_character_
  selection <- "none"
  sel_gap <- NA_real_
  wins <- c(lda = NA_integer_, qda = NA_integer_, tie = NA_integer_)
  if (any(ok)) {
    mr <- setNames(replace(metrics$macro_recall, !ok, NA), metrics$model)
    if (ok[metrics$model == "lda"] && ok[metrics$model == "qda"]) {
      a <- repeat_metrics$macro_recall[repeat_metrics$model == "lda"]
      b <- repeat_metrics$macro_recall[repeat_metrics$model == "qda"]
      d <- b - a
      wins <- c(lda = sum(d < -1e-12), qda = sum(d > 1e-12),
                tie = sum(abs(d) <= 1e-12))
      sel_gap <- mr[["qda"]] - mr[["lda"]]
      if (abs(sel_gap) < 1e-12) {
        leader <- "lda"; selection <- "tie_prefer_lda"; sel_gap <- 0
      } else if (sel_gap > 0) {
        leader <- "qda"; selection <- "qda_better"
      } else {
        leader <- "lda"; selection <- "lda_better"
      }
      sel_gap <- abs(sel_gap)
    } else {
      leader <- metrics$model[ok][1]
      selection <- "only_feasible"
    }
  }
  list(
    status = status, cv_leader = leader, selection = selection,
    selection_gap = sel_gap, selection_wins = wins,
    settings = list(y = y, x = x, folds = folds, repeats = repeats, seed = seed,
      class_levels = lev, prior = "training proportions",
      primary_metric = "macro_recall",
      selection_rule = "higher mean macro recall over repeated stratified 5-fold CV; tie -> LDA",
      impute = impute,
      impute_medians = if (impute == "median") train_medians(seq_len(n)) else NULL,
      missing_cells = sum(is.na(X)),
      R = R.version.string,
      MASS = as.character(utils::packageVersion("MASS"))),
    fold_assignments = data.frame(row_id = rep(seq_len(n), repeats),
      repeat_id = rep(seq_len(repeats), each = n), fold = as.vector(fold_ids)),
    metrics = metrics, repeat_metrics = repeat_metrics,
    class_recall = bind(recalls),
    fold_metrics = bind(fold_metrics),
    cv_predictions = bind(predictions),
    cv_probabilities = bind(probabilities),
    confusion_matrix = bind(confusions), issues = issues,
    imputation = bind(imputation_log),
    candidate_models = full_models,
    final_model = if (!is.na(leader)) full_models[[leader]] else NULL)
}

# ---------------------------------------------------------------------------
# 2. Data checks and assumption diagnostics
# ---------------------------------------------------------------------------

diagnostic_input <- function(data, y, x) {
  if (!is.data.frame(data) || anyDuplicated(names(data)))
    stop("data must be a data.frame with unique column names.")
  if (!is.character(y) || length(y) != 1L || is.na(y) ||
      !is.character(x) || !length(x) || anyNA(x) ||
      anyDuplicated(x) || y %in% x ||
      !all(c(y, x) %in% names(data))) stop("Invalid y or x.")
  if (!all(vapply(data[x], is.numeric, logical(1))))
    stop("X must be numeric.")
  X <- as.matrix(data[x])
  raw <- data[[y]]
  if (!is.atomic(raw) || !is.null(dim(raw)) || anyNA(raw) ||
      any(!nzchar(trimws(as.character(raw)))) ||
      (is.numeric(raw) && any(!is.finite(raw))) ||
      any(!is.finite(X))) stop("Missing, blank, or non-finite values.")
  Y <- droplevels(factor(raw))
  if (nlevels(Y) < 2L) stop("Y needs at least two classes.")
  list(X = X, Y = Y)
}

# Mardia (1970) multivariate skewness and kurtosis, MLE covariance (divisor n).
# For n < 20 the small-sample skewness correction is applied, as in MVN.
mardia_test <- function(z) {
  n <- nrow(z); p <- ncol(z)
  d <- sweep(z, 2, colMeans(z))
  S <- crossprod(d) / n
  D <- d %*% solve(S, t(d))
  b1 <- sum(D^3) / n^2
  b2 <- sum(diag(D)^2) / n
  df <- p * (p + 1) * (p + 2) / 6
  skew <- if (n < 20) {
    k <- (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
    n * k * b1 / 6
  } else n * b1 / 6
  kurt <- (b2 - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  data.frame(
    Test = c("Mardia Skewness", "Mardia Kurtosis"),
    Statistic = c(skew, kurt),
    p.value = c(stats::pchisq(skew, df, lower.tail = FALSE),
                2 * stats::pnorm(-abs(kurt))))
}

# Box (1949) M test with the chi-square approximation.
box_m_test <- function(X, Y) {
  g <- split(seq_len(nrow(X)), Y)
  K <- length(g); p <- ncol(X)
  ni <- vapply(g, length, integer(1))
  Si <- lapply(g, function(i) stats::cov(X[i, , drop = FALSE]))
  Sp <- Reduce(`+`, Map(function(S, m) (m - 1) * S, Si, ni)) / sum(ni - 1)
  logdet <- function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus)
  M <- sum(ni - 1) * logdet(Sp) -
    sum((ni - 1) * vapply(Si, logdet, numeric(1)))
  c1 <- (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (K - 1)) *
    (sum(1 / (ni - 1)) - 1 / sum(ni - 1))
  statistic <- M * (1 - c1)
  df <- p * (p + 1) * (K - 1) / 2
  list(statistic = statistic, df = df,
       p.value = stats::pchisq(statistic, df, lower.tail = FALSE))
}

# Descriptive data-quality checks. Nothing is removed or changed here.
check_data <- function(data, y, x) {
  d <- diagnostic_input(data, y, x)
  X <- d$X; Y <- d$Y
  counts <- table(Y)
  dup <- sum(duplicated(data[c(y, x)]))
  zero_var <- character()
  outlier_rows <- integer()
  for (cl in levels(Y)) {
    ids <- which(Y == cl)
    for (j in seq_len(ncol(X))) {
      z <- X[ids, j]
      if (length(z) > 1L && stats::sd(z) == 0)
        zero_var <- c(zero_var, paste0(cl, ":", colnames(X)[j]))
      q <- stats::quantile(z, c(.25, .75), names = FALSE)
      span <- q[2] - q[1]
      outlier_rows <- union(outlier_rows,
        ids[z < q[1] - 1.5 * span | z > q[2] + 1.5 * span])
    }
  }
  rank_ok <- qr(scale(X, center = TRUE, scale = FALSE))$rank == ncol(X)
  ratio <- max(counts) / min(counts)
  small <- names(counts)[counts <= ncol(X)]
  data.frame(
    check = c("observations", "duplicate_rows", "class_balance_ratio",
              "classes_with_n_le_p", "zero_within_class_variance",
              "full_rank_X", "iqr_outlier_rows"),
    value = c(nrow(X), dup, round(ratio, 3), length(small),
              length(zero_var), as.integer(rank_ok), length(outlier_rows)),
    note = c(
      paste0(nlevels(Y), " classes, ", ncol(X), " predictors"),
      if (dup) "Identical rows on Y and X; kept, not assumed to be the same unit." else "",
      if (ratio > 1.5) "Imbalanced: read macro recall / macro-F1, not accuracy alone." else "",
      if (length(small)) paste("QDA not estimable for:", paste(small, collapse = ", ")) else "",
      paste(zero_var, collapse = ", "),
      if (!rank_ok) "Linear dependence among X; LDA/QDA will be singular." else "",
      "Within-class 1.5 IQR rule per variable; flagged only, all kept."),
    stringsAsFactors = FALSE)
}

diagnose_data <- function(data, y, x) {
  d <- diagnostic_input(data, y, x)
  X <- d$X
  Y <- d$Y
  groups <- split(seq_len(nrow(X)), Y)
  p <- ncol(X)
  summary_rows <- lapply(names(groups), function(cl) {
    z <- X[groups[[cl]], , drop = FALSE]
    do.call(rbind, lapply(seq_len(p), function(j) data.frame(
      class = cl, variable = colnames(z)[j], n = nrow(z),
      mean = mean(z[, j]), sd = stats::sd(z[, j]),
      min = min(z[, j]), median = stats::median(z[, j]),
      max = max(z[, j])
    )))
  })
  tests <- data.frame(
    test = character(), group = character(), n = integer(),
    statistic = numeric(), p_value = numeric(),
    status = character(), note = character())
  add <- function(test, group, n, status, note = "",
                  statistic = NA_real_, p_value = NA_real_) {
    tests <<- rbind(tests, data.frame(
      test = test, group = group, n = n, statistic = statistic,
      p_value = p_value, status = status, note = note))
  }
  capture <- function(expr) {
    warnings <- character()
    answer <- tryCatch(
      withCallingHandlers(expr, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) e)
    list(value = answer, note = paste(warnings, collapse = " | "))
  }
  feasible <- vapply(groups, function(i) {
    z <- X[i, , drop = FALSE]
    nrow(z) > p && qr(scale(z, center = TRUE, scale = FALSE))$rank == p
  }, logical(1))
  raw_results <- list()

  for (cl in names(groups)) {
    z <- X[groups[[cl]], , drop = FALSE]
    if (!feasible[[cl]]) {
      for (nm in c("Mardia Skewness", "Mardia Kurtosis"))
        add(nm, cl, nrow(z), "not_estimable",
            "Insufficient observations or singular covariance.")
      next
    }
    ans <- capture({
      r <- mardia_test(z)
      if (any(!is.finite(r$Statistic)) || any(!is.finite(r$p.value)))
        stop("Invalid Mardia statistics or p-values.")
      r
    })
    if (inherits(ans$value, "error")) {
      for (nm in c("Mardia Skewness", "Mardia Kurtosis"))
        add(nm, cl, nrow(z), "error",
            paste(ans$note, conditionMessage(ans$value)))
    } else {
      r <- ans$value
      raw_results[[cl]] <- r
      for (j in seq_len(nrow(r)))
        add(r$Test[j], cl, nrow(z),
            if (nzchar(ans$note)) "computed_with_warning" else "computed",
            ans$note, r$Statistic[j], r$p.value[j])
    }
  }

  if (!all(feasible)) {
    add("Box M", "ALL", nrow(X), "not_estimable",
        "At least one group has insufficient or singular data.")
  } else {
    ans <- capture({
      r <- box_m_test(X, Y)
      if (!is.finite(r$statistic) || !is.finite(r$p.value))
        stop("Invalid Box M statistic or p-value.")
      r
    })
    if (inherits(ans$value, "error")) {
      add("Box M", "ALL", nrow(X), "error",
          paste(ans$note, conditionMessage(ans$value)))
    } else {
      raw_results[["Box M"]] <- ans$value
      add("Box M", "ALL", nrow(X),
          if (nzchar(ans$note)) "computed_with_warning" else "computed",
          paste0("df = ", ans$value$df,
                 if (nzchar(ans$note)) paste0("; ", ans$note) else ""),
          ans$value$statistic, ans$value$p.value)
    }
  }
  list(
    data_summary = do.call(rbind, summary_rows),
    class_counts = data.frame(class = levels(Y), n = as.integer(table(Y))),
    data_checks = check_data(data, y, x),
    diagnostics = tests, raw = raw_results)
}

# ---------------------------------------------------------------------------
# 3. Figures
# ---------------------------------------------------------------------------

plot_eda <- function(data, y, x) {
  d <- diagnostic_input(data, y, x)
  X <- d$X
  Y <- d$Y
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  for (j in seq_len(ncol(X))) {
    graphics::boxplot(split(X[, j], Y), main = colnames(X)[j],
                      xlab = y, ylab = colnames(X)[j])
  }
  for (cl in levels(Y)) {
    for (j in seq_len(ncol(X))) {
      z <- X[Y == cl, j]
      if (length(z) < 3L || stats::sd(z) == 0) next
      stats::qqnorm(z, main = paste(cl, colnames(X)[j]))
      stats::qqline(z, col = "red")
    }
  }
  if (ncol(X) >= 2L && ncol(X) <= 6L) {
    graphics::par(mfrow = c(1, 1))
    colours <- grDevices::hcl.colors(nlevels(Y), "Dark 3")
    graphics::pairs(
      X, col = colours[as.integer(Y)], pch = 19,
      main = "", oma = c(5, 3, 4, 3))
    graphics::par(fig = c(0, 1, 0, 0.07), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot.new()
    graphics::legend("center", legend = levels(Y), col = colours, pch = 19,
                     ncol = nlevels(Y), bty = "n")
  }
  invisible(NULL)
}

# Writes numbered PNG files; returns the file names that were written.
save_figures <- function(result, data, out) {
  y <- result$settings$y
  x <- result$settings$x
  data <- data[stats::complete.cases(data[x]), , drop = FALSE]
  d <- diagnostic_input(data, y, x)
  X <- d$X; Y <- d$Y
  if (!dir.exists(out) && !dir.create(out, recursive = TRUE))
    stop("Cannot create figure directory.")
  colours <- grDevices::hcl.colors(nlevels(Y), "Dark 3")
  written <- character()
  png_file <- function(name, width = 1600, height = 1000) {
    path <- file.path(out, name)
    grDevices::png(path, width = width, height = height, res = 150)
    written <<- c(written, name)
  }
  p <- ncol(X)
  cols <- min(p, 4L); rows <- ceiling(p / cols)

  png_file("01_class_counts.png", 1200, 800)
  graphics::par(mar = c(4, 4.5, 3, 1), las = 1)
  b <- graphics::barplot(table(Y), col = colours, border = NA,
    main = paste("Class counts:", y), ylab = "n")
  graphics::text(b, table(Y), table(Y), pos = 1, col = "white", font = 2)
  grDevices::dev.off()

  png_file("02_boxplots.png", 400 * cols + 200, 650 * rows + 100)
  graphics::par(mfrow = c(rows, cols), mar = c(4, 4.5, 3, 1), las = 1)
  for (j in seq_len(p))
    graphics::boxplot(split(X[, j], Y), main = colnames(X)[j],
                      col = colours, ylab = colnames(X)[j])
  grDevices::dev.off()

  png_file("03_qq_by_class.png", 380 * p + 200, 380 * nlevels(Y) + 100)
  graphics::par(mfrow = c(nlevels(Y), p), mar = c(3.5, 3.5, 2.5, .8),
                mgp = c(2, .7, 0))
  for (cl in levels(Y)) for (j in seq_len(p)) {
    z <- X[Y == cl, j]
    if (length(z) < 3L || stats::sd(z) == 0) {
      graphics::plot.new(); graphics::title(paste(cl, colnames(X)[j]))
      next
    }
    stats::qqnorm(z, main = paste(cl, colnames(X)[j]), cex.main = .9,
                  pch = 19, cex = .6, col = colours[match(cl, levels(Y))])
    stats::qqline(z, col = "grey30")
  }
  grDevices::dev.off()

  if (p >= 2L && p <= 6L) {
    png_file("04_pairs.png", 1600, 1500)
    graphics::pairs(X, col = colours[as.integer(Y)], pch = 19, cex = .6,
                    oma = c(6, 3, 3, 3))
    graphics::par(fig = c(0, 1, 0, .06), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot.new()
    graphics::legend("center", legend = levels(Y), col = colours, pch = 19,
                     ncol = nlevels(Y), bty = "n")
    grDevices::dev.off()
  }

  lda_full <- result$candidate_models[["lda"]]
  if (!is.null(lda_full)) {
    scores <- predict(lda_full, X)$x
    png_file("05_lda_scores.png", 1500, 1000)
    graphics::par(mar = c(4.5, 4.5, 3, 1), las = 1)
    if (ncol(scores) >= 2L) {
      prop <- lda_full$svd^2 / sum(lda_full$svd^2)
      graphics::plot(scores[, 1], scores[, 2], col = colours[as.integer(Y)],
        pch = 19, cex = .8,
        xlab = sprintf("LD1 (%.1f%%)", 100 * prop[1]),
        ylab = sprintf("LD2 (%.1f%%)", 100 * prop[2]),
        main = "LDA discriminant scores (full data, descriptive)")
      graphics::legend("topright", legend = levels(Y), col = colours,
                       pch = 19, bty = "n")
    } else {
      dens <- lapply(levels(Y), function(cl)
        stats::density(scores[Y == cl, 1]))
      graphics::plot(NULL, xlim = range(sapply(dens, function(z) range(z$x))),
        ylim = c(0, max(sapply(dens, function(z) max(z$y)))),
        xlab = "LD1", ylab = "Density",
        main = "LDA discriminant scores (full data, descriptive)")
      for (k in seq_along(dens))
        graphics::lines(dens[[k]], col = colours[k], lwd = 2)
      graphics::legend("topright", legend = levels(Y), col = colours,
                       lwd = 2, bty = "n")
    }
    grDevices::dev.off()
  }

  cmd <- result$confusion_matrix
  if (is.data.frame(cmd) && nrow(cmd)) {
    models <- unique(cmd$model)
    lev <- result$settings$class_levels
    K <- length(lev)
    png_file("06_confusion_matrix.png", 800 * length(models) + 200, 900)
    graphics::par(mfrow = c(1, length(models)), mar = c(5, 7, 4, 1))
    shades <- grDevices::colorRampPalette(c("#f1f7f8", "#087e8b"))(101)
    for (m in models) {
      z <- cmd[cmd$model == m, ]
      # as.data.frame(table(truth, predicted)) varies truth fastest, so a
      # column-major fill gives rows = truth, columns = predicted.
      # Pooled over all repeats and shown as row percentages: of the rows
      # whose true class is i, the share predicted as j.
      cm <- matrix(z$row_pct, K, K, dimnames = list(lev, lev))
      graphics::plot.new()
      graphics::plot.window(xlim = c(.5, K + .5), ylim = c(.5, K + .5))
      for (i in 1:K) for (j in 1:K) {
        yy <- K + 1 - i
        level <- 1L + round(100 * cm[i, j])
        graphics::rect(j - .49, yy - .49, j + .49, yy + .49,
                       col = shades[level], border = "white")
        graphics::text(j, yy, sprintf("%.1f%%", 100 * cm[i, j]), cex = 1.2,
          col = if (cm[i, j] > .5) "white" else "#173448")
      }
      graphics::axis(1, at = 1:K, labels = lev, tick = FALSE, cex.axis = .85)
      graphics::axis(2, at = K:1, labels = lev, las = 2, tick = FALSE,
                     cex.axis = .85)
      graphics::title(main = paste0(toupper(m), " (", result$settings$repeats,
                                     " repeats pooled, row %)"),
                      xlab = "Predicted class")
      graphics::mtext("True class", side = 2, line = 5.5)
    }
    grDevices::dev.off()
  }
  invisible(written)
}

# ---------------------------------------------------------------------------
# 4. Report
# ---------------------------------------------------------------------------

markdown_table <- function(d) {
  if (!is.data.frame(d) || !ncol(d) || !nrow(d))
    return("無可呈現的紀錄。")
  clean <- function(v) {
    if (is.numeric(v)) {
      v <- ifelse(is.finite(v), formatC(v, digits = 4, format = "g"), "未計算")
    } else {
      v <- as.character(v)
      v[is.na(v)] <- "未提供"
    }
    v <- gsub("[\r\n]+", " ", v)
    v <- gsub("&", "&amp;", v, fixed = TRUE)
    v <- gsub("<", "&lt;", v, fixed = TRUE)
    gsub("|", "&#124;", v, fixed = TRUE)
  }
  row <- function(v) paste0("| ", paste(v, collapse = " | "), " |")
  values <- as.data.frame(lapply(d, clean), check.names = FALSE)
  c(row(clean(names(d))), row(rep("---", ncol(d))),
    apply(values, 1, row))
}

# Plain-language explanation of which model was chosen and why.
selection_text <- function(result) {
  m <- result$metrics
  get <- function(model, col) m[[col]][m$model == model]
  pct <- function(v) sprintf("%.2f%%", 100 * v)
  n <- sum(result$class_counts$n)
  R <- result$settings$repeats
  sel <- if (is.null(result$selection)) "none" else result$selection
  wins <- result$selection_wins
  rule <- paste0("選擇規則：兩種方法用相同的切分，做 ", R,
    " 次分層 5-fold 交叉驗證（每次重新隨機分組），比較 ", R,
    " 次的平均 Macro recall，較高者勝出；完全相同時選較簡單的 LDA。",
    "逐次勝出次數只用來說明結果穩不穩定。Mardia、Box's M 只作診斷參考，不用來選方法。")
  one <- function(model) paste0(toupper(model), " 每次平均答對 ",
    sprintf("%.1f", get(model, "correct_mean")), "／", n, " 筆（範圍 ",
    get(model, "correct_min"), "–", get(model, "correct_max"),
    "，平均 Macro recall ", pct(get(model, "macro_recall")), "）")
  both <- function() paste0(R, " 次重複的結果：", one("lda"), "；", one("qda"),
    "。逐次比較 Macro recall：LDA 較高 ", wins[["lda"]], " 次、QDA 較高 ",
    wins[["qda"]], " 次、相同 ", wins[["tie"]], " 次。")
  caveat <- function(winner) {
    loser <- setdiff(c("lda", "qda"), winner)
    d <- abs(get("lda", "correct_mean") - get("qda", "correct_mean"))
    if (wins[[loser]] >= wins[[winner]] || wins[[loser]] >= ceiling(R * 0.3))
      paste0(toupper(loser), " 也在 ", wins[[loser]], " 次中較高，平均只差 ",
             sprintf("%.1f", d), " 筆，結果不穩定，兩者實際表現接近。")
    else if (wins[[loser]] > 0)
      paste0("多數重複由 ", toupper(winner), " 較高，但也有 ", wins[[loser]],
             " 次反過來；平均差 ", sprintf("%.1f", d), " 筆。")
    else paste0(R, " 次中沒有一次由 ", toupper(loser), " 較高；平均差 ",
                sprintf("%.1f", d), " 筆。")
  }
  verdict <- switch(sel,
    qda_better = paste0(both(), "QDA 平均較高，選 QDA。", caveat("qda")),
    lda_better = paste0(both(), "LDA 平均較高，選 LDA。", caveat("lda")),
    tie_prefer_lda = paste0(both(), "兩者平均相同，選較簡單的 LDA。"),
    only_feasible = paste0("只有 ", toupper(result$cv_leader), " 完成全部 ", R,
      " 次驗證；另一個方法驗證不完整或算不出來（原因見「模型警告與錯誤」），",
      "不參加比較，所以選 ", toupper(result$cv_leader), "。"),
    "本次沒有完整的驗證結果，不選擇方法；請參閱各模型的狀態與問題紀錄。")
  final <- if (!is.null(result$final_model)) paste0("選定後，用全部 ", n,
    " 筆資料重新配適 ", toupper(result$cv_leader), "，作為最終模型（final_model.rds）。") else character()
  c(verdict, final, rule)
}

# Joint reading of Mardia and Box's M. Returns the four-case table (the case
# that applies is marked) and plain-language notes for this data set.
diagnostic_notes <- function(diagnostics, alpha = 0.05) {
  d <- diagnostics
  ok <- d$status %in% c("computed", "computed_with_warning") & is.finite(d$p_value)
  mard <- grepl("^Mardia", d$test)
  box <- d$test == "Box M"
  mardia_groups <- unique(d$group[mard & ok & d$p_value < alpha])
  mardia_done <- any(mard & ok)
  box_done <- any(box & ok)
  box_sig <- box_done && any(d$p_value[box & ok] < alpha)
  cases <- data.frame(
    Mardia = c("不顯著", "不顯著", "顯著", "顯著"),
    `Box's M` = c("不顯著", "顯著", "不顯著", "顯著"),
    解讀 = c(
      "未發現明顯違反假設；仍依交叉驗證結果選方法。",
      "各類散開方式可能不同，QDA 有機會較好。",
      "資料不像鐘形，算出的機率可能不準。",
      "Box's M 可能只是受非鐘形影響，不能直接下結論。"),
    check.names = FALSE)
  current <- if (mardia_done && box_done)
    1L + 2L * (length(mardia_groups) > 0L) + as.integer(box_sig) else NA_integer_
  cases$本資料 <- ifelse(seq_len(4) == current, "← 本資料", "")
  notes <- character()
  if (length(mardia_groups))
    notes <- c(notes, paste0("注意：", paste(mardia_groups, collapse = "、"),
      " 的 Mardia 檢定 p < ", alpha, "，這些類別的資料不太像鐘形（多元常態）。",
      "交叉驗證仍可評估 LDA／QDA 的分類表現，但模型的前提可能不成立，表現可能因此受影響；",
      "算出的機率和 Box's M 也要保守解讀。可另外比較不需要常態假設的方法",
      "（例如羅吉斯迴歸、決策樹），但超出本 skill 範圍，需要另外分析。"))
  if (!mardia_done || !box_done)
    notes <- c(notes, "部分診斷無法計算（見上表狀態），合併解讀只供參考。")
  if (!is.na(current))
    notes <- c(notes, paste0("合併解讀：Mardia ", cases$Mardia[current],
      "、Box's M ", cases$`Box's M`[current], "——", cases$解讀[current]))
  list(table = cases, notes = notes)
}

# Pooled confusion matrix of one model as a wide table of row percentages.
confusion_percent <- function(result, model) {
  z <- result$confusion_matrix
  z <- z[z$model == model, ]
  lev <- result$settings$class_levels
  K <- length(lev)
  cm <- matrix(sprintf("%.1f%%", 100 * z$row_pct), K, K)
  out <- data.frame(lev, cm, stringsAsFactors = FALSE)
  names(out) <- c("真實類別 \\ 預測類別", lev)
  out
}

# Per-repeat results of both models side by side.
repeat_table <- function(result) {
  rm <- result$repeat_metrics
  if (!is.data.frame(rm) || !nrow(rm)) return(data.frame())
  out <- data.frame(repeat_id = sort(unique(rm$repeat_id)))
  for (m in unique(rm$model)) {
    z <- rm[rm$model == m, ]
    out[[paste0(toupper(m), " 答對")]] <- z$correct[match(out$repeat_id, z$repeat_id)]
    out[[paste0(toupper(m), " Macro recall")]] <-
      round(z$macro_recall[match(out$repeat_id, z$repeat_id)], 4)
  }
  if (all(c("LDA Macro recall", "QDA Macro recall") %in% names(out))) {
    d <- out[["QDA Macro recall"]] - out[["LDA Macro recall"]]
    out$較高 <- ifelse(abs(d) <= 1e-12, "相同", ifelse(d > 0, "QDA", "LDA"))
  }
  names(out)[1] <- "第幾次"
  out
}

# Fold-wise medians summarised per variable (the full log is imputation.csv).
imputation_summary <- function(result) {
  im <- result$imputation
  if (!is.data.frame(im) || !nrow(im)) return(data.frame())
  out <- do.call(rbind, lapply(split(im, factor(im$variable, unique(im$variable))),
    function(z) data.frame(變數 = z$variable[1],
      中位數最小 = min(z$training_median), 中位數最大 = max(z$training_median),
      每輪補入練習格數 = round(mean(z$filled_in_training), 1),
      每輪補入考試格數 = round(mean(z$filled_in_validation), 1),
      check.names = FALSE, row.names = NULL)))
  rownames(out) <- NULL
  out
}

# ---- Report tables with Chinese headers, and how-to-read notes -------------

pct_text <- function(v) ifelse(is.na(v), "未計算", sprintf("%.2f%%", 100 * v))

rename_cols <- function(d, labels) {
  hit <- names(d) %in% names(labels)
  names(d)[hit] <- unname(labels[names(d)[hit]])
  d
}

class_counts_display <- function(result) {
  d <- result$class_counts
  data.frame(類別 = d$class, 筆數 = d$n,
             比例 = sprintf("%.1f%%", 100 * d$n / sum(d$n)), check.names = FALSE)
}

variable_info_display <- function(result) {
  d <- result$variable_info
  if (!is.data.frame(d) || !nrow(d)) return(d)
  rename_cols(d, c(column = "欄位", description = "說明", unit = "單位",
                   role = "角色", type = "型態", missing = "缺失數"))
}

# Each check gets a Chinese name and a sentence on what its number means.
data_checks_display <- function(result) {
  d <- result$data_checks
  if (!is.data.frame(d) || !nrow(d)) return(d)
  info <- list(
    observations = c("觀測筆數", "分析用的資料筆數"),
    complete_observations = c("完整觀測筆數",
      "X 沒有空格的筆數；假設診斷、描述統計與圖 1–5 只用這些資料"),
    missing_cells = c("X 空格數", "X 中空白的格數"),
    rows_with_missing = c("有空格的筆數",
      "至少有一格 X 空白的筆數；不刪除，模型在每輪練習資料內補中位數"),
    duplicate_rows = c("重複筆數",
      "Y 和 X 都與其他列完全相同的筆數"),
    class_balance_ratio = c("類別筆數比",
      "最多類 ÷ 最少類；1 代表各類一樣多，大於 1.5 時不要只看正確率"),
    classes_with_n_le_p = c("筆數太少的類別數",
      "筆數 ≤ X 個數的類別數；不是 0 時，QDA 無法估計這些類別"),
    zero_within_class_variance = c("類內沒有變化的組合數",
      "某類別在某個變數上全部同值的組合數；不是 0 時模型可能算不出來"),
    full_rank_X = c("X 沒有完全線性相依",
      "是＝未發現完全線性相依，但變數之間仍可能高度相關；否＝某個 X 可由其他 X 完全算出，LDA／QDA 無法計算"),
    iqr_outlier_rows = c("疑似離群筆數",
      "在各類別內，任一變數超出 1.5 倍四分位距（IQR）的筆數；只標記，全部保留"))
  K <- nrow(result$class_counts)
  p <- length(result$settings$x)
  out <- lapply(seq_len(nrow(d)), function(i) {
    ck <- d$check[i]; note <- d$note[i]; value <- d$value[i]
    if (is.null(info[[ck]]))
      return(data.frame(檢查項目 = ck, 結果 = format(value), 說明 = note,
                        check.names = FALSE))
    detail <- switch(ck,
      observations = paste0(K, " 類、", p, " 個 X"),
      missing_cells = note,
      classes_with_n_le_p = sub("^QDA not estimable for: ", "", note),
      zero_within_class_variance = note,
      duplicate_rows = if (isTRUE(value > 0)) paste0("數值相同不必然代表同一個體或標本；",
        "本次保留，是否為重複登錄仍需依原始紀錄確認") else "本資料沒有",
      "")
    if (is.na(detail)) detail <- ""
    shown <- if (ck == "full_rank_X") ifelse(value == 1, "是", "否") else format(value)
    data.frame(檢查項目 = info[[ck]][1], 結果 = shown,
               說明 = if (nzchar(detail)) paste0(info[[ck]][2], "（", detail, "）")
                      else info[[ck]][2],
               check.names = FALSE)
  })
  do.call(rbind, out)
}

metrics_display <- function(result) {
  metric <- result$metrics
  metric$model <- toupper(metric$model)
  metric$status <- c(complete = "驗證完整", infeasible = "無法估計",
                     cv_failed = "驗證不完整")[metric$status]
  for (nm in intersect(c("coverage", "accuracy", "macro_recall",
                         "macro_recall_sd", "macro_f1"), names(metric)))
    metric[[nm]] <- pct_text(metric[[nm]])
  metric$correct_mean <- ifelse(is.na(metric$correct_mean), "未計算",
                                sprintf("%.1f", metric$correct_mean))
  rename_cols(metric, c(model = "模型", status = "模型狀態",
    coverage = "預測涵蓋率", accuracy = "平均正確率",
    macro_recall = "平均 Macro recall", macro_recall_sd = "Macro recall 標準差",
    macro_f1 = "平均 Macro-F1", correct_mean = "平均答對筆數",
    correct_min = "最少答對", correct_max = "最多答對"))
}

class_recall_display <- function(result) {
  d <- result$class_recall
  if (!is.data.frame(d) || !nrow(d)) return(d)
  d$model <- toupper(d$model)
  for (nm in intersect(c("recall", "precision", "f1"), names(d)))
    d[[nm]] <- pct_text(d[[nm]])
  rename_cols(d, c(model = "模型", class = "類別", recall = "召回率 recall",
                   precision = "精確率 precision", f1 = "F1"))
}

diagnostics_display <- function(result) {
  diag <- result$diagnostics
  labels <- c(
    computed = "已計算", computed_with_warning = "已計算，有警告",
    missing_dependency = "未執行：缺少套件", not_estimable = "無法估計",
    unsupported = "本版未支援", error = "計算失敗")
  known <- diag$status %in% names(labels)
  diag$status[known] <- unname(labels[diag$status[known]])
  names(diag) <- c("診斷", "組別", "樣本數", "統計量", "p 值", "狀態", "備註")
  diag
}

data_summary_display <- function(result) {
  d <- result$data_summary
  if (!is.data.frame(d) || !nrow(d)) return(d)
  rename_cols(d, c(class = "類別", variable = "變數", n = "筆數", mean = "平均數",
                   sd = "標準差", min = "最小值", median = "中位數", max = "最大值"))
}

# One worked example from the confusion matrix of the selected model.
confusion_example <- function(result) {
  cm <- result$confusion_matrix
  if (!is.data.frame(cm) || !nrow(cm)) return(character())
  m <- if (!is.null(result$cv_leader) && result$cv_leader %in% cm$model)
    result$cv_leader else cm$model[1]
  z <- cm[cm$model == m & as.character(cm$truth) != as.character(cm$predicted), ]
  if (!nrow(z) || max(z$row_pct, na.rm = TRUE) <= 0)
    return(paste0(toupper(m), " 沒有任何誤判，對角線都是 100%。"))
  k <- which.max(z$row_pct)
  v <- sprintf("%.1f%%", 100 * z$row_pct[k])
  paste0("例如 ", toupper(m), " 表中，", z$truth[k], " 那一列、", z$predicted[k],
         " 那一欄是 ", v, "：真實是 ", z$truth[k], " 的資料，有 ", v, " 被誤判成 ",
         z$predicted[k], "（這是誤判比例最高的一格）。")
}

# How-to-read notes placed above each table.
report_guides <- function(result) {
  N <- sum(result$class_counts$n)
  R <- result$settings$repeats
  folds <- result$settings$folds
  list(
    class_counts = paste0("下表是每個類別（Y）各有幾筆資料。各類筆數差很多時，",
      "整體正確率會被筆數多的類別主導，所以本報告以 Macro recall 為主要指標。"),
    variable_info = paste0("下表說明每個欄位的意思與單位。「角色」欄：X（預測變數）",
      "是用來判斷的欄位，Y（分類變數）是要分的類別。"),
    data_checks = paste0("以下只做標記與說明，不刪除或轉換任何資料（缺失值處理另見下一節）。",
      "每一列是一項檢查：「結果」是檢查出的數字，「說明」寫這個數字代表什麼。"),
    metrics = paste0("怎麼看這張表：每一列是一個方法，分數都是 ", R, " 次重複的平均。",
      "「平均正確率」是全部資料猜對的比例；",
      "「平均 Macro recall」是先算每一類各自猜對的比例再平均，每一類一樣重要，是本報告選方法的依據；",
      "「Macro recall 標準差」越小，代表換一種隨機分組結果也差不多（它是切分之間的變動，不是信賴區間）；",
      "「平均 Macro-F1」同時考慮漏抓和誤判；",
      "「平均答對筆數」是每次重複 ", N, " 筆中答對幾筆，「最少／最多答對」是 ", R,
      " 次中的範圍；「預測涵蓋率」是成功給出預測的比例，正常是 100%。"),
    repeats = paste0("每一次重複都把資料重新隨機分成 ", folds, " 份；下表看得出分數隨切分變動多少，",
      "以及哪個方法較常勝出。每一列是一次重複：「答對」是 ", N,
      " 筆中答對幾筆，「較高」是該次 Macro recall 較高的方法。"),
    class_recall = paste0("把 ", R, " 次重複的預測合併計算（每筆觀測被預測 ", R,
      " 次）。每一列是一個方法在一個類別的表現：",
      "「召回率」是這一類的資料有多少比例被正確認出；",
      "「精確率」是被預測成這一類的資料中，真的屬於這一類的比例；",
      "「F1」是兩者的綜合。越接近 100% 越好；某一類特別低，代表這一類容易被認錯。"),
    confusion = c(paste0("混淆矩陣以列百分比呈現：每一列是一個真實類別，",
      "數字是被預測成各類別的比例，每一列加起來是 100%。",
      "對角線（列與欄是同一類）就是該類別的召回率，其他格子是誤判的去向。"),
      confusion_example(result)),
    diagnostics = paste0("Mardia 檢查每一類的資料是否接近多元常態（鐘形）：",
      "Skewness（偏態）看分布是否不對稱，Kurtosis（峰態）看尾巴是否太厚或太薄。",
      "Box's M 檢查各類資料的散開程度與方向（共變異數矩陣）是否相同，",
      "是所有類別一起做的檢定（組別 ALL）。",
      "p 值 < 0.05 稱為「顯著」，表示有證據違反該假設；p ≥ 0.05 表示沒有看到明顯違反。"),
    diag_cases = paste0("依 Mardia 和 Box's M 是否顯著分成四種情況，",
      "標「← 本資料」的那一列就是本資料的情況。"),
    data_summary = paste0("每一列是某一類別在某個變數上的摘要：平均數、",
      "標準差（資料散開的程度）、最小值、中位數、最大值。",
      "同一個變數在不同類別的平均數差越多、標準差越小，單看這個變數時越能區分類別；",
      "這是逐一變數的描述，不代表該變數在多變數模型中的獨立貢獻。"),
    figures = "每張圖下方都有怎麼看的說明。"
  )
}

# Title and how-to-read note placed below each figure.
figure_notes <- function(result, figures = result$figures) {
  R <- result$settings$repeats
  imputed <- isTRUE(result$settings$missing_cells > 0)
  dc <- result$data_checks
  ratio <- if (is.data.frame(dc) && any(dc$check == "class_balance_ratio"))
    dc$value[dc$check == "class_balance_ratio"][1] else NA
  ratio_text <- if (is.finite(ratio)) paste0("本資料最多與最少類別的筆數比是 ",
    trimws(formatC(ratio, digits = 3, format = "g")),
    if (ratio <= 1.5) "，算是平衡。" else "，不太平衡，要以 Macro recall 為主、不要只看正確率。")
    else ""
  lda <- result$candidate_models[["lda"]]
  ld_text <- if (is.null(lda)) "" else {
    prop <- lda$svd^2 / sum(lda$svd^2)
    if (length(prop) >= 2L) paste0(
      "LDA 把所有 X 組合成新的座標軸（判別軸 LD1、LD2），讓不同類別盡量分開。",
      "每個點是一筆資料，顏色代表類別。括號裡的百分比是各判別軸的相對判別資訊比例",
      "（該軸的組間對組內變異比，占所有判別軸總和的比例）；本資料 LD1 占 ",
      sprintf("%.1f%%", 100 * prop[1]),
      if (prop[1] >= .9) "，判別資訊幾乎都集中在 LD1 這個方向（橫軸）" else "",
      "。這個百分比不是分類正確率，也不是原始資料的變異解釋率（與主成分分析不同）。",
      "點群分得越開、重疊越少，越容易分類；點群重疊的地方，通常也是容易誤判的類別。")
    else paste0("只有一條判別軸 LD1。LDA 把所有 X 組合成這條軸，讓不同類別盡量分開；",
      "每條曲線是一個類別在 LD1 上的分布。曲線分得越開、重疊越少，越容易分類；",
      "重疊的區域就是容易誤判的資料。")
  }
  models <- if (is.data.frame(result$confusion_matrix))
    toupper(unique(result$confusion_matrix$model)) else character()
  rt <- range_text(result)
  notes <- list(
    "01_class_counts.png" = c("各類別筆數", paste0(
      "長條高度是每個類別的筆數（數字標在長條內），用來看各類是否差不多多。", ratio_text,
      if (imputed) "這張圖只算完整觀測，所以筆數可能比上面的表少。" else "")),
    "02_boxplots.png" = c("各變數的分組箱形圖", paste0(
      "每一小張是一個變數，每個箱子是一個類別。箱子涵蓋中間 50% 的資料，",
      "中間的粗線是中位數，上下延伸的線是一般範圍，線外的圓點是可能的離群值。",
      "同一張小圖中，不同類別的箱子位置差越多、重疊越少，代表單看這個變數時越能把類別分開；",
      "箱子高矮差很多，代表各類散開的程度不同。",
      if (length(rt)) paste0("本資料", if (imputed) "（完整觀測）" else "", "：",
        paste(rt, collapse = "；"), "。") else "",
      "這是逐一變數的比較，不代表該變數在多變數模型中的獨立貢獻最大。")),
    "03_qq_by_class.png" = c("各類別、各變數的常態 Q-Q 圖", paste0(
      "每一列是一個類別、每一欄是一個變數。橫軸是「如果資料是常態分布，理論上該在的位置」，",
      "縱軸是實際資料。點大致貼著灰色直線，代表這個變數在這一類接近常態（鐘形）；",
      "兩端明顯往上翹或往下彎，代表偏態或尾巴較厚。",
      "點排成一階一階，可能與量測精度有關（例如只量到小數一位），不能僅據此判定是否符合常態。",
      "這張圖一次只看一個變數；所有變數合起來是否接近多元常態，",
      "需綜合檢定、圖形與資料特性評估，不能僅由單一檢定判定。")),
    "04_pairs.png" = c("變數兩兩散布圖", paste0(
      "每一小格是兩個變數的散布圖：該格所在列的變數當縱軸、所在欄的變數當橫軸",
      "（變數名稱寫在對角線上），點的顏色代表類別（圖例在下方）。",
      "不同顏色的點群分得越開，代表這兩個變數合起來越能區分類別。",
      "各類點群的形狀、大小、傾斜方向越相似，LDA「各類散開方式相同」的假設越合理；",
      "差很多時，QDA 可能較好。")),
    "05_lda_scores.png" = c("LDA 判別分數圖", paste0(ld_text, "本圖以",
      if (imputed) "完整觀測" else "全部資料",
      "配適，只用來描述類別分得多開，不代表驗證表現。")),
    "06_confusion_matrix.png" = c("混淆矩陣熱圖", paste0(
      "內容與上面的混淆矩陣表相同（", R, " 次重複合併、列百分比），用顏色深淺表示比例",
      if (length(models) > 1L) paste0("；左右依序是 ", paste(models, collapse = "、")) else "",
      "。縱軸是真實類別、橫軸是預測類別；左上到右下的對角線是各類答對的比例，越深越好；",
      "對角線以外有顏色的格子就是誤判的去向。")))
  if (!length(figures))
    return(data.frame(file = character(), title = character(), text = character()))
  titles <- vapply(figures, function(f)
    if (is.null(notes[[f]])) f else notes[[f]][1], "")
  texts <- vapply(figures, function(f)
    if (is.null(notes[[f]])) "" else notes[[f]][2], "")
  # Number by file prefix (01_ -> 圖 1) so "圖 1–5" always means the same plots,
  # even when the pairs plot is skipped.
  num <- ifelse(grepl("^[0-9]{2}_", figures),
                as.integer(substr(figures, 1, 2)), seq_along(figures))
  data.frame(file = figures, title = paste0("圖 ", num, "：", titles),
             text = texts, row.names = NULL, stringsAsFactors = FALSE)
}

# ---- Findings drawn from this data set (notes and conclusion) --------------
# Every sentence below is built from numbers already computed in `result`;
# nothing here fits a new model or runs a new test.

# Dictionary description of a variable ("花瓣長度"), or its name.
var_label <- function(result, v) {
  vi <- result$variable_info
  d <- rep(NA_character_, length(v))
  if (is.data.frame(vi) && all(c("column", "description") %in% names(vi)))
    d <- as.character(vi$description[match(v, vi$column)])
  d <- trimws(sub("[（(].*$", "", d))
  ifelse(is.na(d) | !nzchar(d), v, d)
}
var_label_full <- function(result, v) {
  d <- var_label(result, v)
  ifelse(d == v, v, paste0(d, "（", v, "）"))
}

purpose_text <- function(result) {
  x <- result$settings$x
  y <- result$settings$y
  lev <- as.character(result$class_counts$class)
  xl <- var_label(result, x)
  x_part <- if (length(xl) > 8L) paste0(paste(head(xl, 8L), collapse = "、"), " 等")
    else paste(xl, collapse = "、")
  yd <- var_label(result, y)
  paste0("本分析以 ", length(x), " 個預測變數（", x_part, "）預測 ",
         if (yd == y) paste0(y, " 的 ") else paste0(y, "（", yd, "）的 "), length(lev), " 個類別",
         if (length(lev) <= 8L) paste0("（", paste(lev, collapse = "、"), "）") else "",
         "，並比較 LDA 與 QDA 的分類表現。")
}

# Per-variable ranges (min/max) from the group summaries: which class does
# not overlap the others on a variable, and which class pairs overlap on all.
range_findings <- function(result) {
  s <- result$data_summary
  empty <- list(sep = data.frame(), overlap_all = list(), nvar = 0L)
  if (!is.data.frame(s) || !nrow(s)) return(empty)
  cls <- unique(as.character(s$class))
  vars <- unique(as.character(s$variable))
  sep <- do.call(rbind, lapply(cls, function(cl) do.call(rbind, lapply(vars, function(v) {
    own <- s[s$class == cl & s$variable == v, ]
    oth <- s[s$class != cl & s$variable == v, ]
    if (nrow(own) != 1L || !nrow(oth)) return(NULL)
    if (own$max < min(oth$min))
      data.frame(class = cl, variable = v, side = "below", own = own$max,
                 other = min(oth$min), stringsAsFactors = FALSE)
    else if (own$min > max(oth$max))
      data.frame(class = cl, variable = v, side = "above", own = own$min,
                 other = max(oth$max), stringsAsFactors = FALSE)
    else NULL
  }))))
  pairs <- if (length(cls) >= 2L) utils::combn(cls, 2L, simplify = FALSE) else list()
  overlap_all <- Filter(function(pr) all(vapply(vars, function(v) {
    a <- s[s$class == pr[1] & s$variable == v, ]
    b <- s[s$class == pr[2] & s$variable == v, ]
    a$max >= b$min && b$max >= a$min
  }, logical(1))), pairs)
  list(sep = if (is.null(sep)) data.frame() else sep,
       overlap_all = overlap_all, nvar = length(vars))
}

range_text <- function(result) {
  rf <- range_findings(result)
  num <- function(v) trimws(formatC(v, digits = 4, format = "g"))
  out <- character()
  if (nrow(rf$sep)) for (cl in unique(rf$sep$class)) {
    z <- rf$sep[rf$sep$class == cl, ]
    ex <- z[1, ]
    vl <- paste(var_label_full(result, z$variable), collapse = "、")
    out <- c(out, paste0(cl, " 在", if (grepl("^[A-Za-z0-9]", vl)) " " else "", vl,
      if (grepl("[A-Za-z0-9._]$", vl)) " " else "",
      "的範圍與其他類別完全沒有重疊（例如 ", ex$variable, "：", cl,
      if (ex$side == "below") paste0(" 最大 ", num(ex$own), "，其他類別最小 ", num(ex$other))
      else paste0(" 最小 ", num(ex$own), "，其他類別最大 ", num(ex$other)), "）"))
  }
  if (length(rf$overlap_all)) {
    pr <- head(rf$overlap_all, 3L)
    out <- c(out, paste0(paste(vapply(pr, paste, "", collapse = " 與 "), collapse = "、"),
      " 在全部 ", rf$nvar, " 個變數上的範圍都有重疊"))
  }
  out
}

# Class-level results of the models that completed every repeat.
class_findings <- function(result) {
  m <- result$metrics
  done <- as.character(m$model[m$status == "complete"])
  cr <- result$class_recall
  cr <- if (is.data.frame(cr)) cr[cr$model %in% done, , drop = FALSE] else data.frame()
  cm <- result$confusion_matrix
  cm <- if (is.data.frame(cm)) cm[cm$model %in% done, , drop = FALSE] else data.frame()
  lev <- result$settings$class_levels
  perfect <- lev[vapply(lev, function(cl) {
    z <- cr[cr$class == cl, , drop = FALSE]
    nrow(z) > 0L && all(abs(z$recall - 1) < 1e-12 & abs(z$precision - 1) < 1e-12)
  }, logical(1))]
  list(done = done, cr = cr, cm = cm, perfect = perfect, lev = lev)
}

# Three kinds of sentence: perfectly classified classes, where the errors
# go, and per-class recall differences between the two models.
misclass_text <- function(result) {
  f <- class_findings(result)
  out <- list(perfect = character(), errors = character(), recall = character())
  if (!length(f$done) || !nrow(f$cm)) return(out)
  who <- if (length(f$done) == 2L) "兩種方法" else toupper(f$done)
  if (length(f$perfect))
    out$perfect <- paste0(who, "對 ", paste(f$perfect, collapse = "、"),
                          " 的召回率與精確率皆為 100%。")
  truth <- as.character(f$cm$truth)
  pred <- as.character(f$cm$predicted)
  off <- truth != pred & f$cm$row_pct > 0
  if (!any(off)) {
    out$errors <- paste0(who, "都沒有任何誤判。")
  } else {
    key <- ifelse(truth < pred, paste(truth, pred, sep = "\u0001"),
                  paste(pred, truth, sep = "\u0001"))[off]
    tot <- tapply(f$cm$row_pct[off], key, sum)
    pr <- strsplit(names(tot)[which.max(tot)], "\u0001", fixed = TRUE)[[1]]
    cell <- function(mod, a, b) sprintf("%.1f%%", 100 * f$cm$row_pct[
      f$cm$model == mod & truth == a & pred == b])
    detail <- paste(vapply(f$done, function(mod) paste0(toupper(mod), "：",
      pr[1], " 被誤判為 ", pr[2], " ", cell(mod, pr[1], pr[2]), "，",
      pr[2], " 被誤判為 ", pr[1], " ", cell(mod, pr[2], pr[1])), ""), collapse = "；")
    out$errors <- paste0(if (length(tot) == 1L) "誤判全部發生在 " else "誤判比例最高的是 ",
      pr[1], " 與 ", pr[2], " 之間（", detail, "）",
      if (length(tot) > 1L) "；其他誤判見混淆矩陣" else "", "。")
  }
  if (length(f$done) == 2L) {
    rec <- function(mod, cl) f$cr$recall[f$cr$model == mod & f$cr$class == cl]
    diffs <- unlist(lapply(f$lev, function(cl) {
      a <- rec("lda", cl); b <- rec("qda", cl)
      if (length(a) != 1L || length(b) != 1L || abs(a - b) < 1e-12) return(NULL)
      hi <- if (a > b) c("LDA", "QDA") else c("QDA", "LDA")
      paste0(hi[1], " 對 ", cl, " 的召回率為 ", pct_text(max(a, b)), "，高於 ",
             hi[2], " 的 ", pct_text(min(a, b)))
    }))
    out$recall <- if (length(diffs)) paste0(paste(diffs, collapse = "；"), "。")
      else "兩種方法各類別的召回率相同。"
  }
  out
}

# What the diagnostics say, and how they relate to the validated choice.
diagnostic_outcome <- function(result, alpha = 0.05) {
  d <- result$diagnostics
  ok <- d$status %in% c("computed", "computed_with_warning") & is.finite(d$p_value)
  mard <- grepl("^Mardia", d$test)
  box <- d$test == "Box M"
  fp <- function(p) ifelse(p < 0.001, formatC(p, digits = 2, format = "e"),
                           sprintf("%.3f", p))
  full <- character(); short <- character()
  if (any(mard & ok)) {
    sig <- unique(d$group[mard & ok & d$p_value < alpha])
    if (!length(sig)) {
      r <- range(d$p_value[mard & ok])
      full <- c(full, paste0("Mardia 檢定：", length(unique(d$group[mard & ok])),
        " 個類別的偏態與峰態均未達顯著（p 值 ", fp(r[1]), "–", fp(r[2]),
        "）。這只表示沒有看到明顯違反，不代表已證明多元常態成立。"))
      short <- c(short, "Mardia 檢定在各類別均未達顯著，但這不代表已證明多元常態成立")
    } else {
      full <- c(full, paste0("Mardia 檢定：", paste(sig, collapse = "、"),
        " 達顯著（p < ", alpha, "），這些類別的多元常態假設受到質疑，算出的機率要保守解讀。"))
      short <- c(short, paste0("Mardia 檢定在 ", paste(sig, collapse = "、"),
        " 達顯著，這些類別的多元常態假設受到質疑"))
    }
  } else {
    full <- c(full, "Mardia 檢定未能完成，無法評估多元常態。")
    short <- c(short, "Mardia 檢定未能完成")
  }
  box_sig <- NA
  if (any(box & ok)) {
    i <- which(box & ok)[1]
    box_sig <- d$p_value[i] < alpha
    stat <- paste0("χ² 近似統計量 ", trimws(formatC(d$statistic[i], digits = 4, format = "g")),
      if (grepl("^df = ", d$note[i])) paste0("，", d$note[i]) else "",
      "，p = ", fp(d$p_value[i]))
    full <- c(full, if (box_sig) paste0("Box's M 顯著（", stat,
        "），表示「各類共變異數矩陣相同」的假設受到質疑；但 Box's M 對非常態敏感，",
        "樣本量增加時也可能偵測到較小的差異，應結合驗證表現解讀。")
      else paste0("Box's M 未達顯著（", stat,
        "），沒有看到各類共變異數不同的明顯證據，但也不代表已證明相同。"))
    short <- c(short, if (box_sig) "Box's M 顯著，「各類共變異數相同」的假設受到質疑"
                      else "Box's M 未達顯著")
  } else {
    full <- c(full, "Box's M 未能完成，無法評估各類共變異數是否相同。")
    short <- c(short, "Box's M 未能完成")
  }
  m <- result$metrics
  both <- all(c("lda", "qda") %in% m$model[m$status == "complete"])
  sel <- if (is.null(result$selection)) "none" else result$selection
  link <- character(); link_short <- character()
  if (both && !is.na(box_sig)) {
    l <- pct_text(m$macro_recall[m$model == "lda"])
    q <- pct_text(m$macro_recall[m$model == "qda"])
    K <- length(result$settings$class_levels)
    p <- length(result$settings$x)
    extra <- (K - 1L) * p * (p + 1L) / 2L
    link <- switch(sel,
      lda_better = if (box_sig) paste0("但 QDA 的實際驗證表現沒有更好（平均 Macro recall：QDA ",
          q, "，LDA ", l, "），因此仍依交叉驗證結果選擇 LDA。一種可能的解釋是：QDA 比 LDA 多估計 ",
          extra, " 個參數，允許各類共變異數不同的好處，可能被估計誤差抵消。")
        else paste0("LDA 的驗證表現也較高（平均 Macro recall：LDA ", l, "，QDA ", q,
          "），與等共變異數假設沒有明顯衝突。"),
      qda_better = if (box_sig) paste0("QDA 的驗證表現也較高（平均 Macro recall：QDA ", q,
          "，LDA ", l, "），方向與 Box's M 一致；但選模仍依交叉驗證結果，不依 Box's M。")
        else paste0("雖然 Box's M 未達顯著，QDA 的驗證表現較高（平均 Macro recall：QDA ",
          q, "，LDA ", l, "）；選模依交叉驗證結果。"),
      tie_prefer_lda = paste0("兩種方法的驗證表現相同（平均 Macro recall ", l,
        "），依規則選較簡單的 LDA。"),
      character())
    link_short <- switch(sel,
      lda_better = if (box_sig) "但 QDA 的實際驗證表現沒有更好，因此仍依交叉驗證結果選擇 LDA；診斷只提醒 LDA 的等共變異數假設可能不完全成立"
        else "LDA 的驗證表現也較高，與等共變異數假設沒有明顯衝突",
      qda_better = if (box_sig) "QDA 的驗證表現也較高，方向與 Box's M 一致，但選模仍依交叉驗證"
        else "QDA 的驗證表現較高，選模依交叉驗證",
      tie_prefer_lda = "兩種方法的驗證表現相同，依規則選 LDA",
      character())
  }
  list(full = c(full, link, "診斷用來說明假設是否合理，不改變選模規則。"),
       short = paste0(paste(c(short, link_short), collapse = "。"),
                      "。詳細數字見「假設診斷」一節。"))
}

# Final chapter: performance, why this model, where errors go, what the
# diagnostics mean for the result, and limitations.
conclusion_lines <- function(result) {
  m <- result$metrics
  N <- sum(result$class_counts$n)
  R <- result$settings$repeats
  sel <- if (is.null(result$selection)) "none" else result$selection
  done <- as.character(m$model[m$status == "complete"])
  get <- function(model, col) m[[col]][m$model == model]
  wins <- result$selection_wins
  both <- all(c("lda", "qda") %in% done)
  winner <- if (sel %in% c("lda_better", "qda_better", "only_feasible"))
    result$cv_leader else if (sel == "tie_prefer_lda") "lda" else NA_character_
  one <- function(model) paste0(toupper(model), " 平均答對 ",
    sprintf("%.1f", get(model, "correct_mean")), " 筆（範圍 ", get(model, "correct_min"),
    "–", get(model, "correct_max"), "）")
  perf <- character(); why <- character(); small <- TRUE
  if (both) {
    w <- if (is.na(winner)) "lda" else winner
    o <- setdiff(c("lda", "qda"), w)
    dpp <- 100 * (get(w, "macro_recall") - get(o, "macro_recall"))
    dc <- get(w, "correct_mean") - get(o, "correct_mean")
    small <- abs(dpp) < 2
    perf <- c(
      paste0("LDA 平均 Macro recall 為 ", pct_text(get("lda", "macro_recall")), "，QDA 為 ",
        pct_text(get("qda", "macro_recall")),
        if (abs(dpp) > 1e-9) paste0("，相差 ", sprintf("%.2f", abs(dpp)), " 個百分點") else "，兩者相同",
        if (small) "，兩種方法的表現接近。" else "。"),
      paste0("以答對筆數來看，", one(w), "，", one(o), "；",
        if (dc > 1e-9) paste0(toupper(w), " 每次 ", N, " 筆平均多答對 ", sprintf("%.1f", dc), " 筆。")
        else if (dc < -1e-9) paste0(toupper(w), " 平均答對筆數反而少 ", sprintf("%.1f", -dc),
          " 筆（Macro recall 讓每一類一樣重要，與答對筆數不一定同向）。")
        else "兩者平均答對筆數相同。"),
      paste0(R, " 次重複中，", toupper(w), " 較高 ", wins[[w]], " 次、相同 ", wins[["tie"]],
        " 次、較低 ", wins[[o]], " 次。"))
  } else if (length(done) == 1L) {
    perf <- paste0("只有 ", toupper(done), " 完成全部 ", R, " 次驗證：平均 Macro recall ",
      pct_text(get(done, "macro_recall")), "，", one(done), "。另一個方法驗證不完整或無法估計",
      "（原因見「模型警告與錯誤」），所以無法比較兩種方法。")
  } else {
    perf <- "沒有方法完成全部驗證，本報告不提供分類表現的結論；原因見「模型警告與錯誤」。"
  }
  if (sel %in% c("lda_better", "qda_better")) {
    W <- toupper(winner); O <- toupper(setdiff(c("lda", "qda"), winner))
    why <- c(paste0(W, " 在本次交叉驗證中", if (small) "略優" else "較優",
        "，因此依既定規則（", R, " 次的平均 Macro recall 較高者勝出）選用 ", W, "。"),
      paste0("這是本資料、本次驗證設計下的比較，不是統計檢定；不代表 ", W, " 顯著優於 ", O,
        "，也不代表 ", W, " 在其他資料上普遍較好。",
        if (wins[[tolower(O)]] == 0L) paste0(R, " 次重複中 ", O, " 沒有一次較高，方向一致。")
        else paste0(O, " 也在 ", wins[[tolower(O)]], " 次重複中較高，結果會隨切分改變。")))
  } else if (sel == "tie_prefer_lda") {
    why <- paste0("兩者 ", R, " 次的平均 Macro recall 相同，依規則選較簡單的 LDA；這不代表兩者完全等效。")
  } else if (sel == "only_feasible") {
    why <- paste0("只有 ", toupper(winner), " 完成驗證，所以選 ", toupper(winner),
                  "；這不是兩種方法比較後的結果。")
  } else {
    why <- "本次沒有完整的驗證結果，不選擇方法。"
  }
  if (!is.null(result$final_model) && !is.na(winner))
    why <- c(why, paste0("選定後以全部 ", N, " 筆資料重新配適 ", toupper(winner),
                         "，作為最終模型（final_model.rds）。"))
  mt <- misclass_text(result)
  rt <- range_text(result)
  imputed <- isTRUE(result$settings$missing_cells > 0)
  mis <- c(paste0("- ", c(mt$perfect, mt$errors, mt$recall)),
    if (length(rt)) paste0("- 描述統計", if (imputed) "（以完整觀測計算）" else "", "中，",
      paste(rt, collapse = "；"),
      "。這些是逐一變數的描述，有助於理解誤判集中在哪些類別，但不代表個別變數在模型中的獨立貢獻最大。"))
  if (length(mis) == 0L || !length(done)) mis <- "沒有完整的驗證結果，無法判斷哪些類別容易誤判。"
  dc_rows <- result$data_checks
  dup <- if (is.data.frame(dc_rows) && any(dc_rows$check == "duplicate_rows"))
    dc_rows$value[dc_rows$check == "duplicate_rows"][1] else 0
  lim <- paste0("- ", c(
    if (length(done)) c(
      paste0(R, " 次重複使用的是同一批 ", N, " 筆觀測，只是切分方式不同，不代表有 ",
             format(N * R, big.mark = ","), " 筆獨立樣本。"),
      "Macro recall 標準差是不同切分之間的變動，不是新資料預測表現的信賴區間。",
      "同一組交叉驗證同時用於比較與選擇模型，分數不是獨立測試表現，可能略為樂觀。"),
    "尚無外部資料驗證；推廣到其他資料來源或採集條件，仍需進一步確認。",
    "本報告沒有計算變數重要性或判別係數的不確定性，不宜據此判斷哪個變數最重要。",
    if (isTRUE(dup > 0)) paste0("資料中有 ", dup, " 筆與其他列的 Y 和 X 完全相同，已保留；",
      "數值相同不必然代表同一個體或標本，是否為重複登錄仍需依原始紀錄確認。"),
    if (imputed) "X 的空格以各訓練折的中位數補值；補值方式可能影響結果，本報告未與其他補值方式比較。",
    "本版處理獨立觀測與數值型預測變數，使用重複分層交叉驗證（固定 5 折、10 次）。",
    "各類別診斷有多個檢定，p 值以診斷用途呈現，未進行多重比較校正。"))
  c("## 結論與討論", "",
    "### 兩種方法的分類表現", "", perf, "",
    paste0("### ", if (is.na(winner)) "選模結果" else paste0("為什麼選擇 ", toupper(winner))), "",
    why, "",
    "### 哪些類別容易誤判", "", mis, "",
    "### 假設診斷對結果的影響", "", diagnostic_outcome(result)$short, "",
    "### 限制", "", lim)
}

make_report <- function(result) {
  needed <- c("settings", "metrics", "status", "model_status",
              "diagnostic_status", "diagnostics", "data_summary", "class_counts")
  if (!all(needed %in% names(result))) stop("Use run_pipeline() to assemble results.")
  state <- function(x) switch(x,
    completed = "完成", partial = "部分完成",
    blocked = "待處理", failed = "失敗", x)
  guide <- report_guides(result)
  mt <- misclass_text(result)
  observed_recall <- if (length(c(mt$perfect, mt$recall)))
    paste0("本資料：", paste(c(mt$perfect, mt$recall), collapse = "")) else character()
  leader <- selection_text(result)
  dn <- diagnostic_notes(result$diagnostics)
  issues <- if (is.data.frame(result$issues) && nrow(result$issues)) {
    markdown_table(result$issues)
  } else {
    "模型執行過程沒有記錄到警告或錯誤；診斷未完成原因另見診斷表。"
  }
  intro <- character()
  if (!is.null(result$dataset_description) &&
      nzchar(result$dataset_description))
    intro <- c(intro, paste0("- 資料說明：", result$dataset_description))
  variable_block <- if (is.data.frame(result$variable_info) &&
                        nrow(result$variable_info)) {
    c("", "### 變數介紹", "", guide$variable_info, "",
      markdown_table(variable_info_display(result)))
  } else character()
  check_block <- if (is.data.frame(result$data_checks) &&
                     nrow(result$data_checks)) {
    c("## 資料檢查", "", guide$data_checks, "",
      markdown_table(data_checks_display(result)), "")
  } else character()
  impute_block <- if (identical(result$settings$impute, "median") &&
                       is.data.frame(result$imputation) && nrow(result$imputation)) {
    c("## 缺失值處理", "",
      paste0("X 共有 ", result$settings$missing_cells, " 個空格。依使用者決定，",
             "在每一份練習資料內計算各欄中位數，再用這個中位數補該輪的練習與考試資料；",
             "考試資料不參與計算，避免資料洩漏。最終模型以全部資料的中位數補值。"),
      "假設診斷、描述統計與探索圖（圖 1–5）只使用沒有空格的完整觀測；混淆矩陣來自交叉驗證，包含全部觀測。",
      paste0("下表摘要 ", nrow(result$imputation) / length(result$settings$x),
             " 輪練習資料的中位數：每個變數在各輪用到的中位數最小與最大值，",
             "以及每輪平均補了幾格；每一輪的完整紀錄見 `imputation.csv`。"), "",
      markdown_table(imputation_summary(result)), "")
  } else character()
  repeat_block <- if (is.data.frame(result$repeat_metrics) &&
                      nrow(result$repeat_metrics)) {
    c("## 每一次重複的結果", "", guide$repeats, "",
      markdown_table(repeat_table(result)), "")
  } else character()
  cm_block <- if (is.data.frame(result$confusion_matrix) &&
                  nrow(result$confusion_matrix)) {
    c(guide$confusion[1], "", guide$confusion[-1], misclass_text(result)$errors, "",
      unlist(lapply(unique(result$confusion_matrix$model), function(m) c(
        paste0("### 混淆矩陣：", toupper(m)), "",
        markdown_table(confusion_percent(result, m)), ""))))
  } else character()
  figure_block <- if (length(result$figures)) {
    fn <- figure_notes(result)
    c("## 圖", "", guide$figures, "",
      unlist(lapply(seq_len(nrow(fn)), function(i) c(
        paste0("### ", fn$title[i]), "",
        paste0("- `figures/", fn$file[i], "`"), "",
        fn$text[i], ""))))
  } else character()
  c(
    "# LDA／QDA 分析報告", "",
    purpose_text(result), "",
    paste0("整體狀態：", state(result$status), "。"),
    paste0("模型分析：", state(result$model_status),
           "；假設診斷：", state(result$diagnostic_status), "。"), "",
    "## 資料與設定", "",
    intro,
    paste0("- 分類變數：", gsub("[\r\n]+", " ", result$settings$y)),
    paste0("- 預測變數：",
           paste(gsub("[\r\n]+", " ", result$settings$x), collapse = "、")),
    paste0("- 樣本數：", sum(result$class_counts$n)),
    paste0("- 驗證方式：重複分層 ", result$settings$folds, "-fold 交叉驗證 ",
           result$settings$repeats, " 次（共 ",
           result$settings$folds * result$settings$repeats,
           " 輪練習／考試）；隨機種子：", result$settings$seed),
    "- 主要指標：Macro recall（各類別召回率的等權平均，即平衡正確率）。",
    "- 先驗機率：由各訓練折的類別比例估計。",
    if (identical(result$settings$impute, "median") && result$settings$missing_cells > 0)
      "- 資料處理：未刪除觀測值或轉換變數；缺失值依使用者同意在各練習資料內補中位數。"
    else "- 資料處理：未刪除觀測值、插補缺失值或轉換變數。",
    paste0("- 軟體：", result$settings$R, "；MASS ", result$settings$MASS,
           "；診斷：", result$settings$diagnostic_method), "",
    "### 各類別筆數", "", guide$class_counts, "",
    markdown_table(class_counts_display(result)),
    variable_block, "",
    check_block,
    impute_block,
    "## 模型驗證表現", "", guide$metrics, "",
    markdown_table(metrics_display(result)), "", leader, "",
    paste0("每次重複中，每筆觀測恰好被預測一次。",
           "必須完成全部重複的模型才參加比較，驗證不完整的模型不提供分數。"),
    "這組交叉驗證同時用於比較模型，分數不是獨立測試表現。", "",
    repeat_block,
    "## 各類別表現", "", guide$class_recall, "",
    observed_recall, "",
    markdown_table(class_recall_display(result)), "",
    cm_block,
    "## 假設診斷", "", guide$diagnostics, "",
    markdown_table(diagnostics_display(result)), "",
    "檢定不顯著不代表假設已被證明；無法估計時，不產生替代 p 值。",
    "Box's M 對非常態很敏感；樣本量增加時可能偵測到較小的差異，應結合驗證表現解讀。",
    "診斷結果用於說明適用性，不直接依單一 p 值決定 LDA／QDA。", "",
    "### Mardia 與 Box's M 的合併解讀", "", guide$diag_cases, "",
    markdown_table(dn$table), "", dn$notes, "",
    "### 對照實際驗證結果", "", diagnostic_outcome(result)$full, "",
    "## 分組描述統計", "", guide$data_summary, "",
    markdown_table(data_summary_display(result)), "",
    figure_block,
    "## 模型警告與錯誤", "", issues, "",
    conclusion_lines(result)
  )
}

# ---------------------------------------------------------------------------
# 5. Entry points
# ---------------------------------------------------------------------------

run_pipeline <- function(data, y, x, folds = 5L, repeats = 10L, seed = 2026L,
                         variable_info = NULL, dataset_description = NULL,
                         impute = c("none", "median")) {
  impute <- match.arg(impute)
  if (is.character(data) && length(data) == 1L)
    data <- read.csv(data, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  result <- run_analysis(data, y, x, folds = folds, repeats = repeats,
                         seed = seed, impute = impute)
  complete <- stats::complete.cases(data[x])
  checks <- diagnose_data(data[complete, , drop = FALSE], y, x)
  result$model_status <- result$status
  result$data_summary <- checks$data_summary
  all_y <- droplevels(factor(data[[y]]))
  result$class_counts <- data.frame(class = levels(all_y),
                                    n = as.integer(table(all_y)))
  result$data_checks <- checks$data_checks
  if (any(!complete)) {
    miss <- colSums(is.na(data[x]))
    result$data_checks <- rbind(data.frame(
      check = c("missing_cells", "rows_with_missing"),
      value = c(sum(miss), sum(!complete)),
      note = c(paste(paste0(names(miss)[miss > 0], "=", miss[miss > 0]), collapse = ", "),
               paste0("Diagnostics, summaries and figures use the ", sum(complete),
                      " complete rows; models use fold-wise median imputation.")),
      stringsAsFactors = FALSE), result$data_checks)
    result$data_checks$check[result$data_checks$check == "observations"] <-
      "complete_observations"
  }
  result$diagnostics <- checks$diagnostics
  result$settings$diagnostic_method <-
    "base R (Mardia 1970; Box's M chi-square approximation)"
  result$diagnostic_raw <- checks$raw
  if (!is.null(variable_info)) result$variable_info <- variable_info
  if (!is.null(dataset_description))
    result$dataset_description <- dataset_description
  result$diagnostic_status <- if (
    nrow(checks$diagnostics) > 0L &&
    all(checks$diagnostics$status %in% c("computed", "computed_with_warning"))
  ) "completed" else "partial"
  if (result$model_status == "completed" &&
      result$diagnostic_status != "completed") result$status <- "partial"
  result$report_lines <- make_report(result)
  result
}

# Writes tables, models and report.md. Figures are written by save_figures().
save_analysis <- function(result, out) {
  report <- make_report(result)
  if (!is.character(out) || length(out) != 1L ||
      is.na(out) || !nzchar(out)) stop("out must be one directory path.")
  if (file.exists(out)) stop("Output path already exists; choose a new path.")
  if (!dir.create(out, recursive = TRUE, showWarnings = FALSE))
    stop("Cannot create output directory.")
  manifest <- list(status = "in_progress", files = character(), error = NULL)
  manifest_path <- file.path(out, "export_status.rds")
  tryCatch({
    saveRDS(manifest, manifest_path)
    result$report_lines <- report
    saveRDS(result, file.path(out, "analysis.rds"))
    manifest$files <- c(manifest$files, "analysis.rds")
    saveRDS(result$settings, file.path(out, "settings.rds"))
    manifest$files <- c(manifest$files, "settings.rds")
    if (!is.null(result$final_model)) {
      saveRDS(result$final_model, file.path(out, "final_model.rds"))
      manifest$files <- c(manifest$files, "final_model.rds")
    }
    tables <- intersect(c(
      "metrics", "repeat_metrics", "class_recall", "fold_metrics", "fold_assignments",
      "cv_predictions", "cv_probabilities", "confusion_matrix", "issues",
      "data_checks", "imputation", "data_summary", "class_counts", "diagnostics",
      "variable_info"), names(result))
    for (nm in tables) {
      d <- result[[nm]]
      if (is.data.frame(d) && ncol(d) > 0L) {
        filename <- paste0(nm, ".csv")
        write.csv(d, file.path(out, filename),
                  row.names = FALSE, fileEncoding = "UTF-8", na = "")
        manifest$files <- c(manifest$files, filename)
      }
    }
    writeLines(enc2utf8(report), file.path(out, "report.md"), useBytes = TRUE)
    manifest$files <- c(manifest$files, "report.md")
    manifest$status <- "completed"
    saveRDS(manifest, manifest_path)
  }, error = function(e) {
    manifest$status <- "failed"
    manifest$error <- conditionMessage(e)
    try(saveRDS(manifest, manifest_path), silent = TRUE)
    stop(conditionMessage(e), call. = FALSE)
  })
  invisible(normalizePath(out, winslash = "/", mustWork = TRUE))
}

# ---------------------------------------------------------------------------
# 6. HTML report (written directly, no Quarto needed) and QMD report
# ---------------------------------------------------------------------------

base64_file <- function(path) {
  r <- as.integer(readBin(path, "raw", file.info(path)$size))
  pad <- (3L - length(r) %% 3L) %% 3L
  m <- matrix(c(r, rep(0L, pad)), nrow = 3L)
  v <- m[1, ] * 65536 + m[2, ] * 256 + m[3, ]
  idx <- rbind(v %/% 262144, (v %/% 4096) %% 64, (v %/% 64) %% 64, v %% 64)
  chars <- strsplit(paste0("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz",
                           "0123456789+/"), "")[[1L]]
  s <- chars[as.vector(idx) + 1L]
  if (pad) s[(length(s) - pad + 1L):length(s)] <- "="
  paste(s, collapse = "")
}

report_html_lines <- function(lines, fig_dir) {
  esc <- function(s) gsub("<", "&lt;", gsub("&", "&amp;", s, fixed = TRUE),
                          fixed = TRUE)
  inline <- function(s) gsub("`([^`]+)`", "<code>\\1</code>", s)
  cells <- function(s) trimws(strsplit(sub("^\\|", "", sub("\\|\\s*$", "", s)),
                                       "|", fixed = TRUE)[[1L]])
  html <- character(); i <- 1L; n <- length(lines)
  while (i <= n) {
    s <- lines[[i]]
    if (!nzchar(trimws(s))) { i <- i + 1L; next }
    if (grepl("^#{1,3} ", s)) {
      level <- nchar(sub(" .*", "", s))
      html <- c(html, sprintf("<h%d>%s</h%d>", level,
                              inline(esc(sub("^#+ ", "", s))), level))
      i <- i + 1L
    } else if (startsWith(s, "| ")) {
      block <- character()
      while (i <= n && startsWith(lines[[i]], "| ")) {
        block <- c(block, lines[[i]]); i <- i + 1L
      }
      head <- cells(block[1L])
      body <- block[-(1:2)]
      html <- c(html, "<div class=\"tbl\"><table>",
        paste0("<thead><tr>", paste0("<th>", head, "</th>", collapse = ""),
               "</tr></thead><tbody>"),
        vapply(body, function(b) paste0("<tr>",
          paste0("<td>", cells(b), "</td>", collapse = ""), "</tr>"), ""),
        "</tbody></table></div>")
    } else if (startsWith(s, "- ")) {
      items <- character()
      while (i <= n && startsWith(lines[[i]], "- ")) {
        items <- c(items, sub("^- ", "", lines[[i]])); i <- i + 1L
      }
      figs <- sub("^`figures/(.+)`$", "\\1", items)
      is_fig <- grepl("^`figures/.+`$", items) &
        file.exists(file.path(fig_dir, figs))
      if (all(is_fig)) {
        for (f in figs) html <- c(html, sprintf(
          "<figure><img alt=\"%s\" src=\"data:image/png;base64,%s\"><figcaption>%s</figcaption></figure>",
          f, base64_file(file.path(fig_dir, f)), f))
      } else {
        html <- c(html, "<ul>", paste0("<li>", inline(esc(items)), "</li>"),
                  "</ul>")
      }
    } else {
      para <- character()
      while (i <= n && nzchar(trimws(lines[[i]])) &&
             !grepl("^(#{1,3} |\\| |- )", lines[[i]])) {
        para <- c(para, lines[[i]]); i <- i + 1L
      }
      html <- c(html, paste0("<p>", inline(esc(paste(para, collapse = "<br>"))),
                             "</p>"))
      html <- gsub("&lt;br>", "<br>", html, fixed = TRUE)
    }
  }
  html
}

write_report_html <- function(result, out, fig_dir = file.path(out, "figures")) {
  lines <- if (is.null(result$report_lines)) make_report(result) else
    result$report_lines
  css <- paste(
    "body{font-family:'Segoe UI','Microsoft JhengHei','Noto Sans CJK TC',sans-serif;",
    "max-width:980px;margin:0 auto;padding:32px 24px;color:#173448;line-height:1.75;background:#fff}",
    "h1{font-size:1.9rem;border-top:5px solid #087e8b;padding-top:18px}",
    "h2{font-size:1.35rem;margin-top:2.2rem;border-bottom:1px solid #dbe5eb;padding-bottom:6px}",
    "h3{font-size:1.1rem}code{background:#f0f5f8;padding:1px 5px;border-radius:4px}",
    ".tbl{overflow-x:auto}table{border-collapse:collapse;font-size:.9rem;margin:12px 0}",
    "th,td{border:1px solid #d7e3e9;padding:6px 10px;text-align:left;vertical-align:top}",
    "th{background:#f0f5f8}figure{margin:20px 0}img{max-width:100%;border:1px solid #e2e9ee}",
    "figcaption{font-size:.85rem;color:#607487}")
  page <- c("<!doctype html>", "<html lang=\"zh-Hant\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>LDA／QDA 分析報告</title>", paste0("<style>", css, "</style>"),
    "</head><body>", report_html_lines(lines, fig_dir), "</body></html>")
  path <- file.path(out, "report.html")
  writeLines(enc2utf8(page), path, useBytes = TRUE)
  invisible(path)
}

# QMD with live R chunks. It expects run_analysis.R and input_data.csv next to
# it (the command-line pipeline copies both), so Render recomputes everything.
write_report_qmd <- function(result, out, title = "LDA／QDA 分析報告") {
  q <- function(s) if (is.null(s)) "NULL" else encodeString(s, quote = "\"")
  x_code <- paste0("c(", paste(q(result$settings$x), collapse = ", "), ")")
  chunk <- function(code, opts = character()) c("```{r}", paste0("#| ", opts),
                                               code, "```", "")
  qmd <- c(
    "---", paste0("title: ", q(title)), "lang: zh-TW",
    "output-file: report-quarto.html", "format:", "  html:",
    "    toc: true", "    toc-title: \"目錄\"", "    embed-resources: true",
    "    code-fold: show", "execute:", "  warning: false", "  message: false",
    "---", "",
    "本報告由 LDA／QDA skill 產生。按 **Render** 會重新讀取資料、重新計算、重新畫圖；",
    "每個數字都來自下方顯示的程式碼。需要的檔案（`run_analysis.R`、`input_data.csv`）",
    "已放在同一個資料夾。", "",
    chunk(c("source(\"run_analysis.R\", encoding = \"UTF-8\")",
      "data <- read.csv(\"input_data.csv\", check.names = FALSE, fileEncoding = \"UTF-8-BOM\")",
      "variable_info <- if (file.exists(\"variable_info.csv\"))",
      "  read.csv(\"variable_info.csv\", check.names = FALSE, fileEncoding = \"UTF-8\") else NULL",
      paste0("result <- run_pipeline(data, y = ", q(result$settings$y), ","),
      paste0("  x = ", x_code, ","),
      paste0("  folds = ", result$settings$folds, ", seed = ",
             result$settings$seed, ","),
      "  variable_info = variable_info,",
      paste0("  impute = ", q(if (is.null(result$settings$impute)) "none" else
                               result$settings$impute), ","),
      paste0("  dataset_description = ", q(result$dataset_description), ")"),
      "guide <- report_guides(result)"),
      "label: setup"),
    "`r purpose_text(result)`", "",
    "## 資料與設定", "",
    "`r if (!is.null(result$dataset_description)) result$dataset_description`", "",
    "- 要分類的欄位（Y）：`r result$settings$y`",
    "- 用來判斷的欄位（X）：`r paste(result$settings$x, collapse = '、')`",
    "- 樣本數：`r sum(result$class_counts$n)`",
    "- 驗證方式：把資料分成 `r result$settings$folds` 份，每次 `r result$settings$folds - 1` 份練習、1 份考試，輪流進行（分層交叉驗證）；整個過程重新隨機分組做 `r result$settings$repeats` 次；隨機種子 `r result$settings$seed`",
    "",
    "### 各類別筆數", "", "`r guide$class_counts`", "",
    chunk("knitr::kable(class_counts_display(result))"),
    chunk(c("if (!is.null(result$variable_info)) {",
            "  cat(\"### 變數介紹\\n\\n\", guide$variable_info, \"\\n\\n\")",
            "  print(knitr::kable(variable_info_display(result)))",
            "}"), "results: asis"),
    "## 資料檢查", "", "`r guide$data_checks`", "",
    chunk("knitr::kable(data_checks_display(result))"),
    chunk(c("if (identical(result$settings$impute, \"median\") && nrow(result$imputation)) {",
            "  cat(\"### 缺失值處理\\n\\n在每一份練習資料內計算中位數再補值，考試資料不參與計算。\",",
            "      \"下表是每個變數在各輪用到的中位數最小與最大值，以及每輪平均補了幾格。\\n\\n\")",
            "  print(knitr::kable(imputation_summary(result), digits = 4))",
            "  cat(\"\\n完整紀錄（每一輪、每個變數）見 imputation.csv。\\n\")",
            "}"), "results: asis"),
    "## 假設診斷（參考用）", "",
    "`r guide$diagnostics`", "",
    "這兩個檢定的結果寫進報告當參考，**不用來選擇 LDA 或 QDA**。", "",
    chunk(c("d <- diagnostics_display(result)",
            "d[[\"統計量\"]] <- signif(d[[\"統計量\"]], 4)",
            "d[[\"p 值\"]] <- signif(d[[\"p 值\"]], 4)",
            "knitr::kable(d)")),
    "### Mardia 與 Box's M 的合併解讀", "", "`r guide$diag_cases`", "",
    chunk(c("dn <- diagnostic_notes(result$diagnostics)",
            "print(knitr::kable(dn$table))",
            "cat(\"\\n\", paste(dn$notes, collapse = \"\\n\\n\"), \"\\n\")"),
          "results: asis"),
    "### 對照實際驗證結果", "",
    chunk("cat(diagnostic_outcome(result)$full, sep = \"\\n\\n\")", "results: asis"),
    "## 模型比較與選擇", "",
    "`r guide$metrics`必須完成全部重複的方法才參加比較。", "",
    chunk("knitr::kable(metrics_display(result))"),
    chunk("cat(selection_text(result), sep = \"\\n\\n\")", "results: asis"),
    "### 每一次重複的結果", "", "`r guide$repeats`", "",
    chunk("knitr::kable(repeat_table(result))"),
    "### 各類別表現（`r result$settings$repeats` 次合併）", "", "`r guide$class_recall`", "",
    chunk(c("mt <- misclass_text(result)",
            "if (length(c(mt$perfect, mt$recall)))",
            "  cat(\"本資料：\", mt$perfect, mt$recall, \"\\n\", sep = \"\")"), "results: asis"),
    chunk("knitr::kable(class_recall_display(result))"),
    "### 混淆矩陣（`r result$settings$repeats` 次合併，列百分比）", "",
    "`r paste(c(guide$confusion, misclass_text(result)$errors), collapse = '')`", "",
    chunk(c("for (m in unique(result$confusion_matrix$model)) {",
            "  cat(\"\\n**\", toupper(m), \"**\\n\\n\", sep = \"\")",
            "  print(knitr::kable(confusion_percent(result, m)))",
            "}"), "results: asis"),
    "## 分組描述統計", "", "`r guide$data_summary`", "",
    chunk("knitr::kable(data_summary_display(result), digits = 4)"),
    "## 圖", "", "`r guide$figures`", "",
    chunk(c("figs <- save_figures(result, data, \"figures_qmd\")",
            "fn <- figure_notes(result, figs)",
            "for (i in seq_len(nrow(fn)))",
            "  cat(\"### \", fn$title[i], \"\\n\\n![](figures_qmd/\", fn$file[i], \")\\n\\n\",",
            "      fn$text[i], \"\\n\\n\", sep = \"\")"), "results: asis"),
    "## 模型警告與錯誤", "",
    chunk(c("if (nrow(result$issues)) knitr::kable(result$issues) else",
            "  cat(\"沒有記錄到警告或錯誤。\")"), "results: asis"),
    "## 執行環境", "",
    chunk("sessionInfo()"),
    chunk("cat(conclusion_lines(result), sep = \"\\n\")", "results: asis"))
  path <- file.path(out, "report.qmd")
  writeLines(enc2utf8(qmd), path, useBytes = TRUE)
  invisible(path)
}
