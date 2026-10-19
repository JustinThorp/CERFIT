# Unit tests for the CERFIT package (testthat, 3rd edition).
#
# Sections
#   0. Helpers and shared fixtures
#   1. Input validation
#   2. Fitted-object structure and type detection
#   3. Working response (residualization)
#   4. Resampling and out-of-bag (OOB) storage
#   5. Tree-growing controls
#   6. predict(): shape, naming, and internal consistency
#   7. MinDepth()
#   8. print()
#   9. Internal helpers and the C++ split statistic
#  10. Propensity-score paths
#
# Forests use very few trees so the whole file runs in seconds on CRAN.
# Tests that call twang/CBPS are slow and are skipped on CRAN.

# ---------------------------------------------------------------------------
# 0. Helpers and shared fixtures
# ---------------------------------------------------------------------------

# CERFIT() and predict.CERFIT() print progress with cat()/print().
# Capture that output so the test log stays readable.
quiet <- function(expr) {
  capture.output(res <- expr)
  res
}

warts_df <- as.data.frame(warts) # plain data.frame keeps row names on subsetting
n_warts <- nrow(warts_df)
f_warts <- Result_of_Treatment ~ sex + age + Number_of_Warts + Area + Time + Type | treatment
warts_covars <- c("sex", "age", "Number_of_Warts", "Area", "Time", "Type")

f_educ <- Y ~ SAT_MATH + HSGPA + AGE + GENDER + URM | A

# CERFIT() detects survival outcomes by the bare symbol `Surv` on the formula's
# left-hand side, so bind it locally rather than attaching the survival package.
Surv <- survival::Surv

# Simulated right-censored survival data with a treatment-by-x2 interaction
sim_surv <- function(n = 200, seed = 4) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1, 0.5)
  trt <- stats::rbinom(n, 1, 0.5)
  t_event <- stats::rexp(n, exp(0.5 * x1 + 0.7 * trt * x2))
  t_cens <- stats::rexp(n, 0.3)
  data.frame(
    time = pmin(t_event, t_cens),
    status = as.integer(t_event <= t_cens),
    x1 = x1, x2 = x2, trt = trt, id = seq_len(n)
  )
}
surv_df <- sim_surv()

# Simulated continuous-treatment data (> 15 unique treatment values)
sim_cont <- function(n = 300, seed = 3) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  tr <- round(stats::runif(n, 0, 10) + 0.5 * x1, 2)
  y <- x1 + 0.5 * tr * (x2 > 0) + stats::rnorm(n)
  data.frame(y = y, x1 = x1, x2 = x2, tr = tr)
}
cont_df <- sim_cont()

# Shared fits (grown once, reused across tests)
set.seed(101)
fit_rct <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 10, method = "RCT", mtry = 2))

set.seed(102)
fit_multi <- quiet(CERFIT(f_educ, data = educational, ntrees = 5, method = "RCT"))

set.seed(103)
fit_surv <- quiet(CERFIT(Surv(time, status) ~ x1 + x2 | trt,
                         data = surv_df, ntrees = 5, method = "RCT"))

set.seed(104)
# Continuous-treatment fit for predict() tests. minbucket = 30 gives every leaf
# >= 30 in-bag rows, so the 10-fold cv.glmnet() in predictTree() has >= 3 rows
# per fold and does not warn "grouped=FALSE enforced". This relies on
# oob = FALSE (the default); with oob = TRUE leaves store only out-of-bag rows,
# which minbucket does not bound.
fit_cont <- quiet(CERFIT(y ~ x1 + x2 | tr, data = cont_df, ntrees = 5,
                         method = "RCT", minbucket = 30))

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("CERFIT() errors when the treatment has a single level", {
  one_arm <- warts_df
  one_arm$treatment <- 1
  expect_error(
    quiet(CERFIT(f_warts, data = one_arm, ntrees = 2, method = "RCT")),
    "Only one treatment"
  )
})

test_that("CERFIT() rejects malformed survival specifications", {
  expect_error(
    quiet(CERFIT(Surv(time) ~ x1 | trt, data = surv_df, ntrees = 2, method = "RCT")),
    "Surv\\(time, event\\) or Surv\\(start, stop, event\\)"
  )
  expect_error(
    quiet(CERFIT(Surv(time, status) ~ x1 | trt, data = surv_df, ntrees = 2,
                 method = "RCT", surv_id = 1)),
    "surv_id must be a character string"
  )
  # Valid counting-process data (start < stop on every row), so these checks
  # fail only because of surv_id and survival::Surv() emits no warnings
  cp_df <- surv_df
  cp_df$start <- 0  # simulated times are exponential draws, so time > 0

  # Counting-process format requires surv_id ...
  expect_error(
    quiet(CERFIT(Surv(start, time, status) ~ x1 | trt, data = cp_df,
                 ntrees = 2, method = "RCT")),
    "surv_id must be provided"
  )
  # ... and surv_id must name a column in data
  expect_error(
    quiet(CERFIT(Surv(start, time, status) ~ x1 | trt, data = cp_df,
                 ntrees = 2, method = "RCT", surv_id = "not_a_column")),
    "surv_id must name a column"
  )
})

test_that("CERFIT() rejects invalid argument values", {
  expect_error(
    quiet(CERFIT(f_warts, data = warts_df, ntrees = 2, method = "RCT", sampleMethod = "bogus")),
    "should be one of"
  )
  # HI is only defined for continuous treatments
  expect_error(
    quiet(CERFIT(f_warts, data = warts_df, ntrees = 2, method = "observational", PropForm = "HI")),
    "Not able to estimate"
  )
})

test_that("predict() rejects invalid type and incomplete newdata", {
  expect_error(quiet(predict(fit_rct, type = "bogus")), "should be one of")
  expect_error(quiet(predict(fit_rct, newdata = warts_df[, c("sex", "age")])),
               "not found")
})

# ---------------------------------------------------------------------------
# 2. Fitted-object structure and type detection
# ---------------------------------------------------------------------------

test_that("CERFIT() returns a well-formed CERFIT object", {
  expect_s3_class(fit_rct, "CERFIT")
  expect_named(fit_rct, c("randFor", "trt.type", "response.type", "useRes", "data"))
  expect_length(fit_rct$randFor, 10)
  expect_true(fit_rct$useRes)
  expect_equal(nrow(fit_rct$data), n_warts)

  for (b in fit_rct$randFor) {
    expect_named(b, c("tree", "cases"))
    expect_s3_class(b$tree, "constparty")
    expect_false(is.unsorted(b$cases))
    expect_false(anyDuplicated(b$cases) > 0)
    expect_true(all(b$cases >= 1 & b$cases <= n_warts))
  }
})

test_that("response type is detected correctly", {
  expect_equal(fit_rct$response.type, "continuous")
  expect_equal(fit_surv$response.type, "survival")

  bin_df <- warts_df
  bin_df$Result_of_Treatment <- factor(bin_df$Result_of_Treatment)
  set.seed(1)
  fit_bin <- quiet(CERFIT(f_warts, data = bin_df, ntrees = 2, method = "RCT", mtry = 2))
  expect_equal(fit_bin$response.type, "binary")
})

test_that("treatment type is detected correctly", {
  expect_equal(fit_rct$trt.type, "binary")
  expect_equal(fit_multi$trt.type, "multiple")

  ord_df <- educational
  ord_df$A <- factor(ord_df$A, ordered = TRUE)
  set.seed(1)
  fit_ord <- quiet(CERFIT(f_educ, data = ord_df, ntrees = 2, method = "RCT"))
  expect_equal(fit_ord$trt.type, "ordered")

  set.seed(1)
  fit_cont <- quiet(CERFIT(y ~ x1 + x2 | tr, data = cont_df, ntrees = 2, method = "RCT"))
  expect_equal(fit_cont$trt.type, "continuous")
})

test_that("RCT fits use unit propensity scores and unit weights", {
  expect_true(all(fit_rct$data$prop == 1))
  expect_true(all(fit_rct$data$iptw == 1))
})

# ---------------------------------------------------------------------------
# 3. Working response (residualization)
# ---------------------------------------------------------------------------

test_that("continuous response is replaced by OLS residuals when useRes = TRUE", {
  ols <- stats::lm(Result_of_Treatment ~ sex + age + Number_of_Warts + Area + Time + Type,
                   data = warts_df)
  expect_equal(unname(fit_rct$data$Result_of_Treatment), unname(stats::resid(ols)))
  expect_equal(fit_rct$data$yo, warts_df$Result_of_Treatment) # original kept in yo
})

test_that("binary response is replaced by logistic-regression residuals", {
  bin_df <- warts_df
  bin_df$Result_of_Treatment <- factor(bin_df$Result_of_Treatment)
  set.seed(1)
  fit_bin <- quiet(CERFIT(f_warts, data = bin_df, ntrees = 2, method = "RCT", mtry = 2))
  lr <- stats::glm(Result_of_Treatment ~ sex + age + Number_of_Warts + Area + Time + Type,
                   data = bin_df, family = stats::binomial)
  expected <- (as.numeric(bin_df$Result_of_Treatment) - 1) - stats::fitted(lr)
  expect_equal(unname(fit_bin$data$Result_of_Treatment), unname(expected))
})

test_that("useRes = FALSE leaves the response untouched", {
  set.seed(1)
  fit_raw <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 2, method = "RCT", mtry = 2,
                          useRes = FALSE))
  expect_false(fit_raw$useRes)
  expect_equal(fit_raw$data$Result_of_Treatment, warts_df$Result_of_Treatment)
})

test_that("survival response uses Cox deviance residuals by default", {
  cox <- survival::coxph(survival::Surv(time, status) ~ x1 + x2, data = surv_df)
  expect_true(".CERFIT_y" %in% names(fit_surv$data))
  expect_equal(unname(fit_surv$data$.CERFIT_y),
               unname(stats::resid(cox, type = "deviance")))
  expect_equal(fit_surv$data$yo_time, surv_df$time)
  expect_equal(fit_surv$data$yo_delta, surv_df$status)
  expect_true(all(fit_surv$data$surv_resid_type == "deviance"))
  expect_false(any(fit_surv$data$surv_is_counting))
})

test_that("surv_resid = 'martingale' uses Cox martingale residuals", {
  set.seed(1)
  fit_mart <- quiet(CERFIT(Surv(time, status) ~ x1 + x2 | trt, data = surv_df,
                           ntrees = 2, method = "RCT", surv_resid = "martingale"))
  cox <- survival::coxph(survival::Surv(time, status) ~ x1 + x2, data = surv_df)
  expect_equal(unname(fit_mart$data$.CERFIT_y),
               unname(stats::resid(cox, type = "martingale")))
  expect_true(all(fit_mart$data$surv_resid_type == "martingale"))
})

test_that("counting-process survival data use row-level Cox residuals with id", {
  cp_df <- survival::survSplit(Surv(time, status) ~ ., data = surv_df,
                               cut = c(0.5, 1), episode = "ep")
  set.seed(1)
  fit_cp <- quiet(CERFIT(Surv(tstart, time, status) ~ x1 + x2 | trt, data = cp_df,
                         ntrees = 2, method = "RCT", surv_id = "id"))
  cox <- survival::coxph(survival::Surv(tstart, time, status) ~ x1 + x2,
                         data = cp_df, id = id)
  expect_equal(fit_cp$response.type, "survival")
  expect_equal(unname(fit_cp$data$.CERFIT_y), unname(stats::resid(cox, type = "deviance")))
  expect_true(all(fit_cp$data$surv_is_counting))
  expect_equal(fit_cp$data$yo_tstart, cp_df$tstart)
  expect_equal(fit_cp$data$yo_tstop, cp_df$time)
  expect_equal(fit_cp$data$yo_id, cp_df$id)
})

# ---------------------------------------------------------------------------
# 4. Resampling and out-of-bag (OOB) storage
# ---------------------------------------------------------------------------

test_that("sampleMethod controls which rows each tree is grown on", {
  set.seed(1)
  fit_all <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 3, method = "RCT", mtry = 2,
                          sampleMethod = "allData"))
  for (b in fit_all$randFor) expect_identical(b$cases, seq_len(n_warts))

  set.seed(1)
  fit_sub <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 3, method = "RCT", mtry = 2,
                          sampleMethod = "subsample"))
  for (b in fit_sub$randFor) expect_length(b$cases, round(0.632 * n_warts))
})

test_that("with oob = FALSE each tree stores its in-bag bootstrap sample", {
  # Bootstrap draws n rows with replacement, so the stored sample has n rows
  for (b in fit_rct$randFor) {
    expect_equal(nrow(partykit::data_party(b$tree)), n_warts)
  }
})

test_that("with oob = TRUE each tree stores exactly its out-of-bag rows", {
  set.seed(2)
  fit_oob <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 3, method = "RCT", mtry = 2,
                          oob = TRUE))
  for (b in fit_oob$randFor) {
    stored <- partykit::data_party(b$tree)
    expect_equal(nrow(stored), n_warts - length(b$cases))
    expect_identical(as.integer(rownames(stored)), setdiff(seq_len(n_warts), b$cases))
  }
})

test_that("fits are reproducible under set.seed()", {
  set.seed(9)
  a <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 3, method = "RCT", mtry = 2))
  set.seed(9)
  b <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 3, method = "RCT", mtry = 2))
  expect_identical(quiet(predict(a, type = "ITE")), quiet(predict(b, type = "ITE")))
})

# ---------------------------------------------------------------------------
# 5. Tree-growing controls
# ---------------------------------------------------------------------------

test_that("maxdepth bounds the number of terminal nodes", {
  for (d in 1:2) {
    set.seed(1)
    fit_d <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 5, method = "RCT", mtry = 2,
                          maxdepth = d))
    n_term <- vapply(fit_d$randFor,
                     function(b) length(partykit::nodeids(b$tree, terminal = TRUE)),
                     numeric(1))
    expect_true(all(n_term <= 2^d))
  }
})

test_that("every in-bag terminal node contains at least minbucket rows", {
  set.seed(1)
  fit_mb <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 5, method = "RCT", mtry = 2,
                         minbucket = 15))
  for (b in fit_mb$randFor) {
    node_sizes <- table(stats::fitted(b$tree)[["(fitted)"]])
    expect_true(all(node_sizes >= 15))
  }
})

test_that("continuous-treatment fixture leaves meet the 30-row cv.glmnet floor", {
  for (b in fit_cont$randFor) {
    node_sizes <- table(stats::fitted(b$tree)[["(fitted)"]])
    expect_true(all(node_sizes >= 30))
  }
})

# ---------------------------------------------------------------------------
# 6. predict(): shape, naming, and internal consistency
# ---------------------------------------------------------------------------

test_that("type = 'response' returns one column per treatment level", {
  resp <- quiet(predict(fit_rct, type = "response"))
  expect_true(is.matrix(resp))
  expect_equal(dim(resp), c(n_warts, 2))
  expect_equal(colnames(resp), c("y=0", "y=1"))
  expect_true(all(is.finite(resp)))
})

test_that("type = 'ITE' equals the response contrast against the first level", {
  resp <- quiet(predict(fit_rct, type = "response"))
  ite <- quiet(predict(fit_rct, type = "ITE"))
  expect_equal(dim(ite), c(n_warts, 1))
  expect_equal(colnames(ite), "y1-y0")
  expect_equal(unname(ite[, 1]), unname(resp[, 2] - resp[, 1]))
})

test_that("multi-arm ITEs are contrasts against the reference (lowest) level", {
  resp <- quiet(predict(fit_multi, type = "response"))
  ite <- quiet(predict(fit_multi, type = "ITE"))
  expect_equal(colnames(resp), c("y=0", "y=1", "y=2"))
  expect_equal(colnames(ite), c("y1-y0", "y2-y0"))
  expect_equal(unname(ite), unname(resp[, 2:3] - resp[, 1]))
})

test_that("type = 'opT' returns the arg-max treatment and its predicted response", {
  resp <- quiet(predict(fit_multi, type = "response"))
  opt <- quiet(predict(fit_multi, type = "opT"))
  expect_equal(colnames(opt), c("opTreat", "opResponse"))
  expect_true(all(opt[, "opTreat"] %in% 1:3))
  expect_equal(unname(opt[, "opResponse"]), unname(apply(resp, 1, max)))
  expect_equal(unname(opt[, "opTreat"]), max.col(resp, ties.method = "first"))
})

test_that("predict() on the training data matches predict() with newdata = NULL", {
  expect_equal(quiet(predict(fit_rct)), quiet(predict(fit_rct, newdata = warts_df)))
})

test_that("predictions for newdata are row-wise and ignore the response column", {
  idx <- c(10, 3, 150, 42)
  nd <- warts_df[idx, ]
  p_nd <- quiet(predict(fit_rct, newdata = nd))
  p_all <- quiet(predict(fit_rct))

  expect_equal(nrow(p_nd), length(idx))
  expect_equal(unname(p_nd), unname(p_all[idx, ]))

  # Row order is preserved
  p_rev <- quiet(predict(fit_rct, newdata = nd[rev(seq_along(idx)), ]))
  expect_equal(unname(p_rev), unname(p_nd[rev(seq_along(idx)), ]))

  # Outcome values in newdata do not affect predictions
  nd_y <- nd
  nd_y$Result_of_Treatment <- 999
  expect_equal(quiet(predict(fit_rct, newdata = nd_y)), p_nd)
})

test_that("a single-row newdata returns a one-row matrix", {
  ite1 <- quiet(predict(fit_rct, newdata = warts_df[5, , drop = FALSE], type = "ITE"))
  expect_equal(dim(ite1), c(1, 1))
})

test_that("prediction = 'by iter' ITE converges to the overall ITE (binary treatment)", {
  ite <- quiet(predict(fit_rct, type = "ITE"))
  cum <- quiet(predict(fit_rct, type = "ITE", prediction = "by iter"))
  ntree <- length(fit_rct$randFor)

  expect_type(cum, "list")
  expect_length(cum, 1)
  expect_equal(dim(cum[[1]]), c(n_warts, ntree))
  expect_equal(unname(cum[[1]][, ntree]), unname(ite[, 1]))
})

test_that("prediction = 'by iter' ITE converges to the overall ITE (multi-arm treatment)", {
  ite <- quiet(predict(fit_multi, type = "ITE"))
  cum <- quiet(predict(fit_multi, type = "ITE", prediction = "by iter"))
  ntree <- length(fit_multi$randFor)

  expect_type(cum, "list")
  expect_length(cum, ncol(ite))
  for (j in seq_along(cum)) {
    expect_equal(dim(cum[[j]]), c(nrow(ite), ntree))
    expect_equal(unname(cum[[j]][, ntree]), unname(ite[, j]))
  }
})

test_that("prediction = 'by iter' response converges to the overall response", {
  resp <- quiet(predict(fit_rct, type = "response"))
  cum <- quiet(predict(fit_rct, type = "response", prediction = "by iter"))

  expect_length(cum, ncol(resp))
  for (j in seq_along(cum)) {
    expect_equal(unname(cum[[j]][, ncol(cum[[j]])]), unname(resp[, j]))
  }
})

test_that("prediction = 'by iter' output is labeled by treatment level", {
  resp_cum <- quiet(predict(fit_multi, type = "response", prediction = "by iter"))
  ite_cum <- quiet(predict(fit_multi, type = "ITE", prediction = "by iter"))
  lv <- sub("^y=", "", colnames(quiet(predict(fit_multi, type = "response"))))

  expect_named(resp_cum, paste0("ycum", lv))
  expect_named(ite_cum, paste0("ycum", lv[-1], "-ycum", lv[1]))

  # Each named contrast must hold that arm's cumulative response minus the
  # reference arm's, so a label/value mismatch can't pass on names alone
  for (k in lv[-1]) {
    expect_equal(
      ite_cum[[paste0("ycum", k, "-ycum", lv[1])]],
      resp_cum[[paste0("ycum", k)]] - resp_cum[[paste0("ycum", lv[1])]],
      label = paste0("ITE contrast ycum", k, "-ycum", lv[1])
    )
  }
})

test_that("continuous treatment predicts at the default 21-point quantile grid", {
  set.seed(11); expect_no_warning(resp <- quiet(predict(fit_cont)))
  set.seed(11); expect_no_warning(ite  <- quiet(predict(fit_cont, type = "ITE")))
  expect_equal(dim(resp), c(nrow(cont_df), 21))  # 0%, 5%, ..., 100% quantiles
  expect_equal(ncol(ite), 20)
  expect_equal(unname(ite), unname(resp[, -1] - resp[, 1]))
})

test_that("continuous treatments reject a gridval of 10 or fewer values", {

  # Short grids error for every prediction type
  expect_error(predict(fit_cont, gridval = c(2, 5, 8)),
               "gridval must contain more than 10 values")
  expect_error(predict(fit_cont, gridval = c(2, 5, 8), type = "ITE"),
               "gridval must contain more than 10 values")

  # Boundary: exactly 10 values still errors
  expect_error(predict(fit_cont, gridval = seq(2, 8, length.out = 10)),
               "gridval must contain more than 10 values")

  # 11 values is accepted and predictions are made at exactly that grid.
  # cv.glmnet uses random folds, so reset the seed before each call to make
  # resp and ite come from the same fits
  grid <- seq(2, 8, length.out = 11)
  set.seed(11); expect_no_warning(resp <- quiet(predict(fit_cont, gridval = grid)))
  set.seed(11); expect_no_warning(ite  <- quiet(predict(fit_cont, gridval = grid, type = "ITE")))
  expect_equal(dim(resp), c(nrow(cont_df), length(grid)))
  expect_equal(colnames(resp), paste0("y=", grid))
  expect_equal(colnames(ite), paste0("y", grid[-1], "-y", grid[1]))
  expect_equal(unname(ite), unname(resp[, -1] - resp[, 1]))

  # Shifting the grid must change predictions; the categorical branch
  # ignores gridval values and would return identical output
  set.seed(11); resp_shift <- quiet(predict(fit_cont, gridval = grid + 0.5))
  expect_false(isTRUE(all.equal(unname(resp), unname(resp_shift))))
})

test_that("survival predictions are on the residual scale and internally consistent", {
  resp <- quiet(predict(fit_surv, type = "response"))
  ite <- quiet(predict(fit_surv, type = "ITE"))
  expect_equal(dim(ite), c(nrow(surv_df), 1))
  expect_equal(unname(ite[, 1]), unname(resp[, 2] - resp[, 1]))
  # Deviance residuals are bounded well inside this range
  expect_true(all(is.finite(resp)))
  expect_true(all(abs(resp) < 10))
  # Predicting on the original survival data frame gives the same ITEs
  expect_equal(quiet(predict(fit_surv, newdata = surv_df, type = "ITE")), ite)
})

test_that("opT with by iter errors", {
  expect_error(predict(fit_multi, type = "opT", prediction = "by iter"),
               "not supported")
})

# ---------------------------------------------------------------------------
# 7. MinDepth()
# ---------------------------------------------------------------------------

test_that("MinDepth() returns a named, non-negative importance vector", {
  md <- MinDepth(fit_rct)
  expect_type(md, "double")
  expect_named(md, warts_covars)
  expect_true(all(is.finite(md) & md >= 0))
})

test_that("MinDepth() matches its definition for single-split trees", {
  # With maxdepth = 1 every tree has one split at depth 0; the other p - 1
  # variables are absent and receive treeDepth + 1 = 2. So the per-tree sum is
  # 2 * (p - 1), and averaging over trees preserves that total.
  set.seed(1)
  fit_stump <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 6, method = "RCT", mtry = 2,
                            maxdepth = 1))
  md <- MinDepth(fit_stump)
  expect_equal(sum(md), 2 * (length(warts_covars) - 1))
  expect_true(all(md >= 0 & md <= 2))
})

# ---------------------------------------------------------------------------
# 8. print()
# ---------------------------------------------------------------------------

test_that("print() reports tree count, treatment type, and response type", {
  expect_output(print(fit_rct), "Trees: 10")
  expect_output(print(fit_rct), "Treatment Type: binary")
  expect_output(print(fit_rct), "Response Type: continuous")
})

# ---------------------------------------------------------------------------
# 9. Internal helpers and the C++ split statistic
# ---------------------------------------------------------------------------

test_that("truncquant() winsorizes at the q and 1 - q quantiles", {
  x <- 1:100
  tq <- CERFIT:::truncquant(x, q = 0.9)
  qs <- stats::quantile(x, c(0.1, 0.9), names = FALSE)
  expect_equal(tq, pmin(pmax(x, qs[1]), qs[2]))
})

test_that("findCutpts() returns midpoints that respect minbucket", {
  expect_equal(CERFIT:::findCutpts(1:8, minbucket = 2), c(2.5, 3.5, 4.5, 5.5, 6.5))
  expect_null(CERFIT:::findCutpts(rep(1, 10), minbucket = 3))
})

test_that("ordinalize() codes factor levels by mean response", {
  out <- CERFIT:::ordinalize(factor(c("a", "b", "c", "a")), c(1, 10, 5, 3))
  expect_equal(out$x, c(2, 10, 5, 2))
  expect_equal(as.numeric(out$cutToLvl), c(2, 10, 5))
  # Numeric inputs pass through unchanged
  num <- CERFIT:::ordinalize(c(3, 1, 2), 1:3)
  expect_equal(num$x, c(3, 1, 2))
  expect_null(num$cutToLvl)
})

test_that("find_split() (RCT) returns the squared t-statistic of the interaction", {
  set.seed(11)
  n <- 400
  x <- stats::runif(n)
  trt <- stats::rbinom(n, 1, 0.5)
  y <- 2 * trt * (x <= 0.5) + stats::rnorm(n, sd = 0.5)

  res <- CERFIT:::find_split(y = y, x = x, trt = trt, cutpts = c(0.25, 0.5, 0.75),
                             method = "RCT", propensity = rep(1, n), minbucket = 7,
                             response_type = "continuous")
  expect_equal(res$cutoff, 0.5) # recovers the true change point

  z <- as.numeric(x <= 0.5)
  t_int <- summary(stats::lm(y ~ z * trt))$coefficients["z:trt", "t value"]
  expect_equal(res$stat, t_int^2)
})

test_that("find_split() (observational) returns the HC0 sandwich Wald statistic", {
  set.seed(12)
  n <- 400
  x <- stats::runif(n)
  trt <- stats::rbinom(n, 1, 0.5)
  e <- stats::runif(n, 0.2, 0.8)
  y <- 2 * trt * (x <= 0.5) + stats::rnorm(n, sd = 0.5)

  res <- CERFIT:::find_split(y = y, x = x, trt = trt, cutpts = 0.5,
                             method = "observational", propensity = e, minbucket = 7,
                             response_type = "continuous")

  z <- as.numeric(x <= 0.5)
  X <- cbind(1, z, trt, z * trt, e)
  XtX_inv <- solve(crossprod(X))
  beta <- XtX_inv %*% crossprod(X, y)
  r <- as.vector(y - X %*% beta)
  V <- XtX_inv %*% crossprod(X * r) %*% XtX_inv
  expect_equal(res$stat, beta[4]^2 / V[4, 4])
})

test_that("find_split() returns NA when no cut point satisfies the size rules", {
  set.seed(13)
  n <- 200
  x <- stats::runif(n)
  trt <- stats::rbinom(n, 1, 0.5)
  y <- stats::rnorm(n)
  res <- CERFIT:::find_split(y = y, x = x, trt = trt, cutpts = min(x) + 1e-8,
                             method = "RCT", propensity = rep(1, n), minbucket = 7,
                             response_type = "continuous")
  expect_true(is.na(res$stat))
  expect_true(is.na(res$cutoff))
})

# ---------------------------------------------------------------------------
# 10. Propensity-score paths
# ---------------------------------------------------------------------------

test_that("observational HI generalized propensity scores are valid for continuous treatment", {
  set.seed(1)
  fit_hi <- quiet(CERFIT(y ~ x1 + x2 | tr, data = cont_df, ntrees = 2,
                         method = "observational", PropForm = "HI"))
  expect_equal(fit_hi$trt.type, "continuous")
  expect_true(all(is.finite(fit_hi$data$prop) & fit_hi$data$prop > 0))
  expect_true(all(is.finite(fit_hi$data$iptw) & fit_hi$data$iptw > 0))
})

test_that("observational CBPS fit works for a three-arm treatment", {
  skip_on_cran()
  skip_if_not_installed("CBPS")
  set.seed(1)
  fit_cbps <- quiet(CERFIT(f_educ, data = educational, ntrees = 3,
                           method = "observational", PropForm = "CBPS"))
  expect_equal(fit_cbps$trt.type, "multiple")
  expect_true(all(is.finite(fit_cbps$data$iptw) & fit_cbps$data$iptw > 0))
  expect_equal(dim(quiet(predict(fit_cbps, type = "ITE"))), c(nrow(educational), 2))
})

test_that("observational GBM (twang) fit works for a binary treatment", {
  skip_on_cran() # twang::ps() with 10,000 trees is slow
  skip_if_not_installed("twang")
  set.seed(1)
  fit_gbm <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 3, method = "observational",
                          PropForm = "GBM", mtry = 2))
  prop <- as.numeric(unlist(fit_gbm$data$prop))
  expect_true(all(prop > 0 & prop < 1))
  expect_true(all(is.finite(fit_gbm$data$iptw) & fit_gbm$data$iptw > 0))
})

test_that("observational randomForest propensity scores are valid probabilities", {
  set.seed(1)
  fit_rf <- quiet(CERFIT(f_warts, data = warts_df, ntrees = 3, method = "observational",
                         PropForm = "randomForest", mtry = 2))
  expect_length(fit_rf$data$prop, n_warts)
  expect_true(all(fit_rf$data$prop >= 0 & fit_rf$data$prop <= 1))
  expect_true(all(is.finite(fit_rf$data$iptw) & fit_rf$data$iptw > 0))
  # IPTW must vary; regression guard for the truncquant(Iptw[[1]]) bug
  expect_gt(length(unique(fit_rf$data$iptw)), 1)
})
