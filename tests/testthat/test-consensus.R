# tests/testthat/test-consensus.R
# Tests for validate_clustering() and choose_best_clustering()

# -------------------------------------------------------
# Reusable Helpers
# -------------------------------------------------------

# Perfectly separated partitions (two distinct groups)
make_perfect_parts <- function(n = 20, m = 3) {
  labels <- rep(1:2, each = n / 2)
  replicate(m, labels, simplify = FALSE) |>
    setNames(paste0("imp", seq_len(m)))
}

# Random partitions (pure noise)
make_random_parts <- function(n = 20, m = 3, k = 2, seed = 1) {
  set.seed(seed)
  replicate(m, sample(seq_len(k), n, replace = TRUE), simplify = FALSE) |>
    setNames(paste0("imp", seq_len(m)))
}

# Synthetic validation table for choose_best_clustering
make_val_table <- function() {
  data.frame(
    k                            = 2:4,
    pac                          = c(0.05, 0.20, 0.35),
    silhouette_mean              = c(0.70, 0.50, 0.30),
    ari_mean_between_imputations = c(0.90, 0.70, 0.50),
    ari_consensus_mean           = c(0.88, 0.65, 0.45),
    calinski_harabasz_mean       = c(200,  120,   80),
    davies_bouldin_mean          = c(0.30, 0.60,  0.90),
    dunn_index                   = c(0.50, 0.30,  0.20)
  )
}

# Synthetic consensus results for choose_best_clustering
make_cons_results <- function(n = 20) {
  set.seed(1)
  make_one <- function(k) list(
    k            = k,
    consensus    = sample(seq_len(k), n, replace = TRUE),
    coassignment = matrix(runif(n * n), n, n)
  )
  list(k2 = make_one(2), k3 = make_one(3), k4 = make_one(4))
}

# ============================================================
# validate_clustering — output structure
# ============================================================

test_that("returns a data.frame with one row per evaluated k", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  expect_s3_class(val, "data.frame")
  expect_equal(nrow(val), 1L)
})

test_that("contains all expected columns", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  expected_cols <- c("k", "pac", "silhouette_mean",
                     "ari_mean_between_imputations", "ari_consensus_mean",
                     "calinski_harabasz_mean", "davies_bouldin_mean", "dunn_index")
  expect_true(all(expected_cols %in% names(val)))
})

test_that("multiple k: returns one row for each evaluated k", {
  set.seed(1)
  imp   <- replicate(3, data.frame(x = rnorm(20), y = rnorm(20)), simplify = FALSE)
  parts <- cluster_imputations(imp, method = "ward.D2", k = 2:4)
  cons  <- consensus_clustering(parts)
  val   <- validate_clustering(parts, cons)
  expect_equal(nrow(val), 3L)
  expect_equal(sort(val$k), 2:4)
})

# ============================================================
# validate_clustering — metric ranges and properties
# ============================================================

test_that("pac is in [0, 1]", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  expect_gte(val$pac, 0)
  expect_lte(val$pac, 1)
})

test_that("silhouette_mean is in [-1, 1]", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  if (!is.na(val$silhouette_mean)) {
    expect_gte(val$silhouette_mean, -1)
    expect_lte(val$silhouette_mean, 1)
  }
})

test_that("ari_mean_between_imputations is in [-1, 1]", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  if (!is.na(val$ari_mean_between_imputations)) {
    expect_gte(val$ari_mean_between_imputations, -1)
    expect_lte(val$ari_mean_between_imputations,  1)
  }
})

test_that("ari_consensus_mean is in [-1, 1]", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  if (!is.na(val$ari_consensus_mean)) {
    expect_gte(val$ari_consensus_mean, -1)
    expect_lte(val$ari_consensus_mean,  1)
  }
})

test_that("calinski_harabasz_mean is positive or NA", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  if (!is.na(val$calinski_harabasz_mean)) {
    expect_gt(val$calinski_harabasz_mean, 0)
  }
})

test_that("davies_bouldin_mean is non-negative or NA", {
  parts <- make_random_parts(k = 2)
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  if (!is.na(val$davies_bouldin_mean)) {
    expect_gte(val$davies_bouldin_mean, 0)
  }
})

# ============================================================
# validate_clustering — quality signals with ideal data
# ============================================================

test_that("perfect partitions produce low pac", {
  parts <- make_perfect_parts()
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  # With perfect clusters PAC should be very low
  expect_lt(val$pac, 0.1)
})

test_that("perfect partitions produce high silhouette", {
  parts <- make_perfect_parts()
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  if (!is.na(val$silhouette_mean)) {
    expect_gt(val$silhouette_mean, 0.5)
  }
})

test_that("perfect partitions produce ARI between imputations = 1", {
  parts <- make_perfect_parts()
  cons  <- consensus_clustering(parts, k = 2)
  val   <- validate_clustering(parts, cons)
  expect_equal(val$ari_mean_between_imputations, 1, tolerance = 1e-10)
})

# ============================================================
# validate_clustering — input validations
# ============================================================

test_that("fails if partitions is empty", {
  cons <- consensus_clustering(
    list(imp1 = c(1, 2, 1), imp2 = c(1, 2, 2)), k = 2
  )
  expect_error(validate_clustering(list(), cons), "non-empty")
})

test_that("fails if consensus_results is not a list", {
  parts <- make_random_parts()
  expect_error(validate_clustering(parts, "not_a_list"), "list")
})

# ============================================================
# choose_best_clustering — output structure
# ============================================================

test_that("returns list with all expected elements", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(choose_best_clustering(val, cons))
  expected <- c("best_k", "best_consensus", "best_coassignment",
                "best_consensus_result", "scores_table", "weights", "tie_breaker")
  expect_true(all(expected %in% names(best)))
})

test_that("best_k is one of the evaluated k", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(choose_best_clustering(val, cons))
  expect_true(best$best_k %in% 2:4)
})

test_that("scores_table has score column", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(choose_best_clustering(val, cons))
  expect_true("score" %in% names(best$scores_table))
})

test_that("scores_table is sorted by k ascending", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(choose_best_clustering(val, cons))
  expect_equal(best$scores_table$k, sort(best$scores_table$k))
})

test_that("best_consensus has length n", {
  val  <- make_val_table()
  cons <- make_cons_results(n = 20)
  best <- suppressMessages(choose_best_clustering(val, cons))
  expect_length(best$best_consensus, 20)
})

test_that("best_coassignment is a square matrix", {
  val  <- make_val_table()
  cons <- make_cons_results(n = 20)
  best <- suppressMessages(choose_best_clustering(val, cons))
  expect_true(is.matrix(best$best_coassignment))
  expect_equal(nrow(best$best_coassignment), ncol(best$best_coassignment))
})

# ============================================================
# choose_best_clustering — selection logic
# ============================================================

test_that("selects k=2 when it is clearly better", {
  # k=2 dominates in all metrics of the synthetic table
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(choose_best_clustering(val, cons))
  expect_equal(best$best_k, 2L)
})

test_that("prefer_stability = TRUE uses stability weights", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(
    choose_best_clustering(val, cons, prefer_stability = TRUE)
  )
  expect_equal(best$weights[["pac"]], 2)
  expect_equal(best$weights[["ari_between"]], 2)
})

test_that("prefer_stability = FALSE uses compactness weights", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(
    choose_best_clustering(val, cons, prefer_stability = FALSE)
  )
  expect_equal(best$weights[["silhouette"]], 2)
  expect_equal(best$weights[["ch"]], 2)
})

test_that("prefer_stability = NULL uses equal weights", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(
    choose_best_clustering(val, cons, prefer_stability = NULL)
  )
  expect_true(all(best$weights == 1))
})

test_that("custom weights are respected", {
  val  <- make_val_table()
  cons <- make_cons_results()
  w    <- c(pac = 5, silhouette = 1, ari_between = 1,
            ari_consensus = 1, ch = 1, db = 1, dunn = 1)
  best <- suppressMessages(choose_best_clustering(val, cons, weights = w))
  expect_equal(best$weights[["pac"]], 5)
})

test_that("tie_breaker is recorded in the output", {
  val  <- make_val_table()
  cons <- make_cons_results()
  best <- suppressMessages(
    choose_best_clustering(val, cons, tie_breaker = "dunn")
  )
  expect_equal(best$tie_breaker, "dunn")
})

# ============================================================
# choose_best_clustering — input validations
# ============================================================

test_that("fails if columns are missing in validation_table", {
  val_bad <- data.frame(k = 2:3, pac = c(0.1, 0.2))
  cons    <- make_cons_results()
  expect_error(
    choose_best_clustering(val_bad, cons),
    "missing"
  )
})

test_that("fails with negative weights", {
  val  <- make_val_table()
  cons <- make_cons_results()
  w    <- c(pac = -1, silhouette = 1, ari_between = 1,
            ari_consensus = 1, ch = 1, db = 1, dunn = 1)
  expect_error(
    choose_best_clustering(val, cons, weights = w),
    ">= 0"
  )
})

test_that("fails with all weights at zero", {
  val  <- make_val_table()
  cons <- make_cons_results()
  w    <- c(pac = 0, silhouette = 0, ari_between = 0,
            ari_consensus = 0, ch = 0, db = 0, dunn = 0)
  expect_error(
    choose_best_clustering(val, cons, weights = w),
    "positive"
  )
})

test_that("fails with non-logical prefer_stability", {
  val  <- make_val_table()
  cons <- make_cons_results()
  expect_error(
    choose_best_clustering(val, cons, prefer_stability = "yes"),
    "prefer_stability"
  )
})
