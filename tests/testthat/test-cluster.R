# tests/testthat/test-cluster.R
# Tests para cluster_imputations() y consensus_clustering()

# -------------------------------------------------------
# Reusable Helpers
# -------------------------------------------------------
make_num_imp <- function(n = 20, m = 3, seed = 1) {
  set.seed(seed)
  replicate(m, data.frame(x = rnorm(n), y = rnorm(n)), simplify = FALSE)
}

make_mixed_imp <- function(n = 20, m = 3, seed = 1) {
  set.seed(seed)
  replicate(m, data.frame(
    x    = rnorm(n),
    group = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  ), simplify = FALSE)
}

make_bin_imp <- function(n = 20, m = 3, seed = 1) {
  set.seed(seed)
  replicate(m, data.frame(
    a = sample(0:1, n, replace = TRUE),
    b = sample(0:1, n, replace = TRUE),
    c = sample(0:1, n, replace = TRUE)
  ), simplify = FALSE)
}

# ============================================================
# cluster_imputations — output structure (single k)
# ============================================================

test_that("single k: returns a list of integer vectors", {
  set.seed(1)
  res <- cluster_imputations(make_num_imp(), method = "ward.D2", k = 2)
  expect_type(res, "list")
  expect_true(all(vapply(res, is.integer, logical(1))))
})

test_that("single k: number of partitions equals number of imputations", {
  imp <- make_num_imp(m = 4)
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 2)
  expect_length(res, 4)
})

test_that("single k: each partition has as many elements as observations", {
  imp <- make_num_imp(n = 25, m = 3)
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 3)
  expect_true(all(vapply(res, length, integer(1)) == 25))
})

test_that("single k: cluster values are in {1, ..., k}", {
  imp <- make_num_imp()
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 3)
  expect_true(all(vapply(res, function(z) all(z %in% 1:3), logical(1))))
})

test_that("single k: output names match input list names", {
  imp <- make_num_imp(m = 3)
  names(imp) <- c("imp1", "imp2", "imp3")
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 2)
  expect_equal(names(res), names(imp))
})

# ============================================================
# cluster_imputations — output structure (multiple k)
# ============================================================

test_that("multiple k: returns a list of lists named k2, k3, ...", {
  imp <- make_num_imp()
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 2:4)
  expect_named(res, c("k2", "k3", "k4"))
  expect_true(all(vapply(res, is.list, logical(1))))
})

test_that("multiple k: each sublist has as many partitions as imputations", {
  imp <- make_num_imp(m = 3)
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 2:3)
  expect_true(all(vapply(res, length, integer(1)) == 3))
})

test_that("multiple k: cluster values are correct for each k", {
  imp <- make_num_imp()
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 2:3)
  expect_true(all(vapply(res$k2, function(z) all(z %in% 1:2), logical(1))))
  expect_true(all(vapply(res$k3, function(z) all(z %in% 1:3), logical(1))))
})

# ============================================================
# cluster_imputations — hierarchical methods
# ============================================================

test_that("alternative hierarchical methods work", {
  imp <- make_num_imp()
  for (mth in c("ward.D", "single", "complete", "average")) {
    set.seed(1)
    res <- cluster_imputations(imp, method = mth, k = 2)
    expect_length(res, 3)
  }
})

# ============================================================
# cluster_imputations — kmeans
# ============================================================

test_that("kmeans: single k produces valid partitions", {
  imp <- make_num_imp()
  set.seed(1)
  res <- cluster_imputations(imp, method = "kmeans", k = 2)
  expect_length(res, 3)
  expect_true(all(vapply(res, function(z) all(z %in% 1:2), logical(1))))
})

test_that("kmeans: fails with categorical data", {
  imp <- make_mixed_imp()
  expect_error(
    cluster_imputations(imp, method = "kmeans", k = 2),
    "categorical"
  )
})

# ============================================================
# cluster_imputations — PAM with Gower distance (mixed data)
# ============================================================

test_that("PAM + Gower: works with mixed data", {
  imp <- make_mixed_imp()
  set.seed(1)
  expect_warning(
    res <- cluster_imputations(imp, method = "pam", k = 2, distance = "gower"),
    "ignored"
  )
  expect_length(res, 3)
  expect_true(all(vapply(res, function(z) all(z %in% 1:2), logical(1))))
})

test_that("PAM + Gower: multiple k returns correct structure", {
  imp <- make_mixed_imp()
  set.seed(1)
  expect_warning(
    res <- cluster_imputations(imp, method = "pam", k = 2:3, distance = "gower"),
    "ignored"
  )
  expect_named(res, c("k2", "k3"))
})

# ============================================================
# cluster_imputations — input validations
# ============================================================

test_that("fails with k = 1", {
  imp <- make_num_imp()
  expect_error(cluster_imputations(imp, method = "ward.D2", k = 1), "k must be")
})

test_that("fails with k >= n_obs", {
  imp <- make_num_imp(n = 10)
  expect_error(cluster_imputations(imp, method = "ward.D2", k = 10), "k must be")
})

test_that("fails with unsupported method", {
  imp <- make_num_imp()
  expect_error(
    cluster_imputations(imp, method = "non_existent_method", k = 2),
    "Unsupported"
  )
})

test_that("fails with distance = 'custom' without dist_fun", {
  imp <- make_num_imp()
  expect_error(
    cluster_imputations(imp, method = "ward.D2", k = 2, distance = "custom"),
    "dist_fun"
  )
})

test_that("fails with categorical data and euclidean distance in hierarchical", {
  imp <- make_mixed_imp()
  expect_error(
    cluster_imputations(imp, method = "ward.D2", k = 2, distance = "euclidean"),
    "Categorical data"
  )
})

test_that("fails with jaccard on non-binary data", {
  imp <- make_num_imp()
  expect_error(
    cluster_imputations(imp, method = "ward.D2", k = 2, distance = "jaccard"),
    "binary"
  )
})

# ============================================================
# cluster_imputations — global vs within scaling
# ============================================================

test_that("scale_data = 'global' does not change the number of partitions", {
  imp <- make_num_imp()
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 2, scale_data = "global")
  expect_length(res, 3)
})

test_that("scale_data = 'within' produces output equivalent in structure", {
  imp <- make_num_imp()
  set.seed(1)
  res <- cluster_imputations(imp, method = "ward.D2", k = 2, scale_data = "within")
  expect_length(res, 3)
  expect_true(all(vapply(res, function(z) all(z %in% 1:2), logical(1))))
})

# ============================================================
# consensus_clustering — output structure (single k)
# ============================================================

test_that("single k: returns list with mandatory elements", {
  set.seed(1)
  parts <- list(
    imp1 = sample(1:2, 20, replace = TRUE),
    imp2 = sample(1:2, 20, replace = TRUE),
    imp3 = sample(1:2, 20, replace = TRUE)
  )
  res <- consensus_clustering(parts, k = 2)
  expect_named(res, c("consensus_method", "k", "consensus",
                      "coassignment", "hclust", "weights"),
               ignore.order = TRUE)
})

test_that("single k: consensus has length n", {
  set.seed(1)
  parts <- list(
    imp1 = sample(1:2, 20, replace = TRUE),
    imp2 = sample(1:2, 20, replace = TRUE)
  )
  res <- consensus_clustering(parts, k = 2)
  expect_length(res$consensus, 20)
})

test_that("single k: coassignment is an n x n matrix with values in [0, 1]", {
  set.seed(1)
  parts <- list(
    imp1 = sample(1:2, 20, replace = TRUE),
    imp2 = sample(1:2, 20, replace = TRUE)
  )
  res <- consensus_clustering(parts, k = 2)
  expect_true(is.matrix(res$coassignment))
  expect_equal(dim(res$coassignment), c(20, 20))
  expect_true(all(res$coassignment >= 0 & res$coassignment <= 1))
})

test_that("single k: coassignment diagonal is 1 (each obs. with itself)", {
  set.seed(1)
  parts <- list(
    imp1 = sample(1:3, 15, replace = TRUE),
    imp2 = sample(1:3, 15, replace = TRUE)
  )
  res <- consensus_clustering(parts, k = 3)
  expect_true(all(diag(res$coassignment) == 1))
})

test_that("single k: coassignment is symmetric", {
  set.seed(1)
  parts <- list(
    imp1 = sample(1:2, 15, replace = TRUE),
    imp2 = sample(1:2, 15, replace = TRUE)
  )
  res <- consensus_clustering(parts, k = 2)
  expect_equal(res$coassignment, t(res$coassignment))
})

test_that("single k: weights sum to 1 (classic method)", {
  set.seed(1)
  parts <- list(
    imp1 = sample(1:2, 15, replace = TRUE),
    imp2 = sample(1:2, 15, replace = TRUE),
    imp3 = sample(1:2, 15, replace = TRUE)
  )
  res <- consensus_clustering(parts, k = 2, consensus_method = "classic")
  expect_equal(sum(res$weights), 1, tolerance = 1e-10)
})

test_that("single k: hclust is an object of class hclust", {
  set.seed(1)
  parts <- list(
    imp1 = sample(1:2, 15, replace = TRUE),
    imp2 = sample(1:2, 15, replace = TRUE)
  )
  res <- consensus_clustering(parts, k = 2)
  expect_s3_class(res$hclust, "hclust")
})

# ============================================================
# consensus_clustering — output structure (multiple k)
# ============================================================

test_that("multiple k: returns a list named k2, k3, ...", {
  imp <- make_num_imp()
  set.seed(1)
  parts <- cluster_imputations(imp, method = "ward.D2", k = 2:4)
  res   <- consensus_clustering(parts)
  expect_named(res, c("k2", "k3", "k4"))
})

test_that("multiple k: each element has the single k structure", {
  imp <- make_num_imp()
  set.seed(1)
  parts <- cluster_imputations(imp, method = "ward.D2", k = 2:3)
  res   <- consensus_clustering(parts)
  for (ki in res) {
    expect_true(all(c("consensus", "coassignment", "hclust", "weights") %in% names(ki)))
  }
})

test_that("multiple k: the k field of each element is correct", {
  imp <- make_num_imp()
  set.seed(1)
  parts <- cluster_imputations(imp, method = "ward.D2", k = 2:4)
  res   <- consensus_clustering(parts)
  expect_equal(res$k2$k, 2L)
  expect_equal(res$k3$k, 3L)
  expect_equal(res$k4$k, 4L)
})

# ============================================================
# consensus_clustering — input validations
# ============================================================

test_that("fails with empty list", {
  expect_error(consensus_clustering(list(), k = 2), "non-empty")
})

test_that("fails with single k and no k argument", {
  parts <- list(imp1 = c(1, 2, 1), imp2 = c(1, 1, 2))
  expect_error(consensus_clustering(parts), "k")
})

test_that("fails with k < 2", {
  parts <- list(imp1 = c(1, 2, 1), imp2 = c(1, 1, 2))
  expect_error(consensus_clustering(parts, k = 1), "k must be")
})

test_that("fails with unsupported cluster_method", {
  parts <- list(imp1 = c(1, 2, 1), imp2 = c(1, 1, 2))
  expect_error(
    consensus_clustering(parts, k = 2, cluster_method = "bad_method"),
    "cluster_method"
  )
})

test_that("fails if partitions have different lengths", {
  parts <- list(imp1 = c(1, 2, 1), imp2 = c(1, 2))
  expect_error(consensus_clustering(parts, k = 2), "same length")
})
