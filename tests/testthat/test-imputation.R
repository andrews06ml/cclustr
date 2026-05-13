# tests/testthat/test-imputation.R
# Tests for as_mild_list()

# -------------------------------------------------------
# Helper: valid list of data.frames with no NAs
# -------------------------------------------------------
make_imp <- function(n = 15, m = 3, seed = 1) {
  set.seed(seed)
  replicate(m, data.frame(x = rnorm(n), y = rnorm(n)), simplify = FALSE)
}

# ============================================================
# Output structure
# ============================================================

test_that("returns a named list of data.frames", {
  result <- as_mild_list(make_imp())
  expect_type(result, "list")
  expect_true(all(vapply(result, is.data.frame, logical(1))))
  expect_false(is.null(names(result)))
})

test_that("number of elements matches the number of input imputations", {
  result <- as_mild_list(make_imp(m = 5))
  expect_length(result, 5)
})

test_that("all datasets have the same dimensions", {
  result <- as_mild_list(make_imp(n = 20, m = 4))
  nrows  <- vapply(result, nrow, integer(1))
  ncols  <- vapply(result, ncol, integer(1))
  expect_true(length(unique(nrows)) == 1)
  expect_true(length(unique(ncols)) == 1)
})

test_that("all datasets have identical column names", {
  result    <- as_mild_list(make_imp())
  col_names <- lapply(result, names)
  expect_true(all(vapply(col_names, identical, logical(1), col_names[[1]])))
})

test_that("datasets are base data.frames, not tibbles", {
  result <- as_mild_list(make_imp())
  expect_true(all(vapply(result,
                         function(df) identical(class(df), "data.frame"),
                         logical(1))))
})

# ============================================================
# Long-format input with .imp column
# ============================================================

test_that("accepts a long-format data.frame with .imp column", {
  set.seed(42)
  long_df <- data.frame(
    .imp = rep(1:3, each = 10),
    x    = rnorm(30),
    y    = rnorm(30)
  )
  result <- as_mild_list(long_df)
  expect_length(result, 3)
  expect_true(all(vapply(result, is.data.frame, logical(1))))
})

test_that("excludes rows with .imp == 0 (original incomplete data)", {
  set.seed(42)
  long_df <- data.frame(
    .imp = rep(0:2, each = 10),
    x    = rnorm(30),
    y    = rnorm(30)
  )
  result <- as_mild_list(long_df)
  expect_length(result, 2)
})

test_that("drops .imp and .id columns from the output", {
  set.seed(42)
  long_df <- data.frame(
    .imp = rep(1:2, each = 10),
    .id  = rep(1:10, 2),
    x    = rnorm(20),
    y    = rnorm(20)
  )
  result <- as_mild_list(long_df)
  expect_false(".imp" %in% names(result[[1]]))
  expect_false(".id"  %in% names(result[[1]]))
})

# ============================================================
# Data quality checks
# ============================================================

test_that("fails if any imputation contains NA", {
  imp <- make_imp()
  imp[[2]][3, 1] <- NA
  expect_error(as_mild_list(imp), "NA")
})

test_that("fails if any imputation contains NaN", {
  imp <- make_imp()
  imp[[1]][1, 1] <- NaN
  expect_error(as_mild_list(imp), "NA values")
})

test_that("fails if any imputation contains Inf", {
  imp <- make_imp()
  imp[[1]][2, 2] <- Inf
  expect_error(as_mild_list(imp), "Inf")
})

test_that("fails if a numeric column is constant across all rows", {
  imp <- make_imp()
  imp[[1]]$x <- 5   # constant column
  expect_error(as_mild_list(imp), "constant")
})

test_that("fails if imputations have different numbers of rows", {
  set.seed(1)
  imp <- list(
    data.frame(x = rnorm(10)),
    data.frame(x = rnorm(12))
  )
  expect_error(as_mild_list(imp), "rows")
})

test_that("fails if imputations have different numbers of columns", {
  set.seed(1)
  imp <- list(
    data.frame(x = rnorm(10), y = rnorm(10)),
    data.frame(x = rnorm(10))
  )
  expect_error(as_mild_list(imp), "columns")
})

test_that("fails if imputations have different column names", {
  set.seed(1)
  imp <- list(
    data.frame(x = rnorm(10), y = rnorm(10)),
    data.frame(x = rnorm(10), z = rnorm(10))
  )
  expect_error(as_mild_list(imp), "column names")
})

# ============================================================
# Unsupported input types
# ============================================================

test_that("fails with unsupported input types", {
  expect_error(as_mild_list(42),      "Unsupported")
  expect_error(as_mild_list("text"),  "Unsupported")
  expect_error(as_mild_list(TRUE),    "Unsupported")
})

test_that("fails if the list contains non-data.frame elements", {
  bad_list <- list(data.frame(x = 1:5), matrix(1:5))
  expect_error(as_mild_list(bad_list), "data.frames")
})

test_that("fails with a data.frame that has no .imp column", {
  df <- data.frame(x = 1:10, y = rnorm(10))
  expect_error(as_mild_list(df), "\\.imp")
})
