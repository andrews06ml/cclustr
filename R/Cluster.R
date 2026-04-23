#' Perform clustering on multiple imputed datasets
#'
#' @description
#' Applies a clustering algorithm to each completed dataset in a standardized
#' imputation list (as produced by \code{\link{as_mild_list}}). Supports
#' hierarchical clustering, k-means, PAM, fuzzy c-means, model-based clustering,
#' k-modes (for categorical data), and k-prototypes (for mixed data).
#' Accepts a single value of \code{k} or a range, optionally scales the data
#' prior to clustering, and supports multiple distance metrics including
#' Gower, Jaccard, Simple Matching Coefficient, and user-defined distances.
#'
#' @param imp_list A named list of \code{data.frame} objects, as returned by
#'   \code{\link{as_mild_list}}. All datasets must have identical dimensions
#'   and column names.
#'
#' @param method A character string specifying the clustering algorithm.
#'   For hierarchical clustering, accepted values are \code{"ward.D"},
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"centroid"}, \code{"median"}, and \code{"mcquitty"}.
#'   Additional options are:
#'   \itemize{
#'     \item \code{"kmeans"}: k-means clustering for numeric data.
#'     \item \code{"pam"}: partitioning around medoids; supports mixed data
#'           when used with \code{distance = "gower"}.
#'     \item \code{"fuzzy"}: fuzzy c-means clustering.
#'     \item \code{"mclust"}: model-based clustering via Gaussian mixtures.
#'     \item \code{"kmodes"}: clustering for purely categorical data.
#'     \item \code{"kprototypes"}: clustering for mixed data (numeric and
#'           categorical).
#'   }
#'   Default is \code{"ward.D2"}.
#'
#' @param k A single integer or an integer vector specifying the number of
#'   clusters. If a vector is provided, clustering is performed for each
#'   value of \code{k}.
#'
#' @param scale_data Character string specifying the scaling strategy for
#'   continuous numeric variables (binary variables coded as 0/1 are
#'   automatically excluded from scaling). Options are:
#'   \itemize{
#'     \item \code{"global"}: (default) variables are scaled using pooled
#'           mean and standard deviation computed across all imputations
#'           without stacking them in memory (combinatorial variance formula).
#'           Ensures comparability of distance scales between datasets and
#'           is recommended for consensus clustering.
#'     \item \code{"within"}: each imputed dataset is scaled independently
#'           using its own mean and standard deviation. Useful for
#'           exploratory analyses but may introduce inconsistencies in
#'           distance scales across datasets.
#'     \item \code{"none"}: no scaling is applied.
#'   }
#'   Scaling is automatically disabled for \code{distance = "gower"},
#'   \code{"jaccard"}, or \code{"simple_matching"}, and for
#'   \code{method = "kmodes"} and \code{"kprototypes"}, as these handle
#'   variable scales internally.
#'
#' @param distance A character string specifying the distance metric used for
#'   hierarchical clustering and \code{method = "pam"}. Options are:
#'   \itemize{
#'     \item \code{"euclidean"}: (default) standard Euclidean distance.
#'           Suitable for continuous numeric data after scaling.
#'     \item \code{"manhattan"}: sum of absolute differences. Suitable for
#'           continuous numeric data; more robust to outliers than Euclidean.
#'     \item \code{"gower"}: Gower distance via \code{cluster::daisy()}.
#'           Supports mixed data types (numeric, categorical, binary).
#'           When used with hierarchical clustering, numeric variables are
#'           normalized using global min/max computed across all imputations,
#'           ensuring cross-imputation comparability. Requires package
#'           \pkg{cluster}.
#'     \item \code{"jaccard"}: Jaccard dissimilarity. Requires all variables
#'           to be binary (0/1 or two-level factor). Requires package
#'           \pkg{proxy}.
#'     \item \code{"simple_matching"}: Simple Matching Coefficient (SMC).
#'           Requires all variables to be binary. Requires package
#'           \pkg{proxy}.
#'     \item \code{"custom"}: user-defined distance function supplied via
#'           \code{dist_fun}. Must return an object of class \code{dist}.
#'   }
#'   Ignored for \code{method = "kmeans"}, \code{"fuzzy"}, \code{"mclust"},
#'   \code{"kmodes"}, and \code{"kprototypes"}, which compute distances
#'   internally.
#'
#'   For \code{method = "pam"}, all distance options are fully supported
#'   as PAM operates on a precomputed dissimilarity matrix.
#'
#' @param dist_fun A function that takes a \code{data.frame} and returns a
#'   \code{dist} object. Required when \code{distance = "custom"};
#'   ignored otherwise.
#'
#' @param dist_args A named list of additional arguments passed to the distance
#'   function. For \code{distance = "gower"}, arguments are forwarded to
#'   \code{cluster::daisy()}; for \code{distance = "jaccard"} or
#'   \code{"simple_matching"}, to \code{proxy::dist()}; for
#'   \code{distance = "euclidean"} or \code{"manhattan"}, to
#'   \code{stats::dist()}; also for \code{distance = "custom"}. Default is \code{list()}.
#'
#' @param ... Additional arguments passed to the underlying clustering
#'   function: \code{stats::kmeans}, \code{e1071::cmeans},
#'   \code{mclust::Mclust}, \code{klaR::kmodes}, or
#'   \code{clustMixType::kproto}. Not used for \code{method = "pam"}, which
#'   operates on a precomputed dissimilarity matrix.
#'
#' @return
#' \itemize{
#'   \item If \code{k} is a single value: a named list of integer vectors,
#'         one per imputed dataset, containing the cluster assignment for
#'         each observation.
#'   \item If \code{k} is a vector: a named list of lists, where each
#'         top-level element corresponds to one value of \code{k} (named
#'         \code{"k2"}, \code{"k3"}, etc.), and each inner element is a
#'         named list of cluster assignment vectors, one per imputed dataset.
#' }
#'
#' @details
#' \strong{Scaling strategy and cross-imputation comparability:}
#'
#' For methods that rely on Euclidean or Manhattan distances
#' (\code{"kmeans"}, \code{"fuzzy"}, \code{"mclust"}, hierarchical with
#' \code{distance = "euclidean"} or \code{"manhattan"}), \code{scale_data =
#' "global"} is strongly recommended. It computes the pooled mean and
#' standard deviation across all imputations using the combinatorial variance
#' formula — without creating a stacked copy of the data in memory — and
#' applies the same transformation to every imputation. This ensures that
#' distances are on the same scale across datasets, which is essential for
#' a meaningful co-assignment matrix.
#'
#' Binary variables coded as 0/1 numeric are automatically excluded from
#' scaling because standardizing them loses their binary semantics.
#'
#' \strong{Gower distance and global normalization:}
#'
#' Gower distance normalizes each numeric variable by its range. When
#' computed independently per imputation, different imputations may produce
#' slightly different ranges, leading to incomparable distance scales. When
#' \code{distance = "gower"} is used with hierarchical clustering or PAM,
#' the function pre-normalizes numeric variables using global min/max
#' computed across all imputations (without stacking), so that Gower
#' distances are on a consistent scale across datasets. Categorical and
#' factor variables are left untouched and handled normally by
#' \code{cluster::daisy()}. Additional arguments to \code{cluster::daisy()}
#' can be passed via \code{dist_args}.
#'
#' \strong{Hierarchical clustering optimization:}
#'
#' When \code{method} is hierarchical and \code{k} is a vector, the
#' distance matrix and dendrogram are computed only once per imputed dataset
#' and then cut at each requested \code{k}, which avoids redundant
#' computation.
#'
#' \strong{Data type compatibility:}
#' \itemize{
#'   \item Numeric-only: hierarchical, \code{"kmeans"}, \code{"fuzzy"},
#'         \code{"mclust"}.
#'   \item Categorical-only: \code{"kmodes"}.
#'   \item Mixed (numeric + categorical): \code{"pam"} with
#'         \code{distance = "gower"} or \code{"kprototypes"}.
#'   \item Binary-only: hierarchical with \code{distance = "jaccard"} or
#'         \code{"simple_matching"}.
#' }
#'
#' Character variables are internally converted to factors.
#'
#' @examples
#' # ------------------------------------------------------------
#' # Example 1: Numeric data, hierarchical clustering
#' # ------------------------------------------------------------
#' set.seed(123)
#' imp_list <- replicate(3, {
#'   data.frame(x = rnorm(20), y = rnorm(20))
#' }, simplify = FALSE)
#'
#' # Single k, global scaling, euclidean distance (default)
#' res <- cluster_imputations(imp_list, method = "ward.D2", k = 2)
#' str(res)
#'
#' # Range of k, manhattan distance
#' res_multi <- cluster_imputations(imp_list, method = "ward.D2",
#'                                  k = 2:4, distance = "manhattan")
#' str(res_multi)
#'
#' # ------------------------------------------------------------
#' # Example 2: Mixed data, PAM with Gower distance
#' # ------------------------------------------------------------
#' imp_mixed <- replicate(3, {
#'   data.frame(
#'     x    = rnorm(20),
#'     grup = factor(sample(letters[1:3], 20, replace = TRUE))
#'   )
#' }, simplify = FALSE)
#'
#' res_pam <- cluster_imputations(imp_mixed, method = "pam",
#'                                k = 2, distance = "gower")
#' str(res_pam)
#'
#' # ------------------------------------------------------------
#' # Example 3: Binary data, Jaccard distance
#' # ------------------------------------------------------------
#' imp_bin <- replicate(3, {
#'   data.frame(
#'     a = sample(0:1, 20, replace = TRUE),
#'     b = sample(0:1, 20, replace = TRUE),
#'     c = sample(0:1, 20, replace = TRUE)
#'   )
#' }, simplify = FALSE)
#'
#' \donttest{
#' if (requireNamespace("proxy", quietly = TRUE)) {
#'   res_jac <- cluster_imputations(imp_bin, method = "ward.D2",
#'                                  k = 2, distance = "jaccard")
#'   str(res_jac)
#' }
#'
#' # ------------------------------------------------------------
#' # Example 4: Full pipeline with mice
#' # ------------------------------------------------------------
#' if (requireNamespace("mice", quietly = TRUE)) {
#'
#'   df_pre       <- mice::nhanes
#'   df_pre$age   <- factor(df_pre$age)
#'   df_pre$hyp   <- factor(df_pre$hyp)
#'
#'   imp  <- mice::mice(df_pre, m = 3, printFlag = FALSE)
#'   mild <- as_mild_list(imp)
#'
#'   # PAM + Gower with global range normalization
#'   parts <- cluster_imputations(mild, method = "pam",
#'                                k = 2:4, distance = "gower")
#'   str(parts)
#' }
#' }
#'
#' @seealso \code{\link{as_mild_list}}, \code{\link{consensus_clustering}}
#'
#' @export
cluster_imputations <- function(imp_list,
                                method = "ward.D2",
                                k,
                                scale_data = c("global", "within", "none"),
                                distance = c("euclidean", "manhattan", "gower",
                                             "jaccard", "simple_matching", "custom"),
                                dist_fun = NULL,
                                dist_args = list(), # For distance
                                ...) {    # For clustering function

  scale_data <- match.arg(scale_data)
  distance   <- match.arg(distance)
  cluster_args <- list(...)

  # --------------------------------------------------------
  # Supported methods
  # --------------------------------------------------------
  hierarchical_methods <- c(
    "ward.D", "ward.D2", "single", "complete",
    "average", "centroid", "median", "mcquitty"
  )

  vector_methods <- c("kmeans", "fuzzy", "mclust")
  categorical_methods <- c("kmodes")
  mixed_methods <- c("pam", "kprototypes")

  if (!(method %in% c(hierarchical_methods, vector_methods,
                      categorical_methods, mixed_methods))) {
    stop("Unsupported clustering method.")
  }

  # --------------------------------------------------------
  # Detect variable types
  # --------------------------------------------------------

  imp_list <- lapply(imp_list, function(df) {
    df[] <- lapply(df, function(col) {
      if (is.character(col)) as.factor(col) else col
    })
    df
  })

  .detect_types <- function(df) {
    list(
      numeric = sapply(df, is.numeric),
      factor  = sapply(df, is.factor),
      binary  = sapply(df, function(x) length(unique(stats::na.omit(x))) <= 2)
    )
  }

  types <- .detect_types(imp_list[[1]])
  has_factors <- any(types$factor)
  all_factors <- all(types$factor)
  all_binary  <- all(types$binary)

  # --------------------------------------------------------
  # Validations
  # --------------------------------------------------------
  if (method %in% vector_methods && has_factors) {
    stop("The dataset contains categorical variables but method = '", method, "' requires numeric data only. ")
  }

  if (distance %in% c("jaccard", "simple_matching") && !all_binary) {
    stop(
      "distance = '", distance, "' requires all variables to be binary ",
      "(numeric 0/1 or two-level factor)."
    )
  }

  if (method %in% hierarchical_methods) {
    if (has_factors && distance %in% c("euclidean", "manhattan")) {
      stop("Categorical data requires distance = 'gower', 'jaccard', or 'custom'.")
    }
  }

  if (all_factors && method == "kprototypes") {
    warning("All variables are categorical. Switching method to 'kmodes'.")
    method <- "kmodes"
  }

  if (method == "kmodes" && !all_factors) {
    stop("kmodes requires categorical data.")
  }

  if (method == "pam" && has_factors && distance != "gower") {
    warning("Mixed data detected. Consider distance = 'gower'.")
  }

  if (any(k < 2)) stop("k must be >= 2.")
  n_obs <- nrow(imp_list[[1]])
  if (any(k >= n_obs)) stop("k must be < number of observations.")

  if (method == "kprototypes" && distance != "euclidean") {
    warning("kprototypes handles distance internally. The 'distance' parameter is ignored.")
  }

  if (method %in% vector_methods && distance != "euclidean") {
    warning(paste0("The 'distance' parameter is ignored for method = '", method, "'."))
  }

  if (distance == "custom" && is.null(dist_fun)) {
    stop("distance = 'custom' requires a function supplied via dist_fun.")
  }

  # --------------------------------------------------------
  # Disable scaling when not applicable
  # --------------------------------------------------------
  if (distance %in% c("gower", "jaccard", "simple_matching")) {
    if (scale_data != "none") {
      if (distance == "gower") {
        message(
          "The 'scale_data' parameter is ignored for distance = 'gower'. ",
          "Numeric variables are instead normalized using global min/max ",
          "computed across all imputations to ensure cross-imputation ",
          "comparability before computing Gower distances."
        )
      } else {
        message(
          "The 'scale_data' parameter is ignored for distance = '", distance, "'. ",
          "This metric handles normalization internally."
        )
      }
      scale_data <- "none"
    }
  }

  if (method %in% categorical_methods && scale_data != "none") {
    warning("Scaling disabled for categorical-only method.")
    scale_data <- "none"
  }

  # kprototypes: handles numeric/categorical balance internally
  if (method == "kprototypes" && scale_data != "none") {
    warning(
      "Scaling disabled for method = 'kprototypes'. ",
      "kprototypes handles variable scales internally. ",
      "Variables with larger numeric ranges may still dominate; ",
      "consider pre-scaling manually if needed."
    )
    scale_data <- "none"
  }

  # --------------------------------------------------------
  # Scaling (only numeric)
  # --------------------------------------------------------

  num_cols <- types$numeric & !types$binary

  if (!any(num_cols)) {
    database <- imp_list
  } else if (scale_data == "none") {
    database <- imp_list
  } else if (scale_data == "within") {

    database <- lapply(imp_list, function(df) {
      df[, num_cols] <- scale(df[, num_cols, drop = FALSE])
      df
    })

  } else if (scale_data == "global") {

    # Pooled mean and SD without stacking (combinatorial variance formula)
    # This avoids creating a large stacked object in memory.
    imp_ns    <- sapply(imp_list, nrow)
    imp_means <- lapply(imp_list,
                        function(df) colMeans(df[, num_cols, drop = FALSE]))
    imp_vars  <- lapply(imp_list,
                        function(df) apply(df[, num_cols, drop = FALSE], 2,
                                           stats::var))

    # Grand (pooled) mean: weighted average of per-imputation means
    grand_mean <- Reduce("+",
                         mapply(function(mu, ni) mu * ni,
                                imp_means, imp_ns, SIMPLIFY = FALSE)
    ) / sum(imp_ns)

    # Pooled variance: between- + within-imputation components
    pool_var <- Reduce("+",
                       mapply(function(v, mu, ni) {
                         (ni - 1) * v + ni * (mu - grand_mean) ^ 2
                       }, imp_vars, imp_means, imp_ns, SIMPLIFY = FALSE)
    ) / (sum(imp_ns) - 1)

    grand_sd          <- sqrt(pool_var)
    grand_sd[grand_sd == 0] <- 1   # avoid division by zero

    database <- lapply(imp_list, function(df) {
      df[, num_cols] <- sweep(df[, num_cols, drop = FALSE], 2, grand_mean, "-")
      df[, num_cols] <- sweep(df[, num_cols, drop = FALSE], 2, grand_sd,   "/")
      df
    })
  } else {
    database <- imp_list
  }

  # --------------------------------------------------------
  # Global range normalization for Gower distance
  # Pre-normalize numeric variables using pooled min/max so
  # that Gower distances are comparable across imputations.
  # daisy() will then see numeric variables already in [0, 1]
  # and compute its own range as ~1, leaving the scale intact.
  # Factor/categorical columns are not touched.
  # --------------------------------------------------------

  if (distance == "gower") {

    num_cols_gwr <- types$numeric & !types$binary

    if (any(num_cols_gwr)) {

      # Global min and max without stacking: element-wise pmin/pmax
      global_min <- do.call(pmin,
                            lapply(database,
                                   function(df) apply(df[, num_cols_gwr, drop = FALSE],
                                                      2, min, na.rm = TRUE)))
      global_max <- do.call(pmax,
                            lapply(database,
                                   function(df) apply(df[, num_cols_gwr, drop = FALSE],
                                                      2, max, na.rm = TRUE)))

      # Apply global range normalization to each dataset
      database <- lapply(database, function(df) {
        for (col in names(which(num_cols_gwr))) {
          rng <- global_max[col] - global_min[col]
          if (rng == 0) rng <- 1
          df[[col]] <- (df[[col]] - global_min[col]) / rng
        }
        df
      })
    }
  }

  # --------------------------------------------------------
  # Distance computation
  # --------------------------------------------------------
  .compute_distance <- function(df) {

    if (distance == "gower") {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Package 'cluster' is required for Gower distance.")
      }
      return(do.call(cluster::daisy, c(list(x = df,
                                            metric = "gower"), dist_args)
      ))
    }

    if (distance == "jaccard") {
      if (!all_binary) stop("Jaccard requires binary variables.")
      if (!requireNamespace("proxy", quietly = TRUE)) {
        stop("Package 'proxy' is required.")
      }
      return(do.call(proxy::dist,
                     c(list(x = df,
                            method = "Jaccard"), dist_args)
      ))
    }

    if (distance == "simple_matching") {
      if (!all_binary) stop("simple matching requires binary variables.")
      if (!requireNamespace("proxy", quietly = TRUE)) {
        stop("Package 'proxy' is required.")
      }
      return(do.call(proxy::dist,
                     c(list(x = df,
                            method = "simple matching"), dist_args)
      ))
    }

    if (distance == "custom") {
      result <- do.call(dist_fun, c(list(df), dist_args))
      if (!inherits(result, "dist")) {
        stop("dist_fun must return an object of class 'dist'.")
      }
      return(result)
    }

    # euclidean / manhattan
    do.call(stats::dist,
            c(list(x = df,
                   method = distance), dist_args))
  }

  # --------------------------------------------------------
  # Clustering
  # --------------------------------------------------------
  .run_clustering <- function(df, k_i) {

    if (method %in% hierarchical_methods) {
      d <- .compute_distance(df)
      hc <- stats::hclust(d, method = method)
      return(stats::cutree(hc, k = k_i))
    }

    if (method == "kmeans") {
      return(do.call(stats::kmeans,
                     c(list(x= df, centers = k_i), cluster_args))$cluster)
    }

    if (method == "pam") {
      d <- .compute_distance(df)
      return(cluster::pam(d, k_i, diss = TRUE)$clustering)
    }

    if (method == "fuzzy") {
      return(do.call(e1071::cmeans,
                     c(list(x = df, centers = k_i), cluster_args))$cluster)
    }

    if (method == "mclust") {
      return(do.call(mclust::Mclust,
                     c(list(data = df, G = k_i), cluster_args))$classification)
    }

    if (method == "kmodes") {
      return(do.call(klaR::kmodes,
                     c(list(data = df, modes = k_i), cluster_args))$cluster)
    }

    if (method == "kprototypes") {
      return(do.call(clustMixType::kproto,
                     c(list(x = df, k = k_i), cluster_args))$cluster)
    }
  }

  # --------------------------------------------------------
  # Single k
  # --------------------------------------------------------
  if (length(k) == 1) {
    partitions <- lapply(database, function(df) {
      .run_clustering(df, k)
    })

    names(partitions) <- names(imp_list)
    return(partitions)
  }

  # --------------------------------------------------------
  # Multiple k
  # --------------------------------------------------------

  if (length(k) > 1) {

    if (method %in% hierarchical_methods) {
      # Calcular dendrograma UNA SOLA VEZ por imputación
      hc_list <- lapply(database, function(df) {
        d <- .compute_distance(df)
        stats::hclust(d, method = method)
      })

      partitions_by_k <- lapply(k, function(k_i) {
        parts <- lapply(hc_list, function(hc) stats::cutree(hc, k = k_i))
        names(parts) <- names(imp_list)
        parts
      })

    } else {

      partitions_by_k <- lapply(k, function(k_i) {
        partitions <- lapply(database, function(df) {
          .run_clustering(df, k_i)
        })
        names(partitions) <- names(imp_list)
        partitions
      })
    }

    names(partitions_by_k) <- paste0("k", k)
    return(partitions_by_k)
  }
}
