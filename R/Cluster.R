#' Perform clustering on multiple imputed datasets
#'
#' @description
#' Applies a clustering algorithm to each completed dataset in a standardized
#' imputation list (as produced by \code{\link{as_mild_list}}). Supports
#' hierarchical clustering, k-means, PAM, fuzzy c-means, model-based clustering,
#' k-modes (for categorical data), and k-prototypes (for mixed data).
#' Accepts a single value of \code{k} or a range, and optionally scales the
#' data prior to clustering.
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
#'           when used with \code{metric = "gower"}.
#'     \item \code{"fuzzy"}: fuzzy c-means clustering.
#'     \item \code{"mclust"}: model-based clustering via Gaussian mixtures.
#'     \item \code{"kmodes"}: clustering for purely categorical data.
#'     \item \code{"kprototypes"}: clustering for mixed data (numeric and categorical).
#'   }
#'   Default is \code{"ward.D2"}.
#'
#' @param k A single integer or an integer vector specifying the number of
#'   clusters. If a vector is provided, clustering is performed for each
#'   value of \code{k}.
#'
#' @param scale_data Character string specifying the scaling strategy for
#'   numeric variables. Options are:
#'   \itemize{
#'     \item \code{"global"}: (default) variables are scaled using global
#'           mean and standard deviation computed across all imputations.
#'           Ensures comparability between datasets and is recommended for
#'           consensus clustering.
#'     \item \code{"within"}: each imputed dataset is scaled independently
#'           using its own mean and standard deviation. This option may be
#'           useful for exploratory analyses where the goal is to preserve
#'           the internal structure of each imputation. However, it can
#'           introduce inconsistencies in distance scales across datasets,
#'           potentially affecting the stability and interpretation of
#'           consensus clustering results.
#'     \item \code{"none"}: no scaling is applied.
#'   }
#'
#' @param ... Additional arguments passed to the underlying clustering
#'   function (\code{kmeans}, \code{cluster::pam}, \code{e1071::cmeans},
#'   \code{mclust::Mclust}, \code{klaR::kmodes}, or
#'   \code{clustMixType::kproto}).
#'   For \code{method = "pam"}, passing \code{metric = "gower"} enables
#'   Gower distance for mixed data types.
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
#' When \code{method} is hierarchical and \code{k} is a vector, the
#' dendrogram is computed only once per imputed dataset and then cut at
#' each requested \code{k}, which improves efficiency considerably.
#'
#' For \code{method = "pam"} with \code{metric = "gower"}, Gower distance
#' is computed via \code{cluster::daisy()} and passed as a dissimilarity
#' matrix, making the function suitable for datasets with mixed variable
#' types (numeric and categorical).
#'
#' The choice of clustering method should be consistent with the data types:
#' \itemize{
#'   \item Numeric-only datasets: methods such as \code{"kmeans"},
#'         hierarchical clustering, \code{"fuzzy"}, and \code{"mclust"}.
#'   \item Categorical datasets: \code{"kmodes"}.
#'   \item Mixed datasets (numeric and categorical): \code{"pam"} with
#'         \code{metric = "gower"} or \code{"kprototypes"}.
#' }
#'
#' Character variables are internally converted to factors to ensure
#' compatibility with categorical and mixed-data clustering methods.
#'
#' @examples
#' # ------------------------------------------------------------
#' # Example 1: with numeric data
#' # ------------------------------------------------------------
#' set.seed(123)
#'
#' # Simulate 3 imputed datasets (list of data frames)
#' imp_list <- replicate(3, {
#'   data.frame(
#'     x = rnorm(10),
#'     y = rnorm(10)
#'   )
#' }, simplify = FALSE)
#'
#' # Hierarchical clustering
#' res <- cluster_imputations(imp_list, method = "ward.D2", k = 2)
#' str(res)
#'
#' # k-means clustering
#' res_kmeans <- cluster_imputations(imp_list, method = "kmeans", k = 2)
#' str(res_kmeans)
#'
#' # ------------------------------------------------------------
#' # Example 2: with mixed data (numeric + categorical)
#' # ------------------------------------------------------------
#' imp_list_mixed <- replicate(2, {
#'   data.frame(
#'     x = rnorm(10),
#'     y = factor(sample(letters[1:3], 10, replace = TRUE))
#'   )
#' }, simplify = FALSE)
#'
#' res_pam <- cluster_imputations(
#'   imp_list_mixed,
#'   method = "pam",
#'   k = 2,
#'   metric = "gower"
#' )
#' str(res_pam)
#'
#'\donttest{
#' # ------------------------------------------------------------
#' # Example 3: mixed data (numeric + categorical) with mice
#' # ------------------------------------------------------------
#'if (requireNamespace("mice", quietly = TRUE)) {
#'
#'  # Prepare original data
#'  df_pre <- mice::nhanes
#'  df_pre$age <- factor(df_pre$age)
#'  df_pre$hyp <- factor(df_pre$hyp)
#'
#'  imp   <- mice::mice(df_pre, m = 3, printFlag = FALSE)
#'  mild  <- as_mild_list(imp)
#'
#'  # Single k with PAM with Gower distance (mixed data)
#'  parts_gower <- cluster_imputations(mild, method = "pam", k = 3,
#'                                     metric = "gower")
#'  parts_gower
#'
#'  # Multiple k with Gower distance (mixed data)
#'  parts_multi <- cluster_imputations(mild, method = "pam", k = 2:4,
#'                                     metric = "gower")
#'  parts_multi
#'}
#'}
#'
#' @seealso \code{\link{as_mild_list}}, \code{\link{consensus_clustering}}
#'
#' @export
cluster_imputations <- function(imp_list,
                                method = "ward.D2",
                                k,
                                scale_data = c("global", "within", "none"),
                                ...) {

  scale_data <- match.arg(scale_data)
  # Capture user-specified arguments for downstream clustering functions.
  # These arguments are method-dependent and forwarded internally.
  extra_args <- list(...)


  # --------------------------------------------------------
  # Supported methods
  # --------------------------------------------------------
  hierarchical_methods <- c(
    "ward.D", "ward.D2", "single", "complete",
    "average", "centroid", "median", "mcquitty"
  )

  if (!(method %in% c(hierarchical_methods, "kmeans", "pam", "fuzzy", "mclust",
                      "kmodes", "kprototypes"))) {
    stop("Unsupported clustering method.")
  }

  # --------------------------------------------------------
  # Optional escalation - Scaling strategy
  # --------------------------------------------------------

  # Disable scaling for methods that do not require it
  if (method %in% c("kmodes") && scale_data != "none") {
    warning("Scaling disabled for method = 'kmodes'.")
    scale_data <- "none"
  }

    # Identify numeric columns
    num_cols <- sapply(imp_list[[1]], is.numeric)

    if (!any(num_cols) && scale_data != "none") {
      warning("No numeric variables found. Scaling skipped.")
      database <- imp_list

    } else if (scale_data == "none") {

      database <- imp_list

    } else if (scale_data == "within"){

      # Scale each imputation independently
      database <- lapply(imp_list, function(df) {
        df[, num_cols] <- scale(df[, num_cols])
        df
      })

    } else if (scale_data == "global") {

      # Stack all imputations to compute global parameters
      stacked <- do.call(rbind, imp_list)

      means <- colMeans(stacked[, num_cols, drop = FALSE])
      sds   <- apply(stacked[, num_cols, drop = FALSE], 2, sd)

      # Avoid division by zero
      sds[sds == 0] <- 1

      # Apply same scaling to all imputations
      database <- lapply(imp_list, function(df) {
        df[, num_cols] <- sweep(df[, num_cols, drop = FALSE], 2, means, "-")
        df[, num_cols] <- sweep(df[, num_cols, drop = FALSE], 2, sds, "/")
        df
      })
    }

  # --------------------------------------------------------
  # Ensure categorical variables are factors
  # --------------------------------------------------------
  database <- lapply(database, function(df) {
    df[] <- lapply(df, function(col) {
      if (is.character(col)) as.factor(col) else col
    })
    df
  })

  # -------------------------------------------------------
  # Data type detection
  # --------------------------------------------------------
  has_factors <- any(sapply(database[[1]], is.factor))

  all_factors <- all(sapply(database[[1]], is.factor))

  numeric_methods <- c(hierarchical_methods, "kmeans", "fuzzy", "mclust")
  mixed_methods   <- c("pam", "kprototypes")
  categorical_methods <- c("kmodes")

  # -------------------------------------------------------
  # Validation by method
  # --------------------------------------------------------
  if (all_factors) {

    if (method == "kprototypes") {
      warning(
        "All variables are categorical. Switching to 'kmodes'."
      )
      method <- "kmodes"
    }

    if (method %in% numeric_methods) {
      stop(
        "All variables are categorical. Use 'kmodes' or 'pam' ",
        "(with metric = 'gower') for categorical data."
      )
    }
  }

  if (method == "kprototypes" && scale_data == "none" && !all_factors) {
    warning(
      "No scaling applied. In k-prototypes, variables with larger numeric ",
      "scales may dominate the clustering solution. Consider using ",
      "scale_data = 'global' for improved balance and comparability."
    )
  }

  if (has_factors) {

    if (method %in% numeric_methods) {
      stop(
        "The dataset contains categorical variables. ",
        "Use 'pam' (with metric='gower') or 'kprototypes' for mixed data, ",
        "or 'kmodes' for purely categorical data."
      )
    }

  } else {

    if (method %in% categorical_methods) {
      stop("kmodes requires categorical data.")
    }
  }

  # Warning for mixed data without Gower
  if (method == "pam" &&
      has_factors &&
      (is.null(extra_args$metric) || extra_args$metric != "gower")) {

    warning(
      "Mixed data detected. Consider using metric = 'gower' ",
      "or method = 'kprototypes' for mixed data."
    )
  }

  # Validate the number of clusters (k):
  # Ensures meaningful clustering by enforcing 2 <= k < n

  if (any(k < 2)) {
    stop("k must be greater than or equal to 2.")
  }

  n_obs <- nrow(database[[1]])

  if (any(k >= n_obs)) {
    stop("k must be less than the number of observations.")
  }

  # --------------------------------------------------------
  # Performing clustering for each imputation
  # --------------------------------------------------------
  # --------------------------------------------------------
  # If k is a single value: keep original behavior
  # --------------------------------------------------------

  if (length(k) == 1) {

    partitions <- lapply(database, function(df) {

      if (method %in% hierarchical_methods) {
        d <- stats::dist(df)
        hc <- stats::hclust(d, method = method)
        stats::cutree(hc, k = k)

      } else if (method == "kmeans") {
        stats::kmeans(df, centers = k, ...)$cluster

      } else if (method == "pam") {

        if (!is.null(extra_args$metric) && extra_args$metric == "gower") {

          # Gower distance for mixed data
          d <- cluster::daisy(df, metric = "gower")
          cluster::pam(d, k, diss = TRUE)$clustering

        } else {

          # Classic PAM (numeric data)
          cluster::pam(df, k, ...)$clustering
        }

      } else if (method == "fuzzy") {

        e1071::cmeans(df, centers = k, ...)$cluster

      } else if (method == "mclust") {

        mclust::Mclust(df, G = k, ...)$classification

      } else if (method == "kmodes") {

        klaR::kmodes(df, modes = k, ...)$cluster

      } else if (method == "kprototypes") {

        clustMixType::kproto(df, k = k, ...)$cluster
      }

    })

    names(partitions) <- names(imp_list)
    return(partitions)
  }

  # --------------------------------------------------------
  # If k is a range/vector: compute partitions for each k
  # Return: list by k, each element is the same output as above
  # --------------------------------------------------------

  if (length(k) > 1) {

    # ------------------------------------------
    # Case 1: Hierarchical (optimized)
    # ------------------------------------------
    if (method %in% hierarchical_methods) {

      # Compute dendrogram once per imputation
      hc_list <- lapply(database, function(df) {
        d <- stats::dist(df)
        stats::hclust(d, method = method)
      })

      partitions_by_k <- lapply(k, function(k_i) {

        partitions <- lapply(hc_list, function(hc) {
          stats::cutree(hc, k = k_i)
        })

        names(partitions) <- names(imp_list)
        partitions
      })

    } else {

      # ------------------------------------------
      # Case 2: Non-hierarchical methods
      # ------------------------------------------

      partitions_by_k <- lapply(k, function(k_i) {

        partitions <- lapply(database, function(df) {

          if (method == "kmeans") {

            stats::kmeans(df, centers = k_i, ...)$cluster

          } else if (method == "pam") {

            if (!is.null(extra_args$metric) &&
                extra_args$metric == "gower") {

              d <- cluster::daisy(df, metric = "gower")
              cluster::pam(d, k_i, diss = TRUE)$clustering

            } else {

              cluster::pam(df, k_i, ...)$clustering
            }

          } else if (method == "fuzzy") {

            e1071::cmeans(df, centers = k_i, ...)$cluster

          } else if (method == "mclust") {

            mclust::Mclust(df, G = k_i, ...)$classification

          } else if (method == "kmodes") {

            klaR::kmodes(df, modes = k_i, ...)$cluster

          } else if (method == "kprototypes") {

            clustMixType::kproto(df, k = k_i, ...)$cluster
          }

        })

        names(partitions) <- names(imp_list)
        partitions
      })
    }

    names(partitions_by_k) <- paste0("k", k)
    return(partitions_by_k)
  }
}
