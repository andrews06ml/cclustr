
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
#' @param scale_data Logical. If \code{TRUE} (default), numeric columns are
#'   standardized to zero mean and unit variance prior to clustering.
#'   This option is ignored for methods that do not require scaling,
#'   such as \code{"kmodes"}.
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
#' \dontrun{
#' library(mice)
#'
#' imp   <- mice(nhanes, m = 5, printFlag = FALSE)
#' mild  <- as_mild_list(imp)
#'
#' # Single k with hierarchical clustering
#' parts <- cluster_imputations(mild, method = "ward.D2", k = 3)
#'
#' # Range of k values
#' parts_multi <- cluster_imputations(mild, method = "ward.D2", k = 2:5)
#'
#' # PAM with Gower distance (mixed data)
#' parts_gower <- cluster_imputations(mild, method = "pam", k = 3,
#'                                    metric = "gower")
#'
#' # k-modes for categorical data
#' parts_kmodes <- cluster_imputations(mild, method = "kmodes", k = 3)
#'
#' # k-prototypes for mixed data
#' parts_kproto <- cluster_imputations(mild, method = "kprototypes", k = 3)
#' }
#'
#' @seealso \code{\link{as_mild_list}}, \code{\link{consensus_clustering}}
#'
#' @export
cluster_imputations <- function(imp_list,
                                method = "ward.D2",
                                k,
                                scale_data = TRUE,
                                ...) {
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
  # Optional escalation
  # --------------------------------------------------------

  # Disable scaling for methods that do not require it
  if (method %in% c("kmodes") && scale_data) {
    warning("Scaling disabled for method = 'kmodes'.")
    scale_data <- FALSE
  }

  database <- if (scale_data) {
    lapply(imp_list, function(df) {
      num_cols <- sapply(df, is.numeric)
      df[, num_cols] <- scale(df[, num_cols])
      df
    })
  } else {
    imp_list
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
  if (all_factors && method == "kprototypes") {
    warning(
      "All variables are categorical. Switching to 'kmodes'."
    )
    method <- "kmodes"
  }

  if (method == "kprototypes" && scale_data && !all_factors) {
    warning("Scaling numeric variables may affect k-prototypes clustering.")
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
