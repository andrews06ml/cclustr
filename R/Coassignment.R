#' Build a consensus partition from multiple imputation clustering results
#'
#' @description
#' Constructs a consensus clustering solution from a collection of partitions
#' obtained across multiple imputed datasets (as produced by
#' \code{\link{cluster_imputations}}). The function builds a co-assignment
#' matrix reflecting how frequently each pair of observations is assigned to
#' the same cluster, and derives the final consensus partition via hierarchical
#' clustering on the resulting dissimilarity. Supports both a single value of
#' \code{k} and a range of values.
#'
#' @param partitions A list of cluster assignment vectors (one per imputed
#'   dataset) for a single \code{k}, or a named list of such lists (one per
#'   \code{k} value, named \code{"k2"}, \code{"k3"}, etc.) as returned by
#'   \code{\link{cluster_imputations}} when \code{k} is a vector.
#' @param k An integer specifying the number of consensus clusters. Required
#'   when \code{partitions} is a flat list of vectors (single-k input).
#'   Ignored when \code{partitions} is a list-of-lists with names following
#'   the pattern \code{"k2"}, \code{"k3"}, etc., as \code{k} is inferred
#'   automatically from the names.
#' @param cluster_method A character string specifying the agglomeration
#'   method passed to \code{\link[stats]{hclust}} for the consensus stage.
#'   Accepted values are \code{"ward.D"}, \code{"ward.D2"}, \code{"single"},
#'   \code{"complete"}, \code{"average"}, \code{"centroid"}, \code{"median"},
#'   and \code{"mcquitty"}. Default is \code{"ward.D2"}.
#' @param consensus_method A character string specifying how the
#'   co-assignment matrix is built. Options are:
#'   \itemize{
#'     \item \code{"classic"}: unweighted co-assignment; all partitions
#'           contribute equally.
#'     \item \code{"weighted_ari"}: partitions are weighted by their
#'           pairwise ARI centrality, so more consistent partitions
#'           contribute more to the consensus. Requires \pkg{mclust}.
#'   }
#'   Default is \code{"classic"}.
#'
#' @return
#' \itemize{
#'   \item If \code{partitions} is a flat list (single \code{k}): a named
#'         list with the following elements:
#'     \itemize{
#'       \item \code{consensus_method}: character string with the method used.
#'       \item \code{k}: integer, number of clusters requested.
#'       \item \code{consensus}: integer vector of consensus cluster labels,
#'             one per observation.
#'       \item \code{coassignment}: numeric matrix of size \eqn{n \times n}
#'             with co-assignment probabilities between 0 and 1.
#'       \item \code{hclust}: the \code{hclust} object from the consensus
#'             stage.
#'       \item \code{weights}: named numeric vector of partition weights.
#'     }
#'   \item If \code{partitions} is a list-of-lists (multiple \code{k}):
#'         a named list where each element corresponds to one value of
#'         \code{k} and contains the structure described above.
#' }
#'
#' @details
#' The co-assignment matrix \eqn{C} is defined such that \eqn{C_{ij}}
#' represents the weighted proportion of partitions in which observations
#' \eqn{i} and \eqn{j} are assigned to the same cluster. The final
#' dissimilarity is computed as \eqn{1 - C} and passed to
#' \code{\link[stats]{hclust}}.
#'
#' When \code{consensus_method = "weighted_ari"}, the weight of each
#' partition is proportional to its mean ARI against all other partitions,
#' reflecting its centrality within the ensemble. Negative centrality values
#' are floored at zero. If all centralities are zero, the function falls back
#' to equal weights with a warning.
#'
#' @examples
#' \dontrun{
#' library(mice)
#'
#' imp    <- mice(nhanes, m = 5, printFlag = FALSE)
#' mild   <- as_mild_list(imp)
#' parts  <- cluster_imputations(mild, method = "ward.D2", k = 3)
#'
#' # Single k, classic consensus
#' cons <- consensus_clustering(parts, k = 3)
#'
#' # Single k, ARI-weighted consensus
#' cons_w <- consensus_clustering(parts, k = 3,
#'                                consensus_method = "weighted_ari")
#'
#' # Multiple k values
#' parts_multi <- cluster_imputations(mild, method = "ward.D2", k = 2:5)
#' cons_multi  <- consensus_clustering(parts_multi)
#' }
#'
#' @seealso \code{\link{cluster_imputations}}, \code{\link{validate_clustering}}
#'
#' @export
consensus_clustering <- function(partitions,
                                 k = NULL,
                                 cluster_method = "ward.D2",
                                 consensus_method = c("classic", "weighted_ari")) {

  # --------------------------------------------------------
  # Supported methods for hclust consensus stage
  # --------------------------------------------------------
  hierarchical_methods <- c(
    "ward.D", "ward.D2", "single", "complete",
    "average", "centroid", "median", "mcquitty"
  )

  if (!(cluster_method %in% hierarchical_methods)) {
    stop("cluster_method not supported for hclust().")
  }

  # --------------------------------------------------------
  # Input validation: partitions must be a list
  # --------------------------------------------------------
  if (!is.list(partitions) || length(partitions) == 0) {
    stop("Partitions must be provided as a non-empty list.")
  }

  consensus_method <- match.arg(consensus_method)

  # --------------------------------------------------------
  # Detect input structure
  # Case A: list of vectors  -> single k
  # Case B: list of lists    -> multiple k (output of cluster_imputations when k is a range)
  # --------------------------------------------------------
  is_list_of_lists <- all(sapply(partitions, is.list))

  # --------------------------------------------------------
  # Helper: run consensus for a single k (core logic)
  # --------------------------------------------------------
  .consensus_single_k <- function(partitions_k, k_value) {

    # partitions_k must be a list of vectors
    if (!is.list(partitions_k) || length(partitions_k) == 0) {
      stop("Each k element must be a non-empty list of vectors.")
    }

    # Convert partitions to atomic vectors (defensive)
    partitions_k <- lapply(partitions_k, function(z) {
      if (is.data.frame(z) || is.matrix(z)) z <- as.vector(z)
      if (is.factor(z)) z <- as.character(z)
      if (!is.atomic(z)) stop("Each partition must be an atomic vector.")
      z
    })

    # n observations and m partitions
    n <- length(partitions_k[[1]])
    m <- length(partitions_k)

    if (n == 0) stop("Partitions contain empty vectors.")

    # Consistency check: all partitions must have same length
    if (!all(sapply(partitions_k, length) == n)) {
      stop("All partitions must have the same length.")
    }

    # Validate k
    if (length(k_value) != 1 || !is.numeric(k_value) || is.na(k_value) || k_value < 2) {
      stop("k must be a single numeric value >= 2.")
    }
    k_value <- as.integer(k_value)

    # --------------------------------------------------------
    # Compute weights
    # --------------------------------------------------------
    weights <- rep(1 / m, m)
    names(weights) <- names(partitions_k)

    if (consensus_method == "classic") {
      # do nothing; equal weights
      message(paste0("Detected consensus method: classic (k = ", k_value, ")"))

    } else if (consensus_method == "weighted_ari") {
      message(paste0("Detected consensus method: weighted_ari (k = ", k_value, ")"))

      if (!requireNamespace("mclust", quietly = TRUE)) {
        stop("Package 'mclust' is required for consensus_method = 'weighted_ari'.")
      }

      # Pairwise ARI matrix (m x m)
      ari_mat <- matrix(NA_real_, nrow = m, ncol = m)

      for (i in seq_len(m)) {
        ari_mat[i, i] <- 1
        for (j in seq_len(m)) {
          if (j > i) {
            ari <- mclust::adjustedRandIndex(partitions_k[[i]], partitions_k[[j]])
            ari_mat[i, j] <- ari
            ari_mat[j, i] <- ari
          }
        }
      }

      # Centrality score: mean ARI vs others (exclude self)
      centrality <- sapply(seq_len(m), function(i) mean(ari_mat[i, -i], na.rm = TRUE))
      centrality <- pmax(centrality, 0) # no negative weights

      if (all(centrality == 0)) {
        warning(paste0("All ARI centralities are 0 for k = ", k_value, ". Falling back to equal weights."))
        weights <- rep(1 / m, m)
      } else {
        weights <- centrality / sum(centrality)
      }

      names(weights) <- names(partitions_k)
    }

    # --------------------------------------------------------
    # Build (weighted) co-assignment matrix
    # --------------------------------------------------------
    coassignment <- matrix(0, n, n)

    for (idx in seq_len(m)) {
      z <- partitions_k[[idx]]
      coassignment <- coassignment + weights[idx] * outer(z, z, FUN = "==")
    }

    coassignment <- pmin(pmax(coassignment, 0), 1)

    # --------------------------------------------------------
    # Consensus clustering based on dissimilarity
    # --------------------------------------------------------
    diss <- stats::as.dist(1 - coassignment)
    hc <- stats::hclust(diss, method = cluster_method)
    consensus <- stats::cutree(hc, k = k_value)

    return(list(
      consensus_method = consensus_method,
      k = k_value,
      consensus = consensus,
      coassignment = coassignment,
      hclust = hc,
      weights = weights
    ))
  }

  # --------------------------------------------------------
  # Case A: partitions is a list of vectors (single k)
  # --------------------------------------------------------
  if (!is_list_of_lists) {

    if (is.null(k)) {
      stop("When 'partitions' is a list of vectors (single k), you must provide k.")
    }

    res <- .consensus_single_k(partitions, k)
    message("Consensus clustering completed")
    return(res)
  }

  # --------------------------------------------------------
  # Case B: partitions is a list-of-lists by k (multiple k)
  # Output of cluster_imputations() when k is a range: names like 'k2','k3',...
  # --------------------------------------------------------
  # Infer k values from names if possible; otherwise require k
  k_names <- names(partitions)

  if (!is.null(k_names) && all(grepl("^k[0-9]+$", k_names))) {
    k_values <- as.integer(sub("^k", "", k_names))
  } else {
    if (is.null(k)) {
      stop("For multiple-k input (list-of-lists), provide k or name elements as 'k2','k3',...")
    }
    k_values <- as.integer(k)
  }

  # Run consensus for each k
  results_by_k <- vector("list", length(partitions))
  names(results_by_k) <- names(partitions)

  for (idx in seq_along(partitions)) {
    k_val <- if (!is.null(k_names) && length(k_values) == length(partitions)) k_values[idx] else k_values[idx]
    results_by_k[[idx]] <- .consensus_single_k(partitions[[idx]], k_val)
  }

  message("Consensus clustering completed for multiple k values")
  return(results_by_k)
}
