#' Validate consensus clustering results across multiple imputed datasets
#'
#' @description
#' Computes a comprehensive set of internal and stability validation metrics
#' for consensus clustering results obtained from multiple imputed datasets.
#' All internal metrics are computed on the consensus dissimilarity matrix
#' \eqn{1 - C}, where \eqn{C} is the co-assignment matrix, making the
#' function robust to numeric, categorical, and mixed data types. Supports
#' both single-k and multi-k evaluation.
#'
#' @param partitions A list of cluster assignment vectors (one per imputed
#'   dataset) for a single \code{k}, or a named list of such lists (one per
#'   \code{k} value, named \code{"k2"}, \code{"k3"}, etc.) as returned by
#'   \code{\link{cluster_imputations}}.
#' @param consensus_results A list as returned by
#'   \code{\link{consensus_clustering}}, either for a single \code{k} or
#'   for multiple \code{k} values. Must contain at least \code{$consensus}
#'   and \code{$coassignment} for each \code{k}.
#' @param pac_lower A numeric value between 0 and 1 specifying the lower
#'   bound of the ambiguous assignment region for the PAC metric. Default
#'   is \code{0.1}.
#' @param pac_upper A numeric value between 0 and 1 specifying the upper
#'   bound of the ambiguous assignment region for the PAC metric. Default
#'   is \code{0.9}.
#'
#' @return A \code{data.frame} with one row per evaluated \code{k} and the
#'   following columns:
#'   \itemize{
#'     \item \code{k}: number of clusters evaluated.
#'     \item \code{pac}: Proportion of Ambiguous Clusterings. Lower values
#'           indicate a more stable consensus. Computed from the upper
#'           triangle of the co-assignment matrix.
#'     \item \code{silhouette_mean}: mean silhouette width computed on the
#'           consensus dissimilarity \eqn{1 - C}. Higher values indicate
#'           better-defined clusters.
#'     \item \code{ari_mean_between_imputations}: mean pairwise Adjusted
#'           Rand Index (ARI) across all pairs of imputed partitions.
#'           Higher values indicate greater clustering stability.
#'     \item \code{ari_consensus_mean}: mean ARI between the consensus
#'           partition and each individual imputed partition. Higher values
#'           indicate that the consensus is representative of the ensemble.
#'     \item \code{calinski_harabasz_mean}: Calinski-Harabasz index computed
#'           via \code{fpc::cluster.stats()}. Higher values indicate more
#'           compact and well-separated clusters.
#'     \item \code{davies_bouldin_mean}: Davies-Bouldin index computed using
#'           a medoid-based approach on the consensus dissimilarity. Lower
#'           values indicate better cluster separation.
#'     \item \code{dunn_index}: Dunn index computed via
#'           \code{fpc::cluster.stats()}. Higher values indicate better
#'           cluster compactness and separation.
#'   }
#'
#' @details
#' All internal validation metrics (silhouette, Calinski-Harabasz,
#' Davies-Bouldin, Dunn) are computed exclusively on the consensus
#' dissimilarity matrix \eqn{1 - C} rather than on the original feature
#' space. This design choice ensures compatibility with any data type and
#' avoids recomputing distances from the raw imputed datasets.
#'
#' The Davies-Bouldin index is computed using a custom medoid-based
#' implementation derived from the consensus dissimilarity, as no standard
#' R implementation accepts a precomputed dissimilarity matrix directly.
#'
#' Helper functions prefixed with a dot (e.g., \code{.pac_from_coassignment},
#' \code{.mean_silhouette_from_diss}) are internal and not exported.
#'
#' @examples
#' # ------------------------------------------------------------
#' # Example 1: Basic validation with simulated partitions
#' # ------------------------------------------------------------
#' set.seed(123)
#'
#' # Simulate partitions (3 imputations, k = 2)
#' parts <- list(
#'   imp1 = sample(1:2, 20, replace = TRUE),
#'   imp2 = sample(1:2, 20, replace = TRUE),
#'   imp3 = sample(1:2, 20, replace = TRUE)
#' )
#'
#' cons <- consensus_clustering(parts, k = 2)
#' val  <- validate_clustering(parts, cons)
#'
#' print(val)
#'
#' \donttest{
#' # ------------------------------------------------------------
#' # Example 2: Full pipeline with mice (multiple k)
#' # ------------------------------------------------------------
#' if (requireNamespace("mice", quietly = TRUE)) {
#'
#'   set.seed(123)
#'
#'   imp  <- mice::mice(mice::nhanes, m = 3, printFlag = FALSE)
#'   mild <- as_mild_list(imp)
#'
#'   parts <- cluster_imputations(mild, method = "ward.D2", k = 2:3)
#'   cons  <- consensus_clustering(parts)
#'
#'   val <- validate_clustering(parts, cons)
#'
#'   print(val)
#' }
#' }
#'
#' @references
#' Senbabaoglu, Y., Michailidis, G., & Li, J. Z. (2014).
#' Critical limitations of consensus clustering in class discovery.
#' \emph{Scientific Reports}, \strong{4}, 6207.
#' \doi{10.1038/srep06207}
#'
#' Rousseeuw, P. J. (1987).
#' Silhouettes: A graphical aid to the interpretation and validation
#' of cluster analysis.
#' \emph{Journal of Computational and Applied Mathematics}, \strong{20}, 53-65.
#' \doi{10.1016/0377-0427(87)90125-7}
#'
#' Calinski, T., & Harabasz, J. (1974).
#' A dendrite method for cluster analysis.
#' \emph{Communications in Statistics}, \strong{3}, 1-27.
#' \doi{10.1080/03610927408827101}
#'
#' Davies, D. L., & Bouldin, D. W. (1979).
#' A cluster separation measure.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
#' \strong{PAMI-1}, 224-227.
#' \doi{10.1109/TPAMI.1979.4766909}
#'
#' Dunn, J. C. (1974).
#' Well-separated clusters and optimal fuzzy partitions.
#' \emph{Journal of Cybernetics}, \strong{4}, 95-104.
#' \doi{10.1080/01969727408546059}
#'
#' Hubert, L., & Arabie, P. (1985).
#' Comparing partitions.
#' \emph{Journal of Classification}, \strong{2}, 193-218.
#' \doi{10.1007/BF01908075}
#'
#' @seealso \code{\link{consensus_clustering}}, \code{\link{choose_best_clustering}}
#'
#' @export
validate_clustering <- function(partitions,
                                consensus_results,
                                pac_lower = 0.1,
                                pac_upper = 0.9) {

  # Required packages
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required.")
  }
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required.")
  }
  if (!requireNamespace("fpc", quietly = TRUE)) {
    stop("Package 'fpc' is required.")
  }

  # --------------------------------------------------------
  # First validation
  # --------------------------------------------------------

  if (!is.list(partitions) || length(partitions) == 0) {
    stop("partitions must be a non-empty list.")
  }

  if (!is.list(consensus_results)) {
    stop("consensus_results must be a list returned by consensus_clustering().")
  }

  # --------------------------------------------------------
  # Helpers
  # --------------------------------------------------------

  .pac_from_coassignment <- function(coassignment, lower, upper) {
    vals <- coassignment[upper.tri(coassignment, diag = FALSE)]
    mean(vals > lower & vals < upper)
  }

  .mean_silhouette_from_diss <- function(diss, labels) {
    if (length(unique(labels)) < 2) return(NA_real_)

    sil <- tryCatch(
      cluster::silhouette(as.integer(labels), stats::as.dist(diss)),
      error = function(e) NULL
    )

    if (is.null(sil)) return(NA_real_)

    mean(sil[, "sil_width"])
  }

  .ari_mean_between_partitions <- function(parts_list) {
    m <- length(parts_list)
    if (m < 2) return(NA_real_)
    ari_vals <- c()
    for (i in seq_len(m - 1)) {
      for (j in seq.int(i + 1, m)) {
        ari_vals <- c(ari_vals,
                      mclust::adjustedRandIndex(parts_list[[i]], parts_list[[j]]))
      }
    }
    mean(ari_vals)
  }

  .ari_consensus_mean <- function(consensus, parts_list) {
    vals <- vapply(parts_list, function(z) {
      if (length(unique(z)) < 2 || length(unique(consensus)) < 2) {
        return(NA_real_)
      }
      mclust::adjustedRandIndex(consensus, z)
    }, numeric(1))
    mean(vals, na.rm = TRUE)
  }

  # CH + Dunn from dissimilarity
  .ch_dunn_from_diss <- function(diss, labels) {
    cs <- tryCatch(
      fpc::cluster.stats(d = stats::as.dist(diss),
                         clustering = as.integer(labels)),
      error = function(e) NULL
    )
    if (is.null(cs)) {
      return(list(ch = NA_real_, dunn = NA_real_))
    }
    list(
      ch = as.numeric(cs$ch),
      dunn = as.numeric(cs$dunn)
    )
  }

  # Davies–Bouldin (medoid-based)
  .db_from_diss_medoids <- function(diss, labels) {

    labs <- as.integer(labels)
    cl_ids <- sort(unique(labs))
    k <- length(cl_ids)
    if (k < 2) return(NA_real_)

    medoid <- integer(k)
    S <- numeric(k)

    for (ii in seq_along(cl_ids)) {
      cl <- cl_ids[ii]
      idx <- which(labs == cl)

      if (length(idx) <= 1) {
        medoid[ii] <- idx[1]
        S[ii] <- 0
        next
      }

      Dsub <- diss[idx, idx, drop = FALSE]
      mpos <- which.min(rowSums(Dsub))
      medoid[ii] <- idx[mpos]
      S[ii] <- mean(diss[idx, medoid[ii]])
    }

    d_med <- diss[medoid, medoid, drop = FALSE]

    Rmax <- rep(NA_real_, k)
    for (i in seq_len(k)) {
      ratios <- rep(NA_real_, k)
      for (j in seq_len(k)) {
        if (j == i) next
        denom <- d_med[i, j]
        if (!is.na(denom) && denom > 0) {
          ratios[j] <- (S[i] + S[j]) / denom
        }
      }
      Rmax[i] <- suppressWarnings(max(ratios, na.rm = TRUE))
      if (!is.finite(Rmax[i])) Rmax[i] <- NA_real_
    }

    mean(Rmax, na.rm = TRUE)
  }

  # --------------------------------------------------------
  # Core computation per k
  # --------------------------------------------------------

  .one_k <- function(parts_k, cres_k) {

    if (is.null(cres_k$consensus) || is.null(cres_k$coassignment)) {
      stop("Each consensus result must contain 'consensus' and 'coassignment'.")
    }

    cons  <- cres_k$consensus
    coass <- cres_k$coassignment
    diss  <- 1 - coass

    # Robust protection
    if (length(unique(cons)) < 2) {
      return(data.frame(
        k = cres_k$k,
        pac = NA_real_,
        silhouette_mean = NA_real_,
        ari_mean_between_imputations = NA_real_,
        ari_consensus_mean = NA_real_,
        calinski_harabasz_mean = NA_real_,
        davies_bouldin_mean = NA_real_,
        dunn_index = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    chd <- .ch_dunn_from_diss(diss, cons)

    data.frame(
      k = cres_k$k,
      pac = .pac_from_coassignment(coass, pac_lower, pac_upper),
      silhouette_mean = .mean_silhouette_from_diss(diss, cons),
      ari_mean_between_imputations =
        .ari_mean_between_partitions(parts_k),
      ari_consensus_mean =
        .ari_consensus_mean(cons, parts_k),
      calinski_harabasz_mean = chd$ch,
      davies_bouldin_mean =
        .db_from_diss_medoids(diss, cons),
      dunn_index = chd$dunn,
      stringsAsFactors = FALSE
    )
  }

  # --------------------------------------------------------
  # Single-k vs Multi-k
  # --------------------------------------------------------

  is_multi_k <- all(vapply(partitions, is.list, logical(1)))

  if (!is_multi_k) {
    if (!all(vapply(partitions, is.atomic, logical(1)))) {
      stop("For single-k input, partitions must be a list of atomic vectors.")
    }
    return(.one_k(partitions, consensus_results))
  }

  k_names <- intersect(names(partitions), names(consensus_results))
  k_values <- suppressWarnings(as.integer(sub("^k", "", k_names)))
  k_names <- k_names[order(k_values)]

  rows <- lapply(k_names,
                 function(kk)
                   .one_k(partitions[[kk]],
                          consensus_results[[kk]]))

  do.call(rbind, rows)
}
