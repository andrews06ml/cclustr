#' Select the optimal number of clusters from a validation table
#'
#' @description
#' Identifies the best value of \code{k} from the validation metrics produced
#' by \code{\link{validate_clustering}}, using a weighted rank aggregation
#' strategy. Each metric is ranked independently and then combined into a
#' single score using user-defined weights. The \code{k} with the lowest
#' weighted rank score is selected as optimal. In case of ties, a secondary
#' metric is used as a tiebreaker.
#'
#' @param validation_table A \code{data.frame} as returned by
#'   \code{\link{validate_clustering}}, containing one row per evaluated
#'   \code{k} and columns \code{pac}, \code{silhouette_mean},
#'   \code{ari_mean_between_imputations}, \code{ari_consensus_mean},
#'   \code{calinski_harabasz_mean}, \code{davies_bouldin_mean}, and
#'   \code{dunn_index}.
#' @param consensus_results A list as returned by
#'   \code{\link{consensus_clustering}}, either for a single \code{k} or
#'   for multiple \code{k} values. Used to retrieve the consensus labels
#'   and co-assignment matrix for the selected \code{k}.
#' @param weights A named numeric vector specifying the relative importance
#'   of each metric in the rank aggregation. Expected names are \code{pac},
#'   \code{silhouette}, \code{ari_between}, \code{ari_consensus}, \code{ch},
#'   \code{db}, and \code{dunn}. If \code{NULL} (default), weights are set
#'   automatically based on \code{prefer_stability}. When provided
#'   explicitly, \code{prefer_stability} is ignored.
#' @param prefer_stability Controls the default weighting strategy when
#'   \code{weights} is \code{NULL}. Accepted values are:
#'   \itemize{
#'     \item \code{NULL} (default): all metrics receive equal weight,
#'           providing a balanced selection that does not favor any
#'           particular aspect of clustering quality.
#'     \item \code{TRUE}: stability metrics (\code{pac},
#'           \code{ari_between}, \code{ari_consensus}) receive higher
#'           weights.
#'     \item \code{FALSE}: internal compactness metrics
#'           (\code{silhouette}, \code{ch}, \code{db}, \code{dunn})
#'           receive higher weights.
#'   }
#' @param tie_breaker A character string specifying the secondary metric
#'   used to break ties in the weighted rank score. One of
#'   \code{"silhouette"} (default), \code{"pac"}, \code{"dunn"},
#'   \code{"ch"}, \code{"db"}, \code{"ari_between"}, or
#'   \code{"ari_consensus"}.
#'
#' @return A named list with the following elements:
#'   \itemize{
#'     \item \code{best_k}: integer, the selected optimal number of clusters.
#'     \item \code{best_consensus}: integer vector of consensus cluster
#'           labels for the selected \code{k}, one per observation.
#'     \item \code{best_coassignment}: numeric matrix of size
#'           \eqn{n \times n} with co-assignment probabilities for the
#'           selected \code{k}.
#'     \item \code{best_consensus_result}: the full consensus result object
#'           for the selected \code{k}, as returned by
#'           \code{\link{consensus_clustering}}.
#'     \item \code{scores_table}: the input \code{validation_table} with an
#'           additional \code{score} column containing the weighted rank
#'           score for each \code{k}, sorted by \code{k}.
#'     \item \code{weights}: named numeric vector of weights effectively
#'           used in the aggregation.
#'     \item \code{tie_breaker}: character string indicating the tiebreaker
#'           metric used.
#'   }
#'
#' @details
#' The rank aggregation proceeds as follows:
#' \enumerate{
#'   \item Each metric is ranked from best to worst (rank 1 = best),
#'         accounting for the direction of optimality: lower is better for
#'         \code{pac} and \code{davies_bouldin_mean}; higher is better for
#'         all remaining metrics.
#'   \item For each \code{k}, the weighted average rank is computed using
#'         the provided or default weights. Rows with \code{NA} in some
#'         metrics are handled by renormalizing the weights over available
#'         metrics only.
#'   \item The \code{k} with the lowest weighted rank score is selected.
#'         If multiple \code{k} values share the minimum score, the
#'         tiebreaker metric is used; if still tied, the smallest \code{k}
#'         is preferred.
#' }
#'
#' Default weights when \code{prefer_stability = NULL} (default):
#' \code{pac = 1}, \code{silhouette = 1}, \code{ari_between = 1},
#' \code{ari_consensus = 1}, \code{ch = 1}, \code{db = 1}, \code{dunn = 1}.
#'
#' Default weights when \code{prefer_stability = TRUE}:
#' \code{pac = 2}, \code{silhouette = 1.5}, \code{ari_between = 2},
#' \code{ari_consensus = 2}, \code{ch = 1}, \code{db = 1}, \code{dunn = 1}.
#'
#' Default weights when \code{prefer_stability = FALSE}:
#' \code{pac = 1}, \code{silhouette = 2}, \code{ari_between = 1},
#' \code{ari_consensus = 1}, \code{ch = 2}, \code{db = 2}, \code{dunn = 2}.
#'
#' No automatic selection method replaces domain knowledge. The
#' \code{scores_table} and \code{\link{plot_validation_metrics}} are
#' provided so that the user can inspect individual metrics and make an
#' informed decision.
#'
#' @examples
#' \dontrun{
#' library(mice)
#'
#' imp   <- mice(nhanes, m = 5, printFlag = FALSE)
#' mild  <- as_mild_list(imp)
#' parts <- cluster_imputations(mild, method = "ward.D2", k = 2:5)
#' cons  <- consensus_clustering(parts)
#' val   <- validate_clustering(parts, cons)
#'
#' # Default selection (equal weights)
#' best <- choose_best_clustering(val, cons)
#' best$best_k
#'
#' # Stability-focused selection
#' best2 <- choose_best_clustering(val, cons,
#'                                 prefer_stability = TRUE)
#'
#' # Compactness-focused selection
#' best3 <- choose_best_clustering(val, cons,
#'                                 prefer_stability = FALSE,
#'                                 tie_breaker      = "dunn")
#'
#' # Custom weights
#' best4 <- choose_best_clustering(val, cons,
#'                                 weights = c(pac = 3, silhouette = 1,
#'                                             ari_between = 3, ari_consensus = 3,
#'                                             ch = 1, db = 1, dunn = 1))
#' }
#'
#' @references
#' Pihur, V., Datta, S., & Datta, S. (2007). Weighted rank aggregation
#' of cluster validation measures: a Monte Carlo cross-entropy approach.
#' \emph{Bioinformatics}, \strong{23}(13), 1607–1615.
#' \doi{10.1093/bioinformatics/btm158}
#'
#' @seealso \code{\link{validate_clustering}}, \code{\link{consensus_clustering}},
#'   \code{\link{plot_validation_metrics}}
#'
#' @export
choose_best_clustering <- function(validation_table,
                                   consensus_results,
                                   weights          = NULL,
                                   prefer_stability = NULL,
                                   tie_breaker      = c("silhouette", "pac", "dunn", "ch", "db",
                                                        "ari_between", "ari_consensus")) {

  tie_breaker <- match.arg(tie_breaker)

  # ---- required columns
  req  <- c("k", "pac", "silhouette_mean",
            "ari_mean_between_imputations", "ari_consensus_mean",
            "calinski_harabasz_mean", "davies_bouldin_mean", "dunn_index")
  miss <- setdiff(req, names(validation_table))
  if (length(miss) > 0) {
    stop(paste0("validation_table is missing: ", paste(miss, collapse = ", ")))
  }

  # ---- default weights
  if (is.null(weights)) {
    if (isTRUE(prefer_stability)) {
      message("Weighting strategy: stability-focused")
      weights <- c(pac = 2, silhouette = 1.5,
                   ari_between = 2, ari_consensus = 2,
                   ch = 1, db = 1, dunn = 1)

    } else if (isFALSE(prefer_stability)) {
      message("Weighting strategy: compactness-focused")
      weights <- c(pac = 1, silhouette = 2,
                   ari_between = 1, ari_consensus = 1,
                   ch = 2, db = 2, dunn = 2)

    } else {
      message("Weighting strategy: equal weights")
      weights <- c(pac = 1, silhouette = 1,
                   ari_between = 1, ari_consensus = 1,
                   ch = 1, db = 1, dunn = 1)
    }
  }

  .rank_best1 <- function(x, higher_is_better = TRUE) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    if (higher_is_better) rank(-x, ties.method = "average", na.last = "keep")
    else                  rank( x, ties.method = "average", na.last = "keep")
  }

  vt <- validation_table

  r_pac  <- .rank_best1(vt$pac,                          higher_is_better = FALSE)
  r_sil  <- .rank_best1(vt$silhouette_mean,              higher_is_better = TRUE)
  r_ab   <- .rank_best1(vt$ari_mean_between_imputations, higher_is_better = TRUE)
  r_ac   <- .rank_best1(vt$ari_consensus_mean,           higher_is_better = TRUE)
  r_ch   <- .rank_best1(vt$calinski_harabasz_mean,       higher_is_better = TRUE)
  r_db   <- .rank_best1(vt$davies_bouldin_mean,          higher_is_better = FALSE)
  r_dunn <- .rank_best1(vt$dunn_index,                   higher_is_better = TRUE)

  score <- numeric(nrow(vt))
  for (i in seq_len(nrow(vt))) {
    vals <- c(pac           = r_pac[i],
              silhouette    = r_sil[i],
              ari_between   = r_ab[i],
              ari_consensus = r_ac[i],
              ch            = r_ch[i],
              db            = r_db[i],
              dunn          = r_dunn[i])
    ok <- !is.na(vals)
    if (!any(ok)) { score[i] <- NA_real_; next }
    wi       <- weights[names(vals)][ok]
    wi       <- wi / sum(wi)
    score[i] <- sum(wi * vals[ok])
  }
  vt$score <- score

  tb <- switch(tie_breaker,
               silhouette    =  vt$silhouette_mean,
               pac           = -vt$pac,
               dunn          =  vt$dunn_index,
               ch            =  vt$calinski_harabasz_mean,
               db            = -vt$davies_bouldin_mean,
               ari_between   =  vt$ari_mean_between_imputations,
               ari_consensus =  vt$ari_consensus_mean)

  best_idx <- order(vt$score, -tb, vt$k)[1]
  best_k   <- vt$k[best_idx]

  best_obj <- NULL
  if (!is.null(consensus_results$k) && !is.null(consensus_results$consensus)) {
    best_obj <- consensus_results
  } else {
    key <- paste0("k", best_k)
    if (!is.null(consensus_results[[key]])) {
      best_obj <- consensus_results[[key]]
    } else {
      hit <- which(sapply(consensus_results, function(z) !is.null(z$k) && z$k == best_k))
      if (length(hit) == 1) best_obj <- consensus_results[[hit]]
    }
  }

  if (is.null(best_obj) || is.null(best_obj$consensus) || is.null(best_obj$coassignment)) {
    stop("Could not retrieve consensus labels/coassignment for the selected k from consensus_results.")
  }

  message(paste("✓ Best k selected:", best_k))

  list(
    best_k                = best_k,
    best_consensus        = best_obj$consensus,
    best_coassignment     = best_obj$coassignment,
    best_consensus_result = best_obj,
    scores_table          = vt[order(vt$k), ],
    weights               = weights,
    tie_breaker           = tie_breaker
  )
}
