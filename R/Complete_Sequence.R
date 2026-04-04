#' Run the full multiple-imputation clustering pipeline
#'
#' @description
#' Executes the complete multiple-imputation clustering pipeline in a single
#' call, integrating imputation standardization, per-imputation clustering,
#' consensus construction, validation, and optimal \code{k} selection.
#' Designed as a convenience wrapper around \code{\link{as_mild_list}},
#' \code{\link{cluster_imputations}}, \code{\link{consensus_clustering}},
#' \code{\link{validate_clustering}}, and \code{\link{choose_best_clustering}}.
#' When a single value of \code{k} is supplied, the selection step is skipped
#' and the unique solution is returned directly.
#'
#' @param data An imputation object accepted by \code{\link{as_mild_list}}.
#'   Supported formats include \code{mids} objects from \pkg{mice},
#'   long-format data frames with a \code{.imp} column, \code{amelia} objects,
#'   \code{imputationList} objects from \pkg{mitools}, and plain lists of
#'   \code{data.frame}s.
#' @param method A character string specifying the clustering algorithm applied
#'   to each imputed dataset. For hierarchical clustering, accepted values are
#'   \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"},
#'   \code{"average"}, \code{"centroid"}, \code{"median"}, and
#'   \code{"mcquitty"}. Additional options are \code{"kmeans"}, \code{"pam"},
#'   \code{"fuzzy"}, and \code{"mclust"}. Default is \code{"ward.D2"}.
#' @param k A single integer or an integer vector specifying the number(s) of
#'   clusters to evaluate. If a single value is supplied, the pipeline skips
#'   the \code{\link{choose_best_clustering}} step and returns that solution
#'   directly. If a vector is supplied, all values are evaluated and the
#'   optimal \code{k} is selected via weighted rank aggregation.
#' @param scale_data Logical. If \code{TRUE} (default), numeric columns are
#'   standardized to zero mean and unit variance prior to clustering.
#' @param consensus_method A character string specifying the method used to
#'   build the co-assignment matrix. \code{"classic"} assigns equal weight to
#'   all partitions; \code{"weighted_ari"} weights partitions by their
#'   pairwise ARI centrality. Default is \code{"classic"}.
#' @param cluster_method_consensus A character string specifying the
#'   agglomeration method passed to \code{\link[stats]{hclust}} during the
#'   consensus stage. Accepted values are the same as for \code{method} when
#'   using hierarchical clustering. Default is \code{"ward.D2"}.
#' @param pac_lower A numeric value between 0 and 1 specifying the lower bound
#'   of the ambiguous assignment region for the PAC metric. Default is
#'   \code{0.1}.
#' @param pac_upper A numeric value between 0 and 1 specifying the upper bound
#'   of the ambiguous assignment region for the PAC metric. Default is
#'   \code{0.9}.
#' @param weights A named numeric vector specifying the relative importance of
#'   each metric in the rank aggregation used by
#'   \code{\link{choose_best_clustering}}. Expected names are \code{pac},
#'   \code{silhouette}, \code{ari_between}, \code{ari_consensus}, \code{ch},
#'   \code{db}, and \code{dunn}. If \code{NULL} (default), weights are set
#'   automatically based on \code{prefer_stability}. Ignored when \code{k} is
#'   a single value.
#' @param prefer_stability Logical. If \code{TRUE} (default), stability
#'   metrics (\code{pac}, \code{ari_between}, \code{ari_consensus}) receive
#'   higher weights during \code{k} selection. If \code{FALSE}, internal
#'   compactness metrics (\code{silhouette}, \code{ch}, \code{db},
#'   \code{dunn}) receive higher weights. Ignored when \code{k} is a single
#'   value or when \code{weights} is provided explicitly.
#' @param tie_breaker A character string specifying the secondary metric used
#'   to break ties during \code{k} selection. One of \code{"silhouette"}
#'   (default), \code{"pac"}, \code{"dunn"}, \code{"ch"}, \code{"db"},
#'   \code{"ari_between"}, or \code{"ari_consensus"}. Ignored when \code{k}
#'   is a single value.
#' @param ... Additional arguments passed to the underlying clustering
#'   function via \code{\link{cluster_imputations}}. For example,
#'   \code{metric = "gower"} enables Gower distance when
#'   \code{method = "pam"}.
#'
#' @return An object of class \code{"mi_clustering_result"}: a named list
#'   with the following elements:
#'   \itemize{
#'     \item \code{call}: the matched call.
#'     \item \code{input_k}: the value(s) of \code{k} supplied by the user.
#'     \item \code{clustering_method}: character string with the clustering
#'           method used.
#'     \item \code{consensus_method}: character string with the consensus
#'           method used.
#'     \item \code{scale_data}: logical indicating whether data were scaled.
#'     \item \code{imputations}: standardized list of imputed datasets, as
#'           returned by \code{\link{as_mild_list}}.
#'     \item \code{partitions}: clustering assignments per imputed dataset,
#'           as returned by \code{\link{cluster_imputations}}.
#'     \item \code{consensus_results}: consensus clustering output, as
#'           returned by \code{\link{consensus_clustering}}.
#'     \item \code{validation_table}: validation metrics per \code{k}, as
#'           returned by \code{\link{validate_clustering}}.
#'     \item \code{best_k}: integer, the selected optimal number of clusters.
#'     \item \code{best_consensus}: integer vector of consensus cluster labels
#'           for the selected \code{k}, one per observation.
#'     \item \code{best_coassignment}: numeric matrix of co-assignment
#'           probabilities for the selected \code{k}.
#'     \item \code{best_consensus_result}: the full consensus result object
#'           for the selected \code{k}.
#'     \item \code{scores_table}: validation table with an additional
#'           \code{score} column when \code{k} is a vector; the plain
#'           validation table when \code{k} is a single value.
#'     \item \code{selection}: the full output of
#'           \code{\link{choose_best_clustering}} when \code{k} is a vector,
#'           or a minimal equivalent list when \code{k} is a single value.
#'   }
#'
#' @details
#' The pipeline proceeds as follows:
#' \enumerate{
#'   \item \code{\link{as_mild_list}} standardizes the imputation object and
#'         validates data quality.
#'   \item \code{\link{cluster_imputations}} applies the chosen algorithm to
#'         each imputed dataset, optionally scaling the data beforehand.
#'   \item \code{\link{consensus_clustering}} builds a co-assignment matrix
#'         and derives the consensus partition for each requested \code{k}.
#'   \item \code{\link{validate_clustering}} computes internal and stability
#'         metrics for each candidate \code{k}.
#'   \item \code{\link{choose_best_clustering}} selects the optimal \code{k}
#'         via weighted rank aggregation (only when \code{length(k) > 1}).
#' }
#'
#' When \code{length(k) == 1}, step 5 is bypassed and the unique solution
#' is returned directly, with \code{weights} and \code{tie_breaker} set to
#' \code{NULL} in the \code{selection} element.
#'
#' @examples
#' \donttest{
#' library(mice)
#'
#' imp <- mice(nhanes, m = 5, printFlag = FALSE)
#'
#' # Single k
#' res <- run_mi_clustering(imp, method = "ward.D2", k = 3)
#' res$best_k
#' res$best_consensus
#'
#' # Range of k values (automatic selection)
#' res_multi <- run_mi_clustering(imp, method = "ward.D2", k = 2:5)
#' res_multi$best_k
#' res_multi$validation_table
#'
#' # PAM with Gower distance (mixed data)
#' res_gower <- run_mi_clustering(imp, method = "pam", k = 2:4,
#'                                metric = "gower")
#'
#' # ARI-weighted consensus, compactness-focused selection
#' res_w <- run_mi_clustering(imp, method = "ward.D2", k = 2:5,
#'                            consensus_method  = "weighted_ari",
#'                            prefer_stability  = FALSE,
#'                            tie_breaker       = "dunn")
#' }
#'
#' @seealso
#' \code{\link{as_mild_list}}, \code{\link{cluster_imputations}},
#' \code{\link{consensus_clustering}}, \code{\link{validate_clustering}},
#' \code{\link{choose_best_clustering}}
#'
#' @export
run_mi_clustering <- function(data,
                              method = "ward.D2",
                              k,
                              scale_data = TRUE,
                              consensus_method = c("classic", "weighted_ari"),
                              cluster_method_consensus = "ward.D2",
                              pac_lower = 0.1,
                              pac_upper = 0.9,
                              weights = NULL,
                              prefer_stability = TRUE,
                              tie_breaker = c("silhouette", "pac", "dunn", "ch", "db",
                                              "ari_between", "ari_consensus"),
                              ...) {

  # --------------------------------------------------------
  # Match arguments
  # --------------------------------------------------------
  consensus_method <- match.arg(consensus_method)
  tie_breaker <- match.arg(tie_breaker)

  # --------------------------------------------------------
  # 1) Standardize imputations
  # --------------------------------------------------------
  imp_list <- as_mild_list(data)

  # --------------------------------------------------------
  # 2) Cluster each imputation
  # --------------------------------------------------------
  partitions <- cluster_imputations(
    imp_list = imp_list,
    method = method,
    k = k,
    scale_data = scale_data,
    ...
  )

  # --------------------------------------------------------
  # 3) Consensus clustering
  # --------------------------------------------------------
  consensus_results <- consensus_clustering(
    partitions = partitions,
    k = if (length(k) == 1) k else NULL,
    cluster_method = cluster_method_consensus,
    consensus_method = consensus_method
  )

  # --------------------------------------------------------
  # 4) Validation
  # --------------------------------------------------------
  validation_table <- validate_clustering(
    partitions = partitions,
    consensus_results = consensus_results,
    pac_lower = pac_lower,
    pac_upper = pac_upper
  )

  # --------------------------------------------------------
  # 5) Choose best solution (solo si k es un rango)
  # --------------------------------------------------------
  if (length(k) > 1) {
    best_solution <- choose_best_clustering(
      validation_table = validation_table,
      consensus_results = consensus_results,
      weights = weights,
      prefer_stability = prefer_stability,
      tie_breaker = tie_breaker
    )
  } else {
    # Con k único, la "mejor" solución es la única disponible
    best_solution <- list(
      best_k             = k,
      best_consensus     = consensus_results$consensus,
      best_coassignment  = consensus_results$coassignment,
      best_consensus_result = consensus_results,
      scores_table       = validation_table,
      weights            = NULL,
      tie_breaker        = NULL
    )
  }
  # --------------------------------------------------------
  # Return all results
  # --------------------------------------------------------
  out <- list(
    call = match.call(),
    input_k = k,
    clustering_method = method,
    consensus_method = consensus_method,
    scale_data = scale_data,
    imputations = imp_list,
    partitions = partitions,
    consensus_results = consensus_results,
    validation_table = validation_table,
    best_k = best_solution$best_k,
    best_consensus = best_solution$best_consensus,
    best_coassignment = best_solution$best_coassignment,
    best_consensus_result = best_solution$best_consensus_result,
    scores_table = best_solution$scores_table,
    selection = best_solution
  )

  class(out) <- "mi_clustering_result"
  return(out)
}
