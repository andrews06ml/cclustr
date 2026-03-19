#' Plot a consensus clustering dendrogram
#'
#' @description
#' Displays the hierarchical clustering dendrogram derived from the consensus
#' dissimilarity matrix \eqn{1 - C}, where \eqn{C} is the co-assignment matrix
#' produced by \code{\link{consensus_clustering}}. Optionally draws colored
#' rectangles around the clusters defined by the requested \code{k}.
#'
#' @param consensus_result A named list as returned by
#'   \code{\link{consensus_clustering}} for a single \code{k}. Must contain
#'   at least \code{$hclust} (an \code{hclust} object) and \code{$k} (an
#'   integer specifying the number of clusters).
#' @param rect Logical. If \code{TRUE} (default), colored rectangles are drawn
#'   around the clusters via \code{\link[stats]{rect.hclust}}.
#' @param hang Numeric. The fraction of the plot height by which labels hang
#'   below the rest of the plot. A negative value (default \code{-1}) causes
#'   all labels to hang down from zero.
#'
#' @return Invisibly returns \code{NULL}. The function is called for its
#'   side effect of producing a plot.
#'
#' @details
#' The dendrogram is built from the \code{hclust} object stored in
#' \code{consensus_result$hclust}, which was computed by
#' \code{\link{consensus_clustering}} on the consensus dissimilarity matrix
#' \eqn{1 - C}. The colored rectangles highlight the \code{k} groups
#' obtained by cutting the dendrogram at the level that yields exactly
#' \code{k} clusters.
#'
#' This function is intended to be called on a single \code{k}. To compare
#' dendrograms across multiple values of \code{k}, call the function
#' iteratively and manage the graphics layout externally via
#' \code{\link[graphics]{par}}.
#'
#' @examples
#' \dontrun{
#' library(mice)
#'
#' imp    <- mice(nhanes, m = 5, printFlag = FALSE)
#' mild   <- as_mild_list(imp)
#' parts  <- cluster_imputations(mild, method = "ward.D2", k = 3)
#' cons   <- consensus_clustering(parts, k = 3)
#'
#' # Single dendrogram
#' plot_consensus_dendrogram(cons)
#'
#' # Without rectangles
#' plot_consensus_dendrogram(cons, rect = FALSE)
#'
#' # Compare multiple k values
#' parts_multi <- cluster_imputations(mild, method = "ward.D2", k = 2:4)
#' cons_multi  <- consensus_clustering(parts_multi)
#'
#' par(mfrow = c(1, 3), mar = c(3, 2, 2, 1))
#' for (k_i in 2:4) {
#'   plot_consensus_dendrogram(cons_multi[[paste0("k", k_i)]])
#' }
#' par(mfrow = c(1, 1))
#' }
#'
#' @seealso
#' \code{\link{consensus_clustering}}, \code{\link{plot_consensus_matrix}},
#' \code{\link{plot_validation_metrics}}
#'
#' @export
plot_consensus_dendrogram <- function(consensus_result,
                                      rect = TRUE,
                                      hang = -1) {
  
  if (is.null(consensus_result$hclust)) {
    stop("consensus_result must contain 'hclust'.")
  }
  
  if (is.null(consensus_result$k)) {
    stop("consensus_result must contain 'k'.")
  }
  
  plot(consensus_result$hclust,
       hang = hang,
       main = paste("Dendrograma consenso (k =", consensus_result$k, ")"),
       xlab = "",
       sub = "")
  
  if (rect) {
    rect.hclust(consensus_result$hclust, k = consensus_result$k, border = 2:6)
  }
}
