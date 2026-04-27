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
#' @param k Integer. The number of clusters to display. Only required when
#'   \code{consensus_result} contains results for multiple \code{k} values.
#'   Ignored if a single-k result is passed directly.
#' @param hang Numeric. The fraction of the plot height by which labels hang
#'   below the rest of the plot. A negative value (default \code{-1}) causes
#'   all labels to hang down from zero.
#' @param labels Logical or character vector. Controls the labels displayed
#'   at the leaves of the dendrogram. If \code{NULL} (default), observation
#'   labels are shown. If \code{FALSE}, no labels are displayed — recommended
#'   for large datasets where labels become unreadable. A character vector
#'   of the same length as the number of observations can also be supplied
#'   to use custom labels.
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
#' # ------------------------------------------------------------
#' # Example 1: consensus dendrogram with simulated results
#' # ------------------------------------------------------------
#'
#' hc <- hclust(dist(matrix(rnorm(60), nrow = 20)), method = "ward.D2")
#' cons <- list(hclust = hc, k = 3)
#'
#' plot_consensus_dendrogram(cons)
#' plot_consensus_dendrogram(cons, rect = FALSE)
#'
#' \donttest{
#' # ------------------------------------------------------------
#' # Example 2: Dendrogram with mice
#' # ------------------------------------------------------------
#'
#' if (requireNamespace("mice", quietly = TRUE)) {
#'
#' imp    <- mice::mice(mice::nhanes, m = 5, printFlag = FALSE)
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
#' }
#'
#' @seealso
#' \code{\link{consensus_clustering}}, \code{\link{plot_consensus_matrix}},
#' \code{\link{plot_validation_metrics}}
#'
#' @export
plot_consensus_dendrogram <- function(consensus_result,
                                      k = NULL,
                                      rect = TRUE,
                                      hang = -1,
                                      labels = NULL) {

  # --------------------------------------------------------
  # Auto-extract best_consensus_result from mi_clustering_result
  # --------------------------------------------------------
  if (inherits(consensus_result, "mi_clustering_result")) {
    consensus_result <- consensus_result$best_consensus_result
  }

  # --------------------------------------------------------
  # Auto-extract best_consensus_result from choose_best_clustering()
  # --------------------------------------------------------
  if (is.list(consensus_result) &&
      !is.null(consensus_result$best_consensus_result) &&
      is.null(consensus_result$coassignment)) {
    consensus_result <- consensus_result$best_consensus_result
  }

  # --------------------------------------------------------
  # Detect if multiple k (list of lists)
  # --------------------------------------------------------
  is_multi_k <- is.list(consensus_result) &&
    all(sapply(consensus_result, function(x) is.list(x) && !is.null(x$hclust)))

  # --------------------------------------------------------
  # k: Available options
  # --------------------------------------------------------
  if (is_multi_k && is.null(k)) {
    available_k <- names(consensus_result)
    stop(paste0("Multiple k detected. Available options: ",
                paste(available_k, collapse = ", "),
                ". Please specify k."))
  }

  # --------------------------------------------------------
  # Case 1: multiple k
  # --------------------------------------------------------
  if (is_multi_k) {

    k_name <- paste0("k", k)

    if (!k_name %in% names(consensus_result))
      stop("Requested k not found in consensus_result.")

    consensus_result <- consensus_result[[k_name]]
  }

  # --------------------------------------------------------
  # Validate structure
  # --------------------------------------------------------
  if (is.null(consensus_result$hclust)) {
    stop("consensus_result must contain 'hclust'.")
  }

  if (is.null(consensus_result$k)) {
    stop("consensus_result must contain 'k'.")
  }

  # --------------------------------------------------------
  # Plot
  # --------------------------------------------------------
  # Build plot arguments dynamically
  plot_args <- list(
    x    = consensus_result$hclust,
    hang = hang,
    main = paste("Consensus dendrogram (k =", consensus_result$k, ")"),
    xlab = "",
    sub  = ""
  )

  # Only pass labels if explicitly specified by the user
  if (!is.null(labels)) {
    plot_args$labels <- labels
  }

  do.call(plot, plot_args)

  if (rect) {
    stats::rect.hclust(consensus_result$hclust,
                       k      = consensus_result$k,
                       border = 2:6)
  }
  invisible(NULL)
}
