#' Plot the consensus co-assignment matrix
#'
#' @description
#' Displays a heatmap of the co-assignment matrix \eqn{C} produced by
#' \code{\link{consensus_clustering}}, where each cell \eqn{C_{ij}} represents
#' the weighted proportion of partitions in which observations \eqn{i} and
#' \eqn{j} were assigned to the same cluster. Values close to 1 indicate
#' that two observations are consistently grouped together across all
#' imputed datasets; values close to 0 indicate the opposite. Rows and
#' columns can be reordered by consensus cluster assignment to make the
#' block structure visually apparent.
#'
#' @param consensus_result A named list as returned by
#'   \code{\link{consensus_clustering}} for a single \code{k}. Must contain
#'   at least \code{$coassignment} (a numeric \eqn{n \times n} matrix) and,
#'   when \code{reorder = TRUE}, \code{$consensus} (an integer vector of
#'   cluster labels).
#' @param reorder Logical. If \code{TRUE} (default), rows and columns are
#'   reordered by the consensus cluster assignment stored in
#'   \code{consensus_result$consensus}, making the block structure of the
#'   matrix visually apparent.
#' @param show_labels Logical. If \code{TRUE}, observation indices are
#'   displayed on both axes. Recommended only for small datasets as labels
#'   become unreadable for large \eqn{n}. Default is \code{FALSE}.
#' @param viridis_option A single character specifying the color palette
#'   passed to \code{\link[viridisLite]{viridis}}. Accepted values are
#'   \code{"A"} (magma), \code{"B"} (inferno), \code{"C"} (plasma), and
#'   \code{"D"} (viridis, default). Low co-assignment values are mapped to
#'   dark colors and high values to bright colors.
#'
#' @return Invisibly returns \code{NULL}. The function is called for its
#'   side effect of producing a plot.
#'
#' @details
#' The matrix is rendered using \code{\link[graphics]{image}} on the
#' consensus dissimilarity \eqn{1 - C} reordered by cluster membership.
#' Original graphical parameters are restored on exit via
#' \code{\link[graphics]{par}}.
#'
#' A well-defined consensus solution produces a block-diagonal pattern in
#' the heatmap, where bright blocks along the diagonal correspond to
#' observations that are consistently assigned to the same cluster, and
#' dark off-diagonal regions indicate observations that are rarely grouped
#' together.
#'
#' @examples
#' \donttest{
#' library(mice)
#'
#' imp   <- mice(nhanes, m = 5, printFlag = FALSE)
#' mild  <- as_mild_list(imp)
#' parts <- cluster_imputations(mild, method = "ward.D2", k = 3)
#' cons  <- consensus_clustering(parts, k = 3)
#'
#' # Default plot with viridis palette
#' plot_consensus_matrix(cons)
#'
#' # Without reordering
#' plot_consensus_matrix(cons, reorder = FALSE)
#'
#' # Compare color palettes
#' par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
#' for (opt in c("A", "B", "C", "D")) {
#'   plot_consensus_matrix(cons, viridis_option = opt)
#' }
#' par(mfrow = c(1, 1))
#' }
#'
#' @seealso
#' \code{\link{consensus_clustering}}, \code{\link{plot_consensus_dendrogram}},
#' \code{\link{plot_validation_metrics}}
#'
#' @export
plot_consensus_matrix <- function(consensus_result,
                                  reorder        = TRUE,
                                  show_labels    = FALSE,
                                  viridis_option = "D") {

  if (is.null(consensus_result$coassignment)) {
    stop("consensus_result must contain 'coassignment'.")
  }

  M <- consensus_result$coassignment

  if (!is.matrix(M)) {
    stop("'coassignment' must be a matrix.")
  }

  ord <- seq_len(nrow(M))

  if (reorder) {
    if (!is.null(consensus_result$consensus)) {
      ord <- order(consensus_result$consensus)
    }
  }

  M_ord <- M[ord, ord, drop = FALSE]

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))

  graphics::par(mar = c(4, 4, 2, 2))

  graphics::image(1:nrow(M_ord), 1:ncol(M_ord), t(M_ord[nrow(M_ord):1, ]),
        col   = viridisLite::viridis(100, option = viridis_option),
        axes  = FALSE,
        xlab  = "Observaciones",
        ylab  = "Observaciones",
        main  = paste("Matriz de consenso (k =", consensus_result$k, ")"))

  if (show_labels) {
    graphics::axis(1, at = seq_len(nrow(M_ord)), labels = ord,      las = 2, cex.axis = 0.4)
    graphics::axis(2, at = seq_len(ncol(M_ord)), labels = rev(ord), las = 2, cex.axis = 0.4)
  } else {
    graphics::axis(1)
    graphics::axis(2)
  }

  graphics::box()
}
