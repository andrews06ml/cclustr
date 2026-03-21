#' Plot validation metrics across candidate numbers of clusters
#'
#' @description
#' Produces a multi-panel line plot displaying one or more clustering
#' validation metrics as a function of the number of clusters \code{k},
#' as computed by \code{\link{validate_clustering}} and optionally extended
#' with a \code{score} column by \code{\link{choose_best_clustering}}. For
#' each metric, the optimal value is highlighted in red.
#'
#' @param validation_table A \code{data.frame} as returned by
#'   \code{\link{validate_clustering}} or \code{$scores_table} from
#'   \code{\link{choose_best_clustering}}. Must contain a column \code{k}
#'   and at least one of the requested metric columns.
#' @param metrics A character vector specifying which metrics to plot.
#'   Only metrics present in \code{validation_table} are plotted; others
#'   are silently ignored. Default includes all standard metrics plus
#'   \code{score}:
#'   \code{"pac"}, \code{"silhouette_mean"},
#'   \code{"ari_mean_between_imputations"}, \code{"ari_consensus_mean"},
#'   \code{"calinski_harabasz_mean"}, \code{"davies_bouldin_mean"},
#'   \code{"dunn_index"}, and \code{"score"}.
#' @param ask Logical. If \code{TRUE}, the user is prompted before each new
#'   plot is drawn, useful when displaying many metrics interactively.
#'   Default is \code{FALSE}.
#'
#' @return Invisibly returns \code{NULL}. The function is called for its
#'   side effect of producing a plot.
#'
#' @details
#' The number of panels is determined automatically from the number of
#' available metrics using \code{ceiling(sqrt(n))} rows and
#' \code{ceiling(n / nrow)} columns, producing a layout as square as
#' possible.
#'
#' The direction of optimality is accounted for when highlighting the best
#' value: lower is better for \code{pac}, \code{davies_bouldin_mean}, and
#' \code{score}; higher is better for all remaining metrics. The optimal
#' point for each metric is highlighted in red.
#'
#' Original graphical parameters are restored on exit via
#' \code{\link[graphics]{par}}.
#'
#' When \code{score} is present (i.e., \code{validation_table} comes from
#' \code{\link{choose_best_clustering}}), the panel for \code{score} shows
#' the overall weighted rank aggregation - the \code{k} with the lowest
#' score is the recommended solution.
#'
#' @examples
#' \dontrun{
#' library(mice)
#'
#' imp    <- mice(nhanes, m = 5, printFlag = FALSE)
#' mild   <- as_mild_list(imp)
#' parts  <- cluster_imputations(mild, method = "ward.D2", k = 2:5)
#' cons   <- consensus_clustering(parts)
#' val    <- validate_clustering(parts, cons)
#'
#' # Plot all available metrics (no score column yet)
#' plot_validation_metrics(val)
#'
#' # Include score column from choose_best_clustering
#' best <- choose_best_clustering(val, cons)
#' plot_validation_metrics(best$scores_table)
#'
#' # Plot only stability metrics
#' plot_validation_metrics(val,
#'                         metrics = c("pac",
#'                                     "ari_mean_between_imputations",
#'                                     "ari_consensus_mean"))
#' }
#'
#' @seealso
#' \code{\link{validate_clustering}}, \code{\link{choose_best_clustering}},
#' \code{\link{plot_consensus_dendrogram}}, \code{\link{plot_consensus_matrix}}
#'
#' @export
plot_validation_metrics <- function(validation_table,
                                    metrics = c("pac",
                                                "silhouette_mean",
                                                "ari_mean_between_imputations",
                                                "ari_consensus_mean",
                                                "calinski_harabasz_mean",
                                                "davies_bouldin_mean",
                                                "dunn_index",
                                                "score"),
                                    ask = FALSE) {

  if (!is.data.frame(validation_table)) {
    stop("validation_table must be a data.frame.")
  }

  if (!("k" %in% names(validation_table))) {
    stop("validation_table must contain column 'k'.")
  }

  available_metrics <- intersect(metrics, names(validation_table))

  if (length(available_metrics) == 0) {
    stop("None of the requested metrics are available in validation_table.")
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  n         <- length(available_metrics)
  nrow_plot <- ceiling(sqrt(n))
  ncol_plot <- ceiling(n / nrow_plot)

  par(mfrow = c(nrow_plot, ncol_plot), ask = ask)

  for (metric in available_metrics) {

    y <- validation_table[[metric]]

    plot(validation_table$k, y,
         type = "b",
         pch  = 19,
         xlab = "Number of clusters (k)",
         ylab = metric,
         main = gsub("_", " ", metric))

    best_idx <- which.max(
      if (metric %in% c("pac", "davies_bouldin_mean", "score")) -y else y
    )

    points(validation_table$k[best_idx], y[best_idx],
           pch = 19, cex = 1.4, col = "red")  # <- corregido
  }
}

