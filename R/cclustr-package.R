#' cclustr: Consensus Clustering Methods for Multiple Imputed Data
#'
#' @description
#' The \pkg{cclustr} package provides a comprehensive framework for performing
#' clustering analysis on datasets treated with multiple imputation. It
#' addresses the uncertainty inherent in missing data by finding a stable
#' consensus partition across multiple imputed realizations.
#'
#' @details
#' Standard statistical workflows for imputed data often focus on supervised
#' learning or univariate statistics. \pkg{cclustr} extends this to
#' unsupervised learning by implementing a co-assignment matrix approach.
#'
#' @section Main Workflow:
#' To achieve a robust consensus solution, follow these steps:
#' \enumerate{
#'   \item \strong{Standardization:} Use \code{\link{as_mild_list}} to format
#'   imputed data (e.g., from \pkg{mice}) into a consistent structure.
#'   \item \strong{Clustering:} Apply \code{\link{cluster_imputations}} to run
#'   algorithms (K-means, PAM, K-prototypes, etc.) across all imputations.
#'   \item \strong{Consensus:} Generate the final partition using
#'   \code{\link{consensus_clustering}} based on co-assignment probabilities.
#'   \item \strong{Selection:} Evaluate model quality and optimal cluster
#'   counts with \code{\link{choose_best_clustering}}.
#' }
#'
#' @section Supported Data Types:
#' \itemize{
#'   \item \strong{Numeric:} Hierarchical clustering, K-means, Fuzzy C-means, and Mclust.
#'   \item \strong{Categorical:} K-modes.
#'   \item \strong{Mixed:} PAM with Gower distance or K-prototypes.
#' }
#'
#' @author
#' \strong{Maintainer}: Anhuar Duran Mendoza \email{aduranm@@unbosque.edu.co}
#'
#' Authors:
#' \itemize{
#'   \item Andres Montenegro Lemus \email{afmontenegro@@unbosque.edu.co}
#'   \item Mario Pacheco Lopez \email{mpachecol@@unbosque.edu.co}
#' }
#'
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., & Golub, T. (2003).
#' Consensus clustering: A resampling-based method for class discovery
#' and visualization of gene expression microarray data.
#' \emph{Machine Learning}, \strong{52}(1-2), 91-118.
#' \doi{10.1023/A:1023949509487}
#'
#' Senbabaoglu, Y., Michailidis, G., & Li, J. Z. (2014).
#' Critical limitations of consensus clustering in class discovery.
#' \emph{Scientific Reports}, \strong{4}, 6207.
#' \doi{10.1038/srep06207}
#'
#' Pihur, V., Datta, S., & Datta, S. (2007).
#' Weighted rank aggregation of cluster validation measures.
#' \emph{Bioinformatics}, \strong{23}(13), 1607-1615.
#' \doi{10.1093/bioinformatics/btm158}
#'
#' Rubin, D. B. (1987).
#' \emph{Multiple Imputation for Nonresponse in Surveys}. Wiley.
#'
#' van Buuren, S., & Groothuis-Oudshoorn, K. (2011).
#' mice: Multivariate Imputation by Chained Equations in R.
#' \emph{Journal of Statistical Software}, \strong{45}(3), 1-67.
#' \doi{10.18637/jss.v045.i03}
#'
#' Hubert, L., & Arabie, P. (1985).
#' Comparing partitions.
#' \emph{Journal of Classification}, \strong{2}, 193-218.
#' \doi{10.1007/BF01908075}
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
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/andrews06ml/cclustr}
#'   \item Report bugs at \url{https://github.com/andrews06ml/cclustr/issues}
#' }
#'
#' @docType package
#' @keywords internal
"_PACKAGE"
