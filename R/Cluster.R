
#' Perform clustering on multiple imputed datasets
#'
#' @description
#' Applies a clustering algorithm to each completed dataset in a standardized
#' imputation list (as produced by \code{\link{as_mild_list}}). Supports
#' hierarchical, k-means, PAM, fuzzy c-means, and model-based clustering.
#' Accepts a single value of \code{k} or a range, and optionally scales the
#' data prior to clustering.
#'
#' @param imp_list A named list of \code{data.frame} objects, as returned by
#'   \code{\link{as_mild_list}}. All datasets must have identical dimensions
#'   and column names.
#' @param method A character string specifying the clustering algorithm.
#'   For hierarchical clustering, accepted values are \code{"ward.D"},
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"centroid"}, \code{"median"}, and \code{"mcquitty"}.
#'   Additional options are \code{"kmeans"}, \code{"pam"}, \code{"fuzzy"},
#'   and \code{"mclust"}. Default is \code{"ward.D2"}.
#' @param k A single integer or an integer vector specifying the number of
#'   clusters. If a vector is provided, clustering is performed for each
#'   value of \code{k}.
#' @param scale_data Logical. If \code{TRUE} (default), numeric columns are
#'   standardized to zero mean and unit variance prior to clustering.
#' @param ... Additional arguments passed to the underlying clustering
#'   function (\code{kmeans}, \code{cluster::pam}, \code{e1071::cmeans},
#'   or \code{mclust::Mclust}). For \code{method = "pam"}, passing
#'   \code{metric = "gower"} enables Gower distance for mixed data types.
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
  
  # --------------------------------------------------------
  # Supported methods
  # --------------------------------------------------------
  hierarchical_methods <- c(
    "ward.D", "ward.D2", "single", "complete",
    "average", "centroid", "median", "mcquitty"
  )
  
  if (!(method %in% c(hierarchical_methods, "kmeans", "pam", "fuzzy", "mclust"))) {
    stop("Unsupported clustering method.")
  }
  
  # --------------------------------------------------------
  # Optional escalation
  # --------------------------------------------------------
  database <- if (scale_data) {
    lapply(imp_list, function(df) {
      num_cols <- sapply(df, is.numeric)
      df[, num_cols] <- scale(df[, num_cols])
      df
    })
  } else {
    imp_list
  }
  
  # -------------------------------------------------------
  # Mixed data validation
  # --------------------------------------------------------
  has_factors <- any(sapply(imp_list[[1]], function(col)
    is.factor(col) || is.character(col)))
  
  if (has_factors && !(method == "pam")) {
    stop(
      "The dataset contains categorical variables. ",
      "Use method = 'pam' with metric = 'gower' for mixed data types."
    )
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
        d <- dist(df)
        hc <- hclust(d, method = method)
        cutree(hc, k = k)
        
      } else if (method == "kmeans") {
        kmeans(df, centers = k, ...)$cluster
        
      } else if (method == "pam") {
        
        extra_args <- list(...)
        
        if (!is.null(extra_args$metric) && extra_args$metric == "gower") {
          
          # Gower distance for mixed data
          d <- cluster::daisy(df, metric = "gower")
          cluster::pam(d, k, diss = TRUE)$clustering
          
        } else {
          
          # Clasic pam (numeric)
          cluster::pam(df, k, ...)$clustering
        }
        
      } else if (method == "fuzzy") {
        
        e1071::cmeans(df, centers = k, ...)$cluster
        
      } else if (method == "mclust") {
        
        mclust::Mclust(df, G = k, ...)$classification
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
        d <- dist(df)
        hclust(d, method = method)
      })
      
      partitions_by_k <- lapply(k, function(k_i) {
        
        partitions <- lapply(hc_list, function(hc) {
          cutree(hc, k = k_i)
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
            
            kmeans(df, centers = k_i, ...)$cluster
            
          } else if (method == "pam") {
            
            extra_args <- list(...)
            
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