# R/Imputation.R
utils::globalVariables(c(".imp"))

#' Standardize multiple imputation outputs into a unified list
#'
#' @description
#' Converts multiple imputation outputs from different packages into a
#' standardized list of completed datasets, where each element corresponds
#' to one imputed dataset. Supported formats include \code{mids} objects
#' from \pkg{mice}, long-format data frames, \code{amelia} objects,
#' \code{imputationList} objects from \pkg{mitools}, and plain lists of
#' data frames. After conversion, consistency and data quality checks are
#' applied across all imputed datasets.
#'
#' This function is designed as a preprocessing step for clustering
#' multiple imputed datasets, ensuring strict comparability across imputations.
#'
#' @param x An imputation object. Accepted classes and formats:
#' \itemize{
#'   \item \code{mids}: output of \code{mice::mice()}.
#'   \item \code{data.frame}: long-format output of
#'         \code{mice::complete(x, action = "long")}. Must contain a
#'         \code{.imp} column.
#'   \item \code{amelia}: output of \code{Amelia::amelia()}.
#'   \item \code{imputationList}: output of \code{mitools::imputationList()}.
#'   \item \code{list}: a plain list where every element is a
#'         \code{data.frame} representing one completed dataset.
#' }
#'
#' @return A named list of \code{data.frame} objects, one per imputed
#' dataset. All datasets are guaranteed to have identical dimensions and
#' column names, and to contain no \code{NA}, \code{NaN}, \code{Inf}, or
#' all equal values in a column.
#'
#' @details
#' The function performs the following consistency checks before returning:
#' \itemize{
#'   \item Identical number of rows across all imputed datasets.
#'   \item Identical number of columns across all imputed datasets.
#'   \item Identical column names across all imputed datasets.
#'   \item Absence of \code{NA}, \code{NaN}, and \code{Inf} values.
#'   \item Absence of numeric columns with all values equal,
#'         as such variables provide no information for distance-based
#'         clustering and may lead to degenerate solutions.
#' }
#' All datasets are coerced to base \code{data.frame} to avoid
#' compatibility issues with \code{tibble} or other subclasses.
#'
#' @examples
#' # ------------------------------------------------------------
#' # Example 1: Simple list of completed datasets
#' # ------------------------------------------------------------
#' set.seed(123)
#'
#' imp_list <- replicate(3, {
#'   data.frame(
#'     x = rnorm(10),
#'     y = rnorm(10)
#'   )
#' }, simplify = FALSE)
#'
#' mild <- as_mild_list(imp_list)
#'
#' length(mild)      # number of imputations
#' str(mild[[1]])    # structure of one dataset
#'
#' \donttest{
#' # ------------------------------------------------------------
#' # Example 2: List of completed datasets with mice
#' # ------------------------------------------------------------
#'
#' if (requireNamespace("mice", quietly = TRUE)) {
#'
#'   set.seed(123)
#'
#'   # Example 1: Using a mids object
#'   imp <- mice::mice(mice::nhanes, m = 5, printFlag = FALSE)
#'   mild <- as_mild_list(imp)
#'
#'   length(mild)        # number of imputations
#'   str(mild[[1]])      # structure of one completed dataset
#'
#'   # Example 2: Using long-format data
#'   long_df <- mice::complete(imp, action = "long", include = FALSE)
#'   mild2 <- as_mild_list(long_df)
#'
#'   length(mild2)
#' }
#' }
#'
#' @seealso \code{\link{cluster_imputations}}, \code{\link{consensus_clustering}}
#'
#' @export
as_mild_list <- function(x) {

  #Control verbosity
  verbose <- getOption("cclustr.verbose", FALSE)

  # --------------------------------------------------------
  # Case 1: 'mids' object (output of mice())
  # --------------------------------------------------------
  if (inherits(x, "mids")) {
    if (verbose) message("Detected: mice (mids object)")
    imp_list <- mice::complete(x, action = "all")
  }

  # --------------------------------------------------------
  # Case 2: mice Long format (output of mice::complete(..., "long"))
  # --------------------------------------------------------
  # Each imputation is identified by the '.imp' column
  else if (is.data.frame(x)) {

    if (!(".imp" %in% names(x))) {
      stop("The data.frame does not contain the '.imp' column.")
    }

    if (verbose) message("Detected: long-format imputations")

    # Identify imputation indices (excluding original incomplete data: .imp = 0)
    imp_ids <- sort(unique(x$.imp))
    imp_ids <- imp_ids[imp_ids != 0]  # delete original data

    # Split the long data.frame into a list of completed datasets
    imp_list <- lapply(imp_ids, function(i) {
      subset(x, .imp == i)[ , !(names(x) %in% c(".imp", ".id")), drop = FALSE]
    })

    # Assign standardized names to each imputed dataset
    names(imp_list) <- paste0("imp", imp_ids)
  }

  # --------------------------------------------------------
  # Case 3: 'amelia' object (output of amelia())
  # --------------------------------------------------------
  else if (inherits(x, "amelia")) {
    if (verbose) message("Detected: Amelia object")
    imp_list <- x$imputations
  }

  #-------------------------------------------------------
  # Case 4: 'imputationList' object (output of mitools())
  #-------------------------------------------------------
  else if (inherits(x, "imputationList")) {
    if (verbose) message("Detected: mitools imputationList")
    imp_list <- x$imputations
  }

  # --------------------------------------------------------
  # Case 5: already a list of completed datasets
  # --------------------------------------------------------
  # Each element must be a completed dataset
  else if (is.list(x)) {

    if (!all(sapply(x, is.data.frame))) {
      stop("The input list must contain only data.frames.")
    }

    imp_list <- x
  }

  #-------------------------------------------------------
  # Case 6: Unsupported input type
  #-------------------------------------------------------
  else {
    stop(
      paste(
        "Unsupported imputation format.\n\n",
        "Supported inputs are:\n",
        " - mids object (mice)\n",
        " - long-format data.frame from mice::complete(..., 'long')\n",
        " - amelia object (Amelia package)\n",
        " - imputationList object (mitools package)\n",
        " - list of completed data.frames\n\n"
        )
    )
  }

  #-------------------------------------------------------
  # Final consistency checks
  #-------------------------------------------------------
  # All imputations must have the same number of rows
  if (length(unique(sapply(imp_list, nrow))) != 1) {
    stop("All imputations must have the same number of rows.")
  }

  # All imputations must have the same number of columns
  if (length(unique(sapply(imp_list, ncol))) != 1) {
    stop("All imputations must have the same number of columns.")
  }

  # All imputations must have identical column names
  ref_names <- names(imp_list[[1]])
  if (!all(sapply(imp_list, function(df) identical(names(df), ref_names)))) {
    stop("All imputations must have identical column names.")
  }

  # Convert to data.frame (avoid tibble issues)
  imp_list <- lapply(imp_list, as.data.frame)

  #-------------------------------------------------------
  # check missing, infinite and equal values
  #-------------------------------------------------------

  for (i in seq_along(imp_list)) {

    df <- imp_list[[i]]

    # 1) Check NA values
    if (anyNA(df)) {
      stop(paste("Imputation", i, "contains NA values."))
    }

    # 2) Check NaN values
    if (any(is.nan(as.matrix(df)))) {
      stop(paste("Imputation", i, "contains NaN values."))
    }

    # 3) Check Inf values
    if (any(is.infinite(as.matrix(df)))) {
      stop(paste("Imputation", i, "contains Inf values."))
    }

    # 4) Check columns with the same values
    same_cols <- names(df)[
      sapply(names(df), function(nm) {
        col <- df[[nm]]
        is.numeric(col) && length(unique(col)) == 1
      })
    ]
    if (length(same_cols) > 0) {
      const_vals <- sapply(same_cols, function(nm) unique(df[[nm]]))
      details <- paste(same_cols, "=", const_vals, collapse = ", ")
      stop(paste("Imputation", i,
                 "has constant numeric columns (all values identical):",
                 details))
    }
  }

  if (verbose) {
    message("All datasets passed validation")
    message(paste("Created list of", length(imp_list), "imputed datasets"))
  }

  # Return the standardized list of imputed datasets
  if (is.null(names(imp_list))) {
    names(imp_list) <- paste0("imp", seq_along(imp_list))
  }

  return(imp_list)
}
