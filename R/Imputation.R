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
#' zero values.
#'
#' @details
#' The function performs the following consistency checks before returning:
#' \itemize{
#'   \item Identical number of rows across all imputed datasets.
#'   \item Identical number of columns across all imputed datasets.
#'   \item Identical column names across all imputed datasets.
#'   \item Absence of \code{NA}, \code{NaN}, and \code{Inf} values.
#'   \item Absence of zero values (strict validation for downstream
#'         clustering).
#' }
#' All datasets are coerced to base \code{data.frame} to avoid
#' compatibility issues with \code{tibble} or other subclasses.
#'
#' @examples
#' \dontrun{
#' library(mice)
#'
#' # Example with a mids object
#' imp <- mice(nhanes, m = 5, printFlag = FALSE)
#' mild <- as_mild_list(imp)
#' length(mild)  # 5
#'
#' # Example with long-format data frame
#' long_df <- mice::complete(imp, action = "long", include = FALSE)
#' mild2 <- as_mild_list(long_df)
#' }
#'
#' @seealso \code{\link{cluster_imputations}}, \code{\link{consensus_clustering}}
#'
#' @export
as_mild_list <- function(x) {

  # --------------------------------------------------------
  # Case 1: 'mids' object (output of mice())
  # --------------------------------------------------------
  if (inherits(x, "mids")) {
    message("Detected: mice (mids object)")
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

    message("Detected: long-format imputations")

    # Identify imputation indices (excluding original incomplete data: .imp = 0)
    imp_ids <- sort(unique(x$.imp))
    imp_ids <- imp_ids[imp_ids != 0]  # delete original data

    # Split the long data.frame into a list of completed datasets
    imp_list <- lapply(imp_ids, function(i) {
      subset(x, .imp == i)[ , !(names(x) %in% c(".imp", ".id"))]
    })

    # Assign standardized names to each imputed dataset
    names(imp_list) <- paste0("imp", imp_ids)
  }

  # --------------------------------------------------------
  # Case 3: 'amelia' object (output of amelia())
  # --------------------------------------------------------
  else if (inherits(x, "amelia")) {
    message("Detected: Amelia object")
    imp_list <- x$imputations
  }

  #-------------------------------------------------------
  # Case 4: 'imputationList' object (output of mitools())
  #-------------------------------------------------------
  else if (inherits(x, "imputationList")) {
    message("Detected: mitools imputationList")
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
        "• mids object (mice)\n",
        "• long-format data.frame from mice::complete(..., 'long')\n",
        "• amelia object (Amelia package)\n",
        "• imputationList object (mitools package)\n",
        "• list of completed data.frames\n\n"
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
  if (!all(sapply(imp_list, function(df)
    identical(names(df), names(imp_list[[1]]))))) {
    stop("All imputations must have identical column names.")
  }

  # Convert to base data.frame (avoid tibble issues)
  imp_list <- lapply(imp_list, as.data.frame)

  #-------------------------------------------------------
  # check missing, infinite and zero values
  #-------------------------------------------------------

  for (i in seq_along(imp_list)) {

    df <- imp_list[[i]]

    # 1) Check NA values
    if (anyNA(df)) {
      stop(paste("Imputation", i, "contains NA values."))
    }

    # 2) Check NaN values
    if (any(sapply(df, function(col) any(is.nan(col))))) {
      stop(paste("Imputation", i, "contains NaN values."))
    }

    # 3) Check Inf values
    if (any(sapply(df, function(col) any(is.infinite(col))))) {
      stop(paste("Imputation", i, "contains Inf values."))
    }

    # 4) Check columns fully equal to zero — solo columnas numéricas
    zero_cols <- names(df)[
      sapply(names(df), function(nm) {
        col <- df[[nm]]
        is.numeric(col) && all(col == 0)
      })
    ]
    if (length(zero_cols) > 0) {
      stop(paste("Imputation", i, "has columns with only zeros:",
                 paste(zero_cols, collapse = ", ")))
    }
  }

  message("✓ All datasets passed validation")
  message(paste("✓ Created list of", length(imp_list), "imputed datasets"))

  # Return the standardized list of imputed datasets
  return(imp_list)
}
