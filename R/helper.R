#' Internal Helper Function for Printing Matrices
#'
#' Internal helper function for printing matrices with specified digits.
#'
#' @param x A matrix or numeric object to print.
#' @param digits Number of digits for rounding. Default is 4.
#'
#' @return No return value, called for side effects
#'
#' @keywords internal
#' @noRd
#.mprint <- function(x, digits=4) print(round(x, digits), digits=digits)
.mprint <- function(x, digits=4) print(noquote(format(round(x, digits), digits=digits)))
