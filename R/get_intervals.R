library(data.table)

#' Generate equally spaced intervals between start and end
#'
#' @param start Numeric scalar. Start of the range.
#' @param end Numeric scalar. End of the range.
#' @param n Integer. Number of equally‑spaced chunks desired.
#'
#' @return A data.table with two columns: `from` and `to`, each row representing one interval.
#' @examples
#' get_intervals(0, 10, 4)
#' #    from   to
#' # 1:  0.0  2.5
#' # 2:  2.5  5.0
#' # 3:  5.0  7.5
#' # 4:  7.5 10.0
get_intervals <- function(start, end, n, round=1) {
  if (!is.numeric(start) || !is.numeric(end) || !is.numeric(n)) {
    stop("start, end, and n must be numeric.")
  }
  if (n <= 0 || n %% 1 != 0) {
    stop("n must be a positive integer.")
  }
  if (start == end) stop("start and end must be different values.")
  
  reversed <- FALSE
  if (start > end) {
    temp <- start
    start <- end
    end <- temp
    reversed <- TRUE
  }
  
  breaks <- seq(from = start, to = end, length.out = n + 1)
  
  dt <- data.table(
    from = breaks[-length(breaks)],
    to   = breaks[-1]
  )
  
  if(is.na(round) | is.null(round) ) { NULL } else {
    dt[, from := round(from, round)]
    dt[, to   := round(to,   round)]
  }
  
  if (reversed) {
    dt <- dt[order(-from)]
  }
  
  dt[]
}


#get_intervals(visfrom, visto, n_breaks, 0)

#10^round(log10(visto - visfrom), 0) / n_breaks
