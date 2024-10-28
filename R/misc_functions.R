
#' @export
#'
#'
rownames.from.1st.col <- function(x) {
  x <- data.frame(x[,-1], row.names = x[,1])
  return(x)
}


#' @export
#'
#'
luniq <- function(x) length(unique(x))

#' @export
#'
#'
dup   <- function(x) x[duplicated(x)]


