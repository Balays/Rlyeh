library(data.table)

# Inverse of sum() aggregation: repeat each row 'count' times
unsummarize <- function(DT, count,
                        cols = setdiff(names(DT), deparse(substitute(count))),
                        drop_zero = TRUE, validate = TRUE, keep_count = FALSE) {
  stopifnot(inherits(DT, "data.table"))

  # Accept both unquoted column name and "count" string
  cnt <- if (is.character(count)) count else deparse(substitute(count))

  if (validate) {
    if (!cnt %in% names(DT)) stop("Count column not found: ", cnt)
    if (anyNA(DT[[cnt]]))    stop("Count column contains NA")
    if (any(DT[[cnt]] < 0))  stop("Count column contains negatives")
    if (any(DT[[cnt]] %% 1 != 0)) stop("Count column must be integer-like")
  }

  DT2 <- if (drop_zero) DT[get(cnt) > 0] else copy(DT)
  if (!nrow(DT2)) return(DT[, ..cols][0])

  # Repeat row indices and subset once (fast)
  idx <- rep.int(DT2[, .I], times = as.integer(DT2[[cnt]]))
  cols_to_keep <- if (keep_count) c(cols, cnt) else cols
  out <- DT2[idx, ..cols_to_keep]
  setcolorder(out, cols_to_keep)
  out
}
