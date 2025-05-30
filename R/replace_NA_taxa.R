replace_NA_taxa <- function(dt, pattern = "unknown_", mode = c("column", "left", "last")) {

  require(data.table)
  mode <- match.arg(mode)

  dt_copy <- copy(dt)
  cols <- colnames(dt_copy)

  if (mode == "column") {
    for (col in cols) {
      na_rows <- which(is.na(dt_copy[[col]]))
      set(dt_copy, i = na_rows, j = col, value = paste0(pattern, col))
    }
  }

  if (mode == "left") {
    for (i in seq_len(nrow(dt_copy))) {
      row_vals <- as.character(dt_copy[i, ..cols])
      for (j in seq_along(cols)) {
        if (is.na(row_vals[j])) {
          non_na_left <- if (j > 1) row_vals[1:(j-1)][!is.na(row_vals[1:(j-1)])] else character(0)
          replacement <- if (length(non_na_left)) tail(non_na_left, 1) else cols[j]
          set(dt_copy, i = i, j = cols[j], value = paste0(pattern, replacement))
        }
      }
    }
  }

  if (mode == "last") {
    for (i in seq_len(nrow(dt_copy))) {
      row_vals <- as.character(dt_copy[i, ..cols])
      non_na_vals <- row_vals[!is.na(row_vals)]
      last_non_na <- if (length(non_na_vals)) tail(non_na_vals, 1) else NULL
      replacement <- if (!is.null(last_non_na)) last_non_na else cols
      for (j in seq_along(cols)) {
        if (is.na(row_vals[j])) {
          set(dt_copy, i = i, j = cols[j], value = paste0(pattern, replacement))
        }
      }
    }
  }

  return(dt_copy)
}
