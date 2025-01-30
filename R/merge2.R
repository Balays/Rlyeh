library(data.table)

merge2 <- function(x, y, by, ncol_y, ...) {
  # Convert data frames to data.tables
  x <- as.data.table(x)
  y <- as.data.table(y)

  # Perform the merge
  merged_dt <- merge(x, y, by = by, ...)

  # Get the names of the columns in y that were not used for merging
  y_cols <- setdiff(names(y), by)

  # Get the names of the columns to be moved
  cols_to_move <- y_cols[ncol_y]

  # Reorder the columns in the merged data.table
  new_order <- c(names(x), cols_to_move, setdiff(names(merged_dt), c(names(x), cols_to_move)))
  setcolorder(merged_dt, new_order)

  return(merged_dt)
}

# Example usage
x <- data.table(id = 1:5, x1 = letters[1:5])
y <- data.table(id = 1:5, y1 = LETTERS[1:5], y2 = 6:10, y3 = 11:15)

# Merge and move columns y1 and y2 to the specified positions
result <- merge2(x, y, by = "id", ncol_y = c(2, 3))
print(result)
