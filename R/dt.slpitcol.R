

### Split column into new columns based on a separator, if the number of new columns is unknown using data.table

library(data.table)

dt.splitcol <- function(DT, coltosplit, newcolpattern, sep) {
  # Ensure DT is a data.table
  setDT(DT)
  
  # Split the specified column by the separator
  split_results <- DT[, tstrsplit(get(coltosplit), sep, fixed = TRUE)]
  
  # Determine the number of new columns
  num_columns <- ncol(split_results)
  
  # Create new column names based on the provided pattern
  colnames(split_results) <- paste0(newcolpattern, 1:num_columns)
  
    DT <- cbind(DT, split_results)
  
  # Optionally, return the modified DT if desired
  return(DT)
}

