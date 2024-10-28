

foverlaps4 <- function(DTx, DTy,
                       # Define upstream and downstream distances
                       prime5.upstream = 10,
                       prime5.downstream = 10,
                       prime3.upstream = 10,
                       prime3.downstream = 10) {
  
  require(data.table)
  
  # Check if upstream and downstream distances are provided for prime5
  do_prime5 <- !is.null(prime5.upstream) && !is.null(prime5.downstream)
  do_prime3 <- !is.null(prime3.upstream) && !is.null(prime3.downstream)
  
  # Create copies to avoid modifying original data
  DTx_copy <- copy(DTx)
  DTy_copy <- copy(DTy)
  
  # Add unique identifiers if not present
  if (!'tx_id' %in% names(DTx_copy)) {
    DTx_copy[, tx_id := .I]
  }
  if (!'y_id' %in% names(DTy_copy)) {
    DTy_copy[, y_id := .I]
  }
  
  # Initialize list to hold overlap results
  overlap_results <- list()
  
  ## Overlaps based on 5-prime positions
  if (do_prime5) {
    message('Finding overlaps based on 5-prime positions...')
    
    # Create window start and end in DTy based on strand
    DTy_prime5 <- copy(DTy_copy)
    DTy_prime5[, `:=`(
      win_start = fifelse(strand == '+', prime5 - prime5.upstream, prime5 - prime5.downstream),
      win_end   = fifelse(strand == '+', prime5 + prime5.downstream, prime5 + prime5.upstream)
    )]
    
    # Ensure positions are not negative
    DTy_prime5[, `:=`(
      win_start = pmax(win_start, 0),
      win_end   = pmax(win_end, 0)
    )]
    
    # Set 'start' and 'end' in DTy
    DTy_prime5[, `:=`(start = win_start, end = win_end)]
    
    # Set 'start' and 'end' in DTx
    DTx_prime5 <- copy(DTx_copy)
    DTx_prime5[, `:=`(start = prime5, end = prime5)]
    
    # Set keys
    setkey(DTy_prime5, seqnames, strand, start, end)
    setkey(DTx_prime5, seqnames, strand, start, end)
    
    # Perform overlap
    DTprime5 <- foverlaps(DTx_prime5, DTy_prime5, nomatch = 0)
    DTprime5[, overlap_type := 'prime5']
    
    # Store result
    overlap_results$prime5 <- DTprime5
  } else {
    message('Skipping overlaps based on 5-prime positions.')
  }
  
  ## Overlaps based on 3-prime positions
  if (do_prime3) {
    message('Finding overlaps based on 3-prime positions...')
    
    # Create window start and end in DTy based on strand
    DTy_prime3 <- copy(DTy_copy)
    DTy_prime3[, `:=`(
      win_start = fifelse(strand == '+', prime3 - prime3.upstream, prime3 - prime3.downstream),
      win_end   = fifelse(strand == '+', prime3 + prime3.downstream, prime3 + prime3.upstream)
    )]
    
    # Ensure positions are not negative
    DTy_prime3[, `:=`(
      win_start = pmax(win_start, 0),
      win_end   = pmax(win_end, 0)
    )]
    
    # Set 'start' and 'end' in DTy
    DTy_prime3[, `:=`(start = win_start, end = win_end)]
    
    # Set 'start' and 'end' in DTx
    DTx_prime3 <- copy(DTx_copy)
    DTx_prime3[, `:=`(start = prime3, end = prime3)]
    
    # Set keys
    setkey(DTy_prime3, seqnames, strand, start, end)
    setkey(DTx_prime3, seqnames, strand, start, end)
    
    # Perform overlap
    DTprime3 <- foverlaps(DTx_prime3, DTy_prime3, nomatch = 0)
    DTprime3[, overlap_type := 'prime3']
    
    # Store result
    overlap_results$prime3 <- DTprime3
  } else {
    message('Skipping overlaps based on 3-prime positions.')
  }
  
  ## Merge results
  if (do_prime5 && do_prime3) {
    message('Merging overlaps based on both 5-prime and 3-prime positions...')
    # Merge the two results on 'tx_id' and 'y_id'
    setnames(overlap_results$prime5, old = c('start', 'end'), new = c('start_prime5', 'end_prime5'))
    setnames(overlap_results$prime3, old = c('start', 'end'), new = c('start_prime3', 'end_prime3'))
    
    # Select relevant columns
    cols_to_keep <- c('tx_id.i', 'y_id', 'seqnames', 'strand', 'overlap_type')
    DTprime5_reduced <- overlap_results$prime5[, ..cols_to_keep]
    DTprime3_reduced <- overlap_results$prime3[, ..cols_to_keep]
    
    # Merge on 'tx_id.i' and 'y_id'
    DTmerged <- merge(DTprime5_reduced, DTprime3_reduced, by = c('tx_id.i', 'y_id'))
    
    # Retrieve full details
    DTxy <- DTmerged[
      overlap_results$prime5,
      on = c('tx_id.i', 'y_id'),
      nomatch = 0
    ][
      overlap_results$prime3,
      on = c('tx_id.i', 'y_id'),
      nomatch = 0,
      suffixes = c('.prime5', '.prime3')
    ]
    
  } else if (do_prime5) {
    DTxy <- overlap_results$prime5
  } else if (do_prime3) {
    DTxy <- overlap_results$prime3
  } else {
    stop('No overlaps were performed. Please provide upstream and downstream distances.')
  }
  
  # Return the result
  return(DTxy)
}





