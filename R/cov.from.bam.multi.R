
## this imports the headers of a fast file quickly
get_fasta_headers <- function(fasta.ref) {
  
  headers <- system(paste0("grep '^>' ", fasta.ref), intern = TRUE)
  return(headers)
  
}

## this imports one sequence (according to the header argument) from a fasta file quickly
extract_sequence_by_header <- function(header, fasta.ref) {
  message('importing ', header, ' from ', fasta.ref, '...')
  cmd <- paste0("awk '/^", header, "/ {print; while(getline && !/^>/) print}' ", fasta.ref)
  sequence <- system(cmd, intern = TRUE)
  return(sequence)
  
}

## this imports all sequences from a fasta file quickly
import_seq_multi  <- function(fasta.ref, headers=NULL, nproc=24) {
  
  require(future)
  require(future.apply)
  plan(multicore, workers = nproc)
  
  if(is.null(headers)) {
    headers   <- get_fasta_headers(fasta.ref)
  }
  sequences <- purrr::map(headers, extract_sequence_by_header, fasta.ref)
  names(sequences) <- headers
  
  return(sequences)
  
}

## this counts the sequence lengths f
get_seq_lengths <- function(sequences, nproc=24) {
  
  require(future)
  require(future.apply)
  plan(multicore, workers = nproc)
  
  
  # Compute the sequence lengths (excluding headers)
  sequence_lengths <- future_lapply(sequences, function(seq) {
    sum(nchar(seq[-1]))
  })
  
  # Convert to a data frame
  df <- data.frame(
    seqnames = names(sequence_lengths),
    length   = sequence_lengths,
    stringsAsFactors = FALSE
  )
  
  return(df)
}

cov.from.bam <- function(bam, 
                         fasta.ref, 
                         intron=F, 
                         by=c('seqnames', 'start', 'end'), 
                         samples=NA,
                         nproc=24 ## use 1 for windows
                         ) {
  ## INTRONS ARE NOT CONSIDERED, IF THEY ARE NOT IN THE DATA FRAME !!
  
  #####
  fasta  <- seqinr::read.fasta(fasta.ref)
  
  seqlen.dt <- data.table(
    seqnames = names(fasta),
    length = sapply(fasta, seqinr::getLength)
  )
  
  seqnames <- unique(dt[,seqnames])
  seqname <- seqnames[1]
  
  
  j <- 36
  sample    <- samples[j]
  bam.samp  <- bam.filt[bam.filt$sample == sample, ]
   
  get.cov <- function(seqname, seqlen.dt, bam.cov) {
  
    dt <- seqlen.dt[seqlen.dt$seqnames == seqname, ]
    
    long_format_dt <- dt[, .(position = seq(1, length)), by = seqnames]
    
    # Add start and end columns
    long_format_dt[, c("start", "end") := .(position, position)]
    
    # Create + and - strand tables and rbind them together
    genome_dt <- rbind(copy(long_format_dt)[, strand := '+'], 
                       copy(long_format_dt)[, strand := '-'], 
                       use.names=TRUE)
    
    
    # Set the order of columns
    setcolorder(genome_dt, c('seqnames', 'start', 'end', 'strand'))
    
     count.df  <- data.frame(genome.df, count=GenomicRanges::countOverlaps(GRanges(genome.df), GRanges(bam.samp)), sample=sample)
      
  }
  
  purrr::map(seqnames[1], )
  
  
  
  library(data.table)
  
  # Assuming your data is in a data.table called dt
  setDT(dt)
  
  # For each row in the original data, create a list of positions that the read covers
  dt_expanded <- dt[, .(sample = sample, 
                        seqnames = seqnames, 
                        position = list(start:end), 
                        strand = strand), 
                    by = 1:nrow(dt)]
  
  # Convert the list of positions into long format
  dt_expanded <- dt_expanded[, .(sample, seqnames, position = unlist(position), strand)]
  
  # Group by the relevant columns and count the number of reads
  dt_summary <- dt_expanded[, .N, by = .(sample, seqnames, position, strand)]
  
  
  
  
  
  n_iter    <- nrow(dt)
  # Create a list of data tables
  dt_list <- lapply(1:n_iter, function(i) {
    setTxtProgressBar(pb, i) 
    data.table(
      seqnames = dt$seqnames[i],
      position = seq_len(dt$length[i])
    )
  })
  
  # Combine into one data.table
  long_format_dt <- rbindlist(dt_list)
  
  # Check the result
  head(long_format_dt)
  
  
  
  
  n_iter    <- nrow(dt)
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  # Loop through each sequence in dt
  long_format_dt <- data.frame(NULL)
  for(i in 1:n_iter) {
    setTxtProgressBar(pb, i, title = sample)
    temp_dt <- data.table(
      seqnames = dt$seqnames[i],
      position = seq_len(dt$length[i])
    )
    long_format_dt <- rbindlist(list(long_format_dt, temp_dt))
  }
  
  
  
  
  ##
  long_format_df$start <- long_format_df$position
  long_format_df$end   <- long_format_df$position
  
  genome.df <- plyr::rbind.fill(data.frame(long_format_df, strand='+'),
                                data.frame(long_format_df, strand='-')) [,c(by, 'strand')]
    
  
   
  
  
  ## filter for samples
  if (all(!is.na(samples))) {
    bam <- bam[is.element(bam$sample, samples), ]
    
    
    n_iter    <- length(samples)
    
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = n_iter, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    
    count.all <- data.frame(NULL)
    #i <- 2
    for (i in 1:n_iter){
      setTxtProgressBar(pb, i, title = sample)
      sample    <- samples[i]
      bam.samp  <- bam[bam$sample == sample, ]
      count.df  <- data.frame(genome.df, count=GenomicRanges::countOverlaps(GRanges(genome.df), GRanges(bam.samp)), sample=sample)
      count.all <- plyr::rbind.fill(count.df, count.all)
    }
    
  } else {
    count.all  <- data.frame(genome.df, count=GenomicRanges::countOverlaps(GRanges(genome.df), GRanges(bam)), sample=NA)
  }
  
  return(count.all)
}


##