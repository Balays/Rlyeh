#library(misc)
library(GenomicRanges)


## The pileup function from the Rsamtools package generates coverage information for each position in the genome, 
## reporting coverage for each nucleotide type separately, that is
## it outputs the number of reads that align to that position with each type of nucleotide
require(Rsamtools)

cov.from.bam <- function(bamfile, pattern='.bam', pileup=F, param=ScanBamParam(what=scanBamWhat(), flag=scanBamFlag()), ...) {
  bam.cov <- data.table(NULL)
  try({
    #bamfile <- bamfiles[i]
    bamname <- gsub('.*\\/', '', bamfile)
    bamname <- gsub(pattern, '', bamname)
    message('Start analyzing ', bamfile, '...')
    bam.cov <- pileup(bamfile, scanBamParam=param, ...)
    setDT(bam.cov)
    
    # Summarize the counts for each position
    if( pileup == F) {
      total_coverage <- bam.cov[, .(count = sum(count)), by = .(seqnames, pos, strand)]
      bam.cov <- total_coverage
    }
    
    bam.cov[,sample := bamname]
  })
  
  #bam.cov <- bam.cov[order(sample, seqnames, qname, -mapq, flag, cigar, strand)]
  #bam.cov[ , aln_ID  := .GRP, by = .(sample, seqnames, qname, cigar, flag, mapq, strand)]
  #bam[ , aln_nr  := rank(aln_ID, ties.method="min"), by = qname]
  #if(length(unique(bam[,aln_ID])) != length(unique(bam[,aln_ID]))) { message('something is wrong with the IDs !') } 
  
  return(bam.cov)
}



### coverage from bam data frame

cov.from.bam.arch <- function(bam, fasta.ref, intron=F, by=c('seqnames', 'start', 'end'), samples=NA) {
  ## this is depracated, use cov.from.bam()
  
  ## INTRONS ARE NOT CONSIDERED, IF THEY ARE NOT IN THE DATA FRAME !!
  ## WORKS FOR THE FIRST CONTIG OF THE .FASTA file ONLY !!
  #gff::check_cnames(bam, c(by, 'sample'))
  fasta  <- seqinr::read.fasta(fasta.ref)
  l_genome <- length(fasta[[1]])
  #genome <- gsub('.fasta', fasta.ref)
  genome <- names(fasta)[1]
  genome.df <- plyr::rbind.fill(data.frame(seqnames = genome, start=1:l_genome, end=1:l_genome, strand='+'),
                                data.frame(seqnames = genome, start=1:l_genome, end=1:l_genome, strand='-'))

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


### Summarise coverages by windows

window.cov <- function(cov.df, fasta.ref=NA, seq_nr=1, window_size = 100, window_step = 100, by=c('seqnames', 'start', 'end'),
                       genome=NULL, l_genome=NULL) {

  if(!is.na(fasta.ref) & !is.na(seq_nr)) {
    fasta  <- seqinr::read.fasta(fasta.ref)
    l_genome <- length(fasta[[seq_nr]])
    #genome <- gsub('.fasta', fasta.ref)
    genome <- names(fasta)[seq_nr]
  } else {
    stop('either reference fasta file and nr of the sequence (contig) in it (if mulitifasta); or genome name and length must be provided!')
  }

  ### Get windows
  window_df           <- data.frame(seqnames=genome, start=seq(1, l_genome, window_step), end=seq(1, l_genome, window_step) + window_size - 1)
  window_df$window_id <- paste0('window_', 1:nrow(window_df))
  window_df           <- plyr::rbind.fill(data.frame(window_df, strand='+'), data.frame(window_df, strand='-'))

  ### Merge windows with coverages
  win.cov      <- genome_join(window_df, cov.df, by=by)
  win.cov      <- win.cov[win.cov$seqnames.x == win.cov$seqnames.y,]
  win.cov      <- win.cov %>% dplyr::select(!any_of('seqnames.y'))
  colnames(win.cov)[1] <- by[1]
  colnames(win.cov)[grep('\\.x', colnames(win.cov))] <- gsub('.x', '.window', colnames(win.cov)[grep('\\.x', colnames(win.cov))])
  colnames(win.cov)[grep('\\.y', colnames(win.cov))] <- gsub('.y', '', colnames(win.cov)[grep('\\.y', colnames(win.cov))])
  win.cov      <- win.cov[win.cov$strand == win.cov$strand.window, ]
  win.cov      <- win.cov %>% dplyr::select(!any_of('strand.window'))

  ### Calculate coverage per window
  win.cov.sum  <- win.cov %>%
    group_by(seqnames, window_id, start.window, end.window, strand, sample) %>%
    summarise(window_mean.cov=mean(count), window_median.cov=median(count), window_min.cov=min(count), window_max.cov=max(count), window_sum.cov=sum(count))

  return(win.cov.sum)

}

