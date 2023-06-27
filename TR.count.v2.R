#' This function finds reference transcripts in the query transcript list
#' 
#'
#' @export
#'
#'

require(tidyr)
require(fuzzyjoin)

TR.count <- function(ex.ref, TR.EX,
                     #tr.counts=NULL,
                     #bam.filt=NA,
                     maxgap=35, minoverlap=0, type="equal",
                     intron.gap=1,
                     by = c('seqnames', 'start', 'end')) {

  #bam.test <- read.delim('bam.test.tsv'); TR.EX  <- bam.test

  ### Dropping non-mandatory columns
  TR.EX  <- TR.EX[,c('TR_ID', 'EX_ID', "strand", by)]

  ex.ref <- ex.ref[,c('TR_ID', 'EX_ID', "strand", by)]

  ### order exons
  order.exons  <- function(TR.EX, nr=T) {
    TR.plus   <- TR.EX[TR.EX$strand == '+', ]
    TR.minus  <- TR.EX[TR.EX$strand == '-', ]
    TR.plus   <- TR.plus[order(TR.plus$start), ]
    TR.minus  <- TR.minus[order(TR.minus$end, decreasing = T), ]
    TR.EX     <- data.frame( plyr::rbind.fill(TR.plus, TR.minus)  )
    if (nr) {
      TR.EX     <- TR.EX %>% group_by(TR_ID) %>% 
        summarise(EX_ID, seqnames, start, end, strand, EX_nr=seq_along(EX_ID), EX_count=max(EX_nr))
    }

  }

  TR.EX  <- order.exons(TR.EX)
  ex.ref <- order.exons(ex.ref)


  #### find overlaps of features and Transcripts
  tr.ov <- genome_join(ex.ref[,],
                       TR.EX[,],
                       type=type, by=by, maxgap=maxgap)
  colnames(tr.ov)[grep('.x', colnames(tr.ov))] <- gsub('.x', '.ref',   colnames(tr.ov)[grep('.x', colnames(tr.ov))])
  colnames(tr.ov)[grep('.y', colnames(tr.ov))] <- gsub('.y', '.query', colnames(tr.ov)[grep('.y', colnames(tr.ov))])

  ### Overlapping Seqnames
  tr.ov <- tr.ov[tr.ov$seqnames.ref == tr.ov$seqnames.query, ]
  tr.ov <- tr.ov  %>% dplyr::select(!('seqnames.ref'))
  colnames(tr.ov)[grep('seqnames', colnames(tr.ov))] <- 'seqnames'

  ### Overlapping Strands
  tr.ov   <- tr.ov[tr.ov$strand.ref == tr.ov$strand.query, ]

  ### Other exons of query TRs
  TR.found <- unique(tr.ov$TR_ID.query)
  TR.found <- TR.EX[is.element(TR.EX$TR_ID, TR.found),]
  tr.ov    <- merge(tr.ov, TR.found,
                    by.x=c("EX_ID.query", "TR_ID.query",  "EX_nr.query", "EX_count.query", "seqnames", "start.query", "end.query", "strand.query"),
                    by.y=c("EX_ID","TR_ID", "EX_nr", "EX_count", by, 'strand'),
                    all=T)

  ### Other exons of reference TRs
  ## List of refernce transcripts
  tr.refs  <- na.omit(unique(tr.ov$TR_ID.ref))
  ##
  TR.ref   <- ex.ref[is.element(ex.ref$TR_ID, tr.refs),]
  tr.ov    <- merge(tr.ov, TR.ref,
                    by.x=c("EX_ID.ref", "TR_ID.ref",  "EX_nr.ref", "EX_count.ref", "seqnames", "start.ref", "end.ref", "strand.ref"),
                    by.y=c("EX_ID","TR_ID", "EX_nr", "EX_count", by, 'strand'),
                    all=T)

  ## this contains all overlaps, including splice variants

  ### from these, we keep only those TRs that have all exons in the ref and in the query as well

  tr.ov.sum <- tr.ov %>% group_by(TR_ID.ref, TR_ID.query) %>% summarise(EX_sum.ref=sum(EX_nr.ref), EX_sum.query=sum(EX_nr.query))

  tr.ov <- merge(tr.ov, tr.ov.sum, by=c('TR_ID.ref', 'TR_ID.query'))

  tr.ov$splice.var <- T

  i <- 'TR_505' #

  for (i in tr.refs) {

    ex.ref.sum <- sum(ex.ref$EX_nr[ex.ref$TR_ID == i])

    tr.ov$splice.var[tr.ov$TR_ID.ref == i &
                       tr.ov$EX_count.ref == tr.ov$EX_count.query &
                       tr.ov$EX_sum.ref   == ex.ref.sum &
                       tr.ov$EX_sum.query == ex.ref.sum  ] <- F

  }


  tr.ov.splice.var <- tr.ov[tr.ov$splice.var == T, ]

  tr.ov   <- tr.ov[tr.ov$splice.var == F, ]

  ### calculate distance betwwen ref and query hit
  nulldf <- data.frame(NULL)
  tr.ov$dist <- NA
  for (i in tr.refs) {
    tr.ref <- tr.ov[tr.ov$TR_ID.ref == i, ]

    if (all(tr.ref$strand.ref == '+')) {
      tr.ref$dist <- tr.ref$start.ref - tr.ref$start.query
    } else if (all(tr.ref$strand.ref == '-')) {
      tr.ref$dist <- tr.ref$end.query - tr.ref$end.ref
    } else {stop('wrong strand')}

    nulldf <- plyr::rbind.fill(nulldf, tr.ref)
  }
  tr.ov   <- nulldf

  ### Get the best hit for each query transcript
  tr.query.sum <- unique.data.frame(tr.ov[,c('TR_ID.query', 'TR_ID.ref')])
  ##
  tr.dups      <- dup(tr.query.sum$TR_ID.query)

  nulldf <- data.frame(NULL)
  i <- 'TR_108435'
  for (i in tr.dups) {
    tr.dup <- tr.ov [tr.ov$TR_ID.query ==  i, ]
    tr.dup <- tr.dup[tr.dup$EX_nr.ref  ==  1 & tr.dup$EX_nr.query == 1 & abs(tr.dup$dist) == min(abs(tr.dup$dist)), ]
    nulldf <- plyr::rbind.fill(nulldf,
                               unique.data.frame(tr.dup[,c('TR_ID.query', 'TR_ID.ref')]) )
  }

  tr.query.sum <- plyr::rbind.fill(tr.query.sum[!is.element(tr.query.sum$TR_ID.query, tr.dups),], nulldf)

  ## If there are hits with equal distance, consider the query as a longer variant of the shorter ref
  tr.dups      <- dup(tr.query.sum$TR_ID.query)

  nulldf <- data.frame(NULL)
  i <- 'TR_119432'
  for (i in tr.dups) {
    tr.dup <- tr.ov [tr.ov$TR_ID.query ==  i, ]
    tr.dup <- tr.dup[tr.dup$EX_nr.ref  ==  1 & tr.dup$EX_nr.query == 1 & tr.dup$dist == max(tr.dup$dist), ]
    nulldf <- plyr::rbind.fill(nulldf,
                               unique.data.frame(tr.dup[,c('TR_ID.query', 'TR_ID.ref')]) )
  }

  tr.query.sum <- plyr::rbind.fill(tr.query.sum[!is.element(tr.query.sum$TR_ID.query, tr.dups),], nulldf)

  ## remove dupliacted hits
  tr.ov.all.hits <- tr.ov
  tr.ov <- merge(tr.ov, tr.query.sum, by=c('TR_ID.query', 'TR_ID.ref'))

  tr.query.sum <- unique.data.frame(tr.ov[,c('TR_ID.query', 'TR_ID.ref')])
  tr.dups      <- dup(tr.query.sum$TR_ID.query)
  if(length(tr.dups) != 0) {
    message('The following query TRs are duplicated, as they can be associated with more than one reference TR: ', paste(tr.dups, collapse = ', '))
  }

  ### Introns on exact-matching
  
  ## one query transcript with multiple ref introns? ##
  
  ## introns in ref
  mrna.ex.ref <- ex.ref %>% group_by(TR_ID, seqnames, strand) %>% summarise(start=min(start), end=max(end))
  #nulldf <- data.frame(NULL)
  i <- 'asNOIRanalouge-SP2'
  #for (i in mrna.ex.ref$TR_ID) {
  introns.ref <- genome_subtract(mrna.ex.ref, #[mrna.ex.ref$TR_ID == i, ],
                                 ex.ref, #[ex.ref$TR_ID == i, ],
                                 by=c('TR_ID', 'start', 'end')) # by) #
  #  nulldf <- plyr::rbind.fill(nulldf, introns.ref)
  #}
  #introns.ref <- nulldf
  introns.ref.sum <- unique.data.frame(introns.ref[,c(by, 'strand')])
  introns.ref.sum <- data.frame(introns.ref.sum, intron_ID=paste0('intron_', 1:nrow(introns.ref.sum)))
  introns.ref     <- merge(introns.ref, introns.ref.sum, by=c(by, 'strand'))

  ## introns in query
  mrna.ex.query <- TR.EX %>% group_by(TR_ID, seqnames, strand) %>% summarise(start=min(start), end=max(end))
  #nulldf <- data.frame(NULL)
  #tr.found <- unique(tr.ov$TR_ID.query)
  #l_tr <- length(tr.found)
  #tr <- 'TR_105440'
  #for (i in 1:l_tr) { ### mrna.ex.query$TR_ID
  #tr <- tr.found[i]
  #message(i, ' of ', l_tr)
  introns.query <- genome_subtract(mrna.ex.query, # [mrna.ex.query$TR_ID == tr, ],
                                   TR.EX, # [TR.EX$TR_ID == tr, ],
                                   by=c('TR_ID', 'start', 'end')) # by) #
  #nulldf <- plyr::rbind.fill(nulldf, introns.query)
  #}
  #introns.query <- nulldf
  introns.query.sum <- unique.data.frame(introns.query[,c(by, 'strand')])
  introns.query.sum <- data.frame(introns.query.sum, intron_ID=paste0('intron_', 1:nrow(introns.query.sum)))
  introns.query     <- merge(introns.query, introns.query.sum, by=c(by, 'strand'))


  ## Overlapping introns
  introns.found <- genome_join(introns.ref.sum, introns.query.sum, by=by, maxgap=intron.gap, minoverlap=0, type="equal")

  colnames(introns.found)[grep('.x', colnames(introns.found))] <- gsub('.x', '.ref',   colnames(introns.found)[grep('.x', colnames(introns.found))])
  colnames(introns.found)[grep('.y', colnames(introns.found))] <- gsub('.y', '.query', colnames(introns.found)[grep('.y', colnames(introns.found))])

  # Overlapping Seqnames
  introns.found <- introns.found[introns.found$seqnames.ref == introns.found$seqnames.query, ]
  introns.found <- introns.found  %>% dplyr::select(!('seqnames.ref'))
  colnames(introns.found)[grep('seqnames', colnames(introns.found))] <- 'seqnames'

  # Overlapping Strands
  introns.found   <- introns.found[introns.found$strand.ref == introns.found$strand.query, ]

  ## Query transcripts with ref introns
  tr.ov.intron.var  <- tr.ov
  tr.found <- introns.query$TR_ID[is.element(introns.query$intron_ID, introns.found$intron_ID.query)]
  ## Non-spliced transcripts
  tr.found <- c(tr.found, TR.EX$TR_ID[TR.EX$EX_count == 1])

  tr.ov    <- tr.ov[is.element(tr.ov$TR_ID.query, tr.found),]

  ## OK

  ##
  tr.ov <- merge(tr.ov, tr.gt[,c("TR_ID", 'sample', "read_count")], by.x='TR_ID.query', by.y='TR_ID')

  ##
  tr.found.counts <- unique.data.frame(tr.ov[,c('TR_ID.query', 'TR_ID.ref', 'sample', 'read_count')])

  ## sum counts
  tr.found.sum.counts <- tr.found.counts %>% group_by(TR_ID.ref, sample) %>% summarise(sum.count=sum(read_count)) %>% spread(sample, sum.count, fill=0)

  return(list(tr.ov=tr.ov, tr.found.counts=tr.found.counts, tr.found.sum.counts=tr.found.sum.counts ))

}
