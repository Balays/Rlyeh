### what is the ratio of the mapped part of the whole read?
calc.alignment.ratio <- function(bam) {
  
  #aln.sum <- bam %>% group_by(sample, qname) %>% summarise(qwidth, sum.aln.width=sum(width), width.ratio=sum.aln.width/qwidth)
  #aln.sum <- unique.data.frame(aln.sum)
  
  ##
  setDT(bam)
  # Group by sample and qname, and then compute summaries
  aln.sum <- bam[, .(qwidth = mean(qwidth),  # assuming qwidth is the same for each group; using mean as a placeholder
                     sum.aln.width = sum(qwidth),
                     width.ratio = sum(qwidth) / mean(qwidth)), 
                 by = .(sample, qname)]
  
  # Remove duplicate rows
  aln.sum <- unique(aln.sum)
  ##
  
  bam <- merge(aln.sum, bam, by=c('sample', 'qname', 'qwidth'), all=T)
  return(bam)
}


### how many contigs each read has been mapped onto? (detection of chimeras)
calc.n_seq <- function(bam) {
  
  #bam.uni.seq <- unique.data.frame(bam %>% select(sample, qname, seqnames))
  
  ##
  setDT(bam)
  # Select desired columns and remove duplicates
  bam.uni.seq <- unique(bam[, .(sample, qname, seqnames)])
  
  # Group by sample and qname, and then compute summaries
  #aln.sum <- bam.uni.seq[, .(seqnames = list(seqnames), n_seq = .N, max_n = max(.N)),
  #                       by = .(sample, qname)]
  aln.sum <- bam.uni.seq %>% group_by(sample, qname) %>% summarise(seqnames, n_seq=seq_along(seqnames), max_n=max(n_seq))
  # (dplyr was much faster here)
  setDT(aln.sum)
  bam <- merge(aln.sum, bam, by=c('sample', 'qname', 'seqnames'), all=T)
  return(bam)
}


### how many organims were each read mapped onto?
calc.n_org <- function(bam) {
  #bam.uni.seq <- unique.data.frame(bam %>% select(sample, qname, seqnames, org))
  
  #aln.sum     <- bam.uni.seq %>% group_by(sample, qname, org) %>% summarise(n_seq_org=n())
  #bam         <- merge(aln.sum, bam, by=c('sample', 'qname', 'org'), all=T)
  
  #org.sum     <- aln.sum     %>% group_by(sample, qname, org) %>% summarise(n_org=n())
  #org.sum     <- org.sum     %>% group_by(sample, qname) %>% summarise(n_org=sum(n_org))
  
  ##
  setDT(bam)
  
  # Select desired columns and remove duplicates
  bam.uni.seq <- unique(bam[, .(sample, qname, seqnames, org)])
  
  # Group by sample, qname, and org, then compute the summary
  aln.sum <- bam.uni.seq[, .N, by = .(sample, qname, org)]
  setnames(aln.sum, "N", "n_seq_org")
  
  # Merge the summary with the original data
  bam <- merge(bam, aln.sum, by = c("sample", "qname", "org"), all = TRUE)
  
  # Compute further summaries based on the previous results
  org.sum <- aln.sum[, .N, by = .(sample, qname, org)]
  setnames(org.sum, "N", "n_org")
  org.sum <- org.sum[, .(n_org = sum(n_org)), by = .(sample, qname)]
  ##
  
  bam         <- merge(org.sum, bam, by=c('sample', 'qname'), all=T)
  
  return(bam)
}


### does the reads have supplementary or secondary alignments?

has_secondary <- function(flag) {
  
  result <- names(bam_flag_analyser(flag))
  
  if (any(grepl('secondary', result))) { result <- T } else { result <- F }
  
  return(result)
  
}

has_supplementary <- function(flag) {
  
  result <- names(bam_flag_analyser(flag))
  
  if (any(grepl('supplementary', result))) { result <- T } else { result <- F }
  
  return(result)
  
}






### LEGACY
map.eval <- function(bam, bam.flags) {
  
  ## 
  sec.flag   <- 256
  sec.flags  <- c(sec.flag, bam.flags$Decimal + sec.flag)
  
  qname.sec <- bam$qname[is.element(bam$flag, c(sec.flags)) ]
  bam$has.secondary <- F
  bam$has.secondary[is.element(bam$qname, qname.sec) ] <- T
  
  
  supp.flag   <- 2048
  supp.flags  <- c(supp.flag, bam.flags$Decimal + supp.flag)
  
  qname.supp  <- bam$qname[is.element(bam$flag, c(supp.flags)) ]
  bam$has.supplementary <- F
  bam$has.supplementary[is.element(bam$qname, qname.supp) ] <- T
  
  return(bam)
}


