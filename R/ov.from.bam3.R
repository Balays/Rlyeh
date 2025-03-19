#' Import alignments from a .BAM file.
#' Optional filtering can be applied for the alignments, based on bam flags, mapping quality, and reference sequences.
#' LoRTIA ouptuts can be imported, along with LoRTIA tags for further filtering.
#' Gaps (Ns) in the alignments can be removed into distinct alignments for the same read using the CIGAR values.
#'
#' @export

ov.from.bam3 <- function (bamfile,
                          what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand"),
                          flag = scanBamFlag(),
                          is.lortia = F, lortia.tags = c("l3", "l5", "r3", "r5"),
                          rm.non.correct=F, rm.false.exons=F,
                          crop.na.cigar = T,
                          rm.gaps.in.aln = F,
                          add.primes=T
                          ) {
  require(GenomicAlignments)
  require(Rsamtools)
  require(data.table)
  #bamfile <- 'I:/data/SARS-CoV2/mapped_v8/.bam.fastq.pych.bam/Hpi10_A.merged_pychopped.bam'

  ## Is this LoRTIA output?
  if (is.lortia) {
    params <- ScanBamParam(what = what, tag = lortia.tags, flag = flag)
    
    bam <- as.data.frame(scanBam(bamfile, param = params))
    setDT(bam)
    
    # Define column names dynamically
    tag_cols <- paste0("tag.", lortia.tags)
    
    # Create a temporary cleaned version of the tag columns without modifying the original
    cleaned_tags <- bam[, lapply(.SD, function(x) gsub(".*,", "", x)), .SDcols = tag_cols]
    
    # Paste the cleaned values together into a new 'tags' column
    bam[, tags := do.call(paste, c(cleaned_tags, sep = ","))]
      
    
    ## Remove alignments designated as having no correct adapter by LoRTIA
    if (rm.non.correct) {
      bam <- bam[grepl("correct", tags), ]
    }
    ## Remove alignments designated as having false exons by LoRTIA
    if (rm.false.exons) {
      bam <- bam[!grepl("false", tags), ]
    }
    
  ## Non-LoRTIA input  
  } else {
    params <- ScanBamParam(what = what, flag = flag)
    bam <- as.data.frame(scanBam(bamfile, param = params))
    setDT(bam)
  }



  ## separate unmapped reads or reads with NA CIGAR
  bam.na <- bam[ is.na(rname) |  is.na(cigar),  ]
  bam    <- bam[!is.na(rname) & !is.na(cigar),  ]

  ## Separate alignments with gaps into sub_alignments i.e. exons ?
  if (rm.gaps.in.aln) {
    bam$group <- seq(1, nrow(bam))
    aln  <- as.data.frame(extractAlignmentRangesOnReference(bam$cigar, pos = bam$pos, f = NULL))
    stopifnot(luniq(aln$group) == nrow(bam))
  } else {
    bam$group <- seq(1, nrow(bam))
    aln <- as.data.frame(cigarRangesAlongReferenceSpace(bam$cigar, pos = bam$pos, f = NULL, reduce.ranges=T))
    stopifnot(nrow(aln) == nrow(bam))
  }
  aln <- data.table(dplyr::select(aln, -group_name))
  aln <- merge(aln, bam, by = "group")[, -1]
  bam <- aln

  ## put back unmapped reads
  bam <- rbind(bam, bam.na, use.names=T, fill=T)
  setDT(bam)
  ##

  ##
  if (is.lortia) {
    cnames <- c(what, "start", "end", paste0("tag.", lortia.tags), "tags")
    bam <- bam[, ..cnames]
  } else {
    cnames <- c(what, "start", "end")
    bam <- bam[, ..cnames]
  }
  colnames(bam)[1] <- "seqnames"
  

  ## add 3' and 5' end positions
  if (add.primes) {
    bam[strand == '+', prime5 := start ]
    bam[strand == '+', prime3 := end   ]
    bam[strand == '-', prime5 := end   ]
    bam[strand == '-', prime3 := start ]
  }


  return(bam)

}

