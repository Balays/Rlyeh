#' 
#' 
#' @export
#'
#'

require(tidyr)
require(fuzzyjoin)
require(data.table)


##
add.correct.ends <- function(bam.lortia, correct.strings.5p='correct', correct.strings.3p='correct') {
  bam.lortia$correct_tes  <- F
  bam.lortia$correct_tss  <- F
  bam.lortia$correct_tss[grepl(paste(correct.strings.5p, collapse='|'), bam.lortia$tag.l5) |
                           grepl(paste(correct.strings.5p, collapse='|'), bam.lortia$tag.r5)] <- T
  bam.lortia$correct_tes[grepl(paste(correct.strings.3p, collapse='|'), bam.lortia$tag.l3) |
                           grepl(paste(correct.strings.3p, collapse='|'), bam.lortia$tag.r3)] <- T
  return(bam.lortia)
}
