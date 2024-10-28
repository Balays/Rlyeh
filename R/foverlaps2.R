

foverlaps2 <- function(DTx, DTy,
                       by=c('seqnames', 'strand', 'start', 'end'),
                       by.x=by, by.y=by,
                       type='any',
                       minoverlap=20,
                       maxgap=0
                        ) {
  setkeyv(DTx, by.x)
  setkeyv(DTy, by.y)

  DTxy <- foverlaps(DTx, DTy, minoverlap = 1, maxgap = maxgap, type = type)  
  
  DTxy[, width_x      := end - start + 1]  [, width_y := i.end - i.start + 1]
  #DTxy[, overlap_size := end - i.start ]
  
  DTxy[, overlap_size := pmax(0, pmin(end, i.end) - pmax(start, i.start) + 1)]
  
  DTxy <- DTxy[overlap_size >= minoverlap]
  
  ## those transcript that are shorter than the minoverlap will be omitted !
  
  return(DTxy)
}




