


###
#######

foverlaps3 <- function(DTx, DTy,
         by=c('seqnames', 'strand', 'start', 'end'),
         by.x=by, by.y=by,
         
         #type='any',
         #minoverlap=20,
         #maxgap=0,
         
         prime5.win.start =  10, ## distance of prime5s of reference.TR and read.TR
         prime5.win.end   = -10,
         
         prime3.win.start = -10,
         prime3.win.end   =  10,
         
         intron.junction.wobble = 0
         
) {
  
  setkeyv(DTx, by.x)
  setkeyv(DTy, by.y)
  
  DTxy <- foverlaps(DTx, DTy, minoverlap = 1, maxgap = 0, type = 'any')
  ## .i colimns are from DTx (read-transcripts)
  ## coumns without .i are from DTy (ref transcripts)
  
  ## calculate overlap sizes
  DTxy[, overlap_size := pmax(0, pmin(end, i.end) - pmax(start, i.start) + 1)]
  
  ## calculate distance of read-TR from reference TR
  DTxy[, dist_prime5 := prime5 - i.prime5]
  DTxy[, dist_prime3 := prime3 - i.prime3]
  
  ### filter for distances
  
  ## prime5
  if (!is.null(prime5.win.start) & !is.null(prime5.win.end)) {
    message('Filtering overlap results based on 5-prime distances...')
    DTxy.filt <- DTxy[dist_prime5 <= prime5.win.start &
                      dist_prime5 >= prime5.win.end   
                         ]
  } else if (!is.null(prime5.win.start) &  is.null(prime5.win.end)) { stop('either prime5.win.start AND prime5.win.end should be NULL or NEITHER !') 
  } else if ( is.null(prime5.win.start) & !is.null(prime5.win.end)) { stop('either prime5.win.start AND prime5.win.end should be NULL or NEITHER !')
  } else (message('Overlap results were not filtered based on 5-prime distances.'))
  
  ## prime3
  if (!is.null(prime3.win.start) & !is.null(prime3.win.end)) {
    message('Filtering overlap results based on 3-prime distances...')
    DTxy.filt <- DTxy[dist_prime3 >= prime3.win.start &
                      dist_prime3 <= prime3.win.end   
    ]
  } else if (!is.null(prime3.win.start) &  is.null(prime3.win.end)) { stop('either prime5.win.start AND prime5.win.end should be NULL or NEITHER !') 
  } else if ( is.null(prime3.win.start) & !is.null(prime3.win.end)) { stop('either prime5.win.start AND prime5.win.end should be NULL or NEITHER !') 
  } else (message('Overlap results were not filtered based on 3-prime distances.'))
  
  
  
  return(DTxy)
}





