## ------------------------------------------------------------
## Compute the genomic interval to visualise for one gene
##   • gene_id      – gene (or cluster) to plot
##   • end_type     – "TSS" or "TES"
##   • clusters_dt  – data.table with the columns used below
## Returns a named list: seqnames, strand, start, end
## ------------------------------------------------------------
get_vis_interval <- function(gene_id,
                             end_type   = c("TSS", "TES"),
                             clusters_dt) {
  
  end_type <- match.arg(end_type)
  
  ## -- anchor: seqname / strand / initial window -------------
  visregion <- unique(
    clusters_dt[gene == gene_id,
                .(seqnames, strand, visfrom, visto)]
  )
  if (nrow(visregion) != 1L)
    stop("Gene must map to exactly one cluster / strand combination")
  
  visfrom <- visregion$visfrom
  visto   <- visregion$visto
  strand  <- visregion$strand          # kept only for completeness
  
  ## -- expand toward alternative region ----------------------
  if (end_type == "TSS") {
    
    alt <- as.integer(
      clusters_dt[gene == gene_id,
                  .(TSS.region.start, TSS.region.end)]
    )
    
    visfrom <- min(visfrom, alt[1])    # **FIX** – use min, not max
    visto   <- max(visto,   alt[2])    # **FIX** – use max, not min
    
  } else {  # TES ------------------------------------------------
    
    alt <- as.integer(
      clusters_dt[gene == gene_id,
                  .(TES.region.start, TES.region.end)]
    )
    
    visfrom <- min(visfrom,
                   min(clusters_dt[gene_cluster == gene_id, visfrom]),
                   alt[1])             # **unchanged**
    visto   <- max(visto,
                   max(clusters_dt[gene_cluster == gene_id, visto]),
                   alt[2])             # **unchanged**
  }
  
  data.table(seqnames = visregion$seqnames,
       strand   = strand,
       start    = visfrom,
       end      = visto)
}
