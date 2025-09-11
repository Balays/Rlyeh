

cover_dt_weighted <- function(dt,
                              weight_col = NULL,
                              by = c("seqnames", "strand", "sample"),
                              return = c("rle", "pos"),
                              nworkers = 1L,
                              # zero-fill options
                              include_zero = FALSE,
                              contigs_dt = NULL,
                              contig_seq_col = "seqnames",
                              contig_start_col = "start",
                              contig_end_col   = "end") {
  stopifnot(all(c("seqnames","strand","start","end","sample") %in% names(dt)))
  return <- match.arg(return)
  setDT(dt)
  
  # Coerce to ints & numeric weight
  dt[, `:=`(start = as.integer(start), end = as.integer(end))]
  if (!is.null(weight_col)) {
    stopifnot(weight_col %in% names(dt))
    set(dt, j = weight_col, value = as.numeric(dt[[weight_col]]))
  }
  if (is.null(weight_col)) {
    dt[, `__w__` := 1.0]
    on.exit(dt[, `__w__` := NULL], add = TRUE)
    weight_col <- "__w__"
  }
  
  # Drop 0 / NA / non-finite weights early
  dt <- dt[is.finite(get(weight_col)) & get(weight_col) != 0]
  
  # Safe +1 with overflow guard
  add1 <- function(x) {
    y <- x + 1L
    y[is.na(x)] <- NA_integer_
    y[ y < x ] <- .Machine$integer.max
    y
  }
  
  make_runs <- function(x, bycols, wcol) {
    # Weighted events per group: +w at start, -w at end+1
    ev <- x[, {
      w <- get(wcol)
      .(pos   = c(start, add1(end)),
        delta = c(w,     -w))
    }, by = bycols]
    
    # Combine same-position events
    ev <- ev[, .(delta = sum(delta)), by = c(bycols, "pos")]
    
    # Order by group columns and position
    do.call(setorder, c(list(ev), as.list(c(bycols, "pos"))))
    
    # Weighted coverage via cumsum of deltas per group
    ev[, cov := cumsum(delta), by = bycols]
    ev[, next_pos := shift(pos, type = "lead"), by = bycols]
    
    runs <- ev[!is.na(next_pos) & cov > 0,
               .(start = pos, end = next_pos - 1L, count = cov),
               by = bycols]
    runs[]
  }
  
  # Compute positive-coverage runs
  if (nrow(dt) == 0L) {
    runs <- data.table(matrix(ncol = length(by) + 3L, nrow = 0L))
    setnames(runs, c(by, "start", "end", "count"))
  } else if (nworkers > 1L) {
    requireNamespace("future.apply", quietly = TRUE)
    requireNamespace("future", quietly = TRUE)
    split_dt <- split(dt, by = by, drop = TRUE, keep.by = TRUE)
    
    oplan <- NULL
    if (!inherits(future::plan(), "multiprocess")) {
      oplan <- future::plan(future::multisession, workers = nworkers)
      on.exit({ if (!is.null(oplan)) future::plan(oplan) }, add = TRUE)
    }
    
    parts <- future.apply::future_lapply(
      split_dt, make_runs, bycols = by, wcol = weight_col, future.seed = TRUE
    )
    runs <- rbindlist(parts, use.names = TRUE)
  } else {
    runs <- make_runs(dt, by, weight_col)
  }
  
  # === Zero-fill RLE based on contig spans ===
  if (isTRUE(include_zero)) {
    stopifnot(!is.null(contigs_dt),
              all(c(contig_seq_col, contig_start_col, contig_end_col) %in% names(contigs_dt)))
    
    g0 <- unique(as.data.table(contigs_dt)[,
                                           .(seqnames = get(contig_seq_col),
                                             g_start  = as.integer(get(contig_start_col)),
                                             g_end    = as.integer(get(contig_end_col)))
    ])
    
    # Expand contigs across any missing by-columns using dt's uniques
    comb_list <- list(seqnames = unique(g0$seqnames))
    for (col in setdiff(by, "seqnames")) {
      if (col %in% names(contigs_dt)) comb_list[[col]] <- unique(contigs_dt[[col]])
      else                            comb_list[[col]] <- unique(dt[[col]])
    }
    comb <- do.call(CJ, c(comb_list, list(unique = TRUE)))
    contigs_exp <- g0[comb, on = "seqnames", nomatch = 0L]
    
    # Groups with no positive runs -> full-length zero run
    runs_keys <- unique(runs[, ..by])
    missing_groups <- contigs_exp[!runs_keys, on = by]
    zero_missing <- if (nrow(missing_groups)) {
      missing_groups[, .(start = g_start, end = g_end, count = 0.0), by = by]
    } else data.table(matrix(ncol = length(by) + 3L, nrow = 0L))
    if (nrow(zero_missing)) setnames(zero_missing, c(by, "start", "end", "count"))
    
    # For groups with positives, add leading/trailing and between-run zeros
    if (nrow(runs)) {
      runs2 <- contigs_exp[runs, on = by]
      # g_start/g_end now available per group
      zeros_between <- runs2[, {
        sd <- copy(.SD)
        setorder(sd, start, end)
        gs <- g_start[1L]; ge <- g_end[1L]
        out <- vector("list", 3L)
        k <- 0L
        if (nrow(sd) && sd$start[1L] > gs) { k <- k + 1L; out[[k]] <- data.table(start = gs, end = sd$start[1L] - 1L, count = 0.0) }
        if (nrow(sd) > 1L) {
          gaps_s <- sd$end[-.N] + 1L
          gaps_e <- sd$start[-1L] - 1L
          sel <- gaps_s <= gaps_e
          if (any(sel)) { k <- k + 1L; out[[k]] <- data.table(start = gaps_s[sel], end = gaps_e[sel], count = 0.0) }
        }
        if (nrow(sd) && sd$end[.N] < ge) { k <- k + 1L; out[[k]] <- data.table(start = sd$end[.N] + 1L, end = ge, count = 0.0) }
        if (k) rbindlist(out[1:k]) else data.table(start = integer(), end = integer(), count = numeric())
      }, by = by]
      
      runs <- rbindlist(list(runs, zeros_between, zero_missing), use.names = TRUE, fill = TRUE)
    } else {
      runs <- zero_missing
    }
    
    if (nrow(runs)) do.call(setorder, c(list(runs), as.list(c(by, "start", "end"))))
  }
  
  if (return == "rle") {
    setcolorder(runs, c(by, setdiff(names(runs), by)))
    return(runs[])
  }
  
  # Per-position expansion (row-wise, length-safe). This will also include zeros if include_zero=TRUE.
  runs[, {
    pos_list <- Map(seq.int, start, end)
    .(pos   = unlist(pos_list, use.names = FALSE),
      count = rep.int(count, lengths(pos_list)))
  }, by = by]
}
