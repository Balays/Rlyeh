#' Windowed coverage over RLE runs (no per-base expansion)
#'
#' @param runs data.table with columns: by..., start, end, count (from cover_dt_weighted(return='rle'))
#' @param contigs_dt optional data.table with seqnames,start,end (and optionally strand/sample) to define spans
#' @param by grouping columns in runs (default c('seqnames','strand','sample'))
#' @param window_size integer window width (bp)
#' @param step integer step between windows (bp)
#' @param stat one of 'sum','mean','median' (default 'sum')
#' @param tail whether to include a partial tail window at the contig end: 'drop' or 'keep' (default 'drop')
#' @param nworkers parallel groups via future.apply multisession (Windows-safe)
#' @return data.table with by..., w_start, w_end, value
window_coverage <- function(runs,
                            contigs_dt = NULL,
                            by = c("seqnames","strand","sample"),
                            window_size = 1000L,
                            step = 1000L,
                            stat = c("sum","mean","median"),
                            tail = c("drop","keep"),
                            nworkers = 1L,
                            contig_seq_col = "seqnames",
                            contig_start_col = "start",
                            contig_end_col   = "end") {
  stopifnot(all(c(by, "start","end","count") %in% names(runs)))
  stat <- match.arg(stat)
  tail <- match.arg(tail)
  setDT(runs)
  # ensure integer coordinates
  runs[, `:=`(start = as.integer(start), end = as.integer(end))]
  
  # Build per-group spans
  if (!is.null(contigs_dt)) {
    stopifnot(all(c(contig_seq_col, contig_start_col, contig_end_col) %in% names(contigs_dt)))
    g0 <- unique(as.data.table(contigs_dt)[,
                                           .(seqnames = get(contig_seq_col),
                                             g_start  = as.integer(get(contig_start_col)),
                                             g_end    = as.integer(get(contig_end_col)))
    ])
    # Expand across missing by columns using runs' uniques
    comb_list <- list(seqnames = unique(g0$seqnames))
    for (col in setdiff(by, "seqnames")) {
      if (col %in% names(contigs_dt)) comb_list[[col]] <- unique(contigs_dt[[col]])
      else                            comb_list[[col]] <- unique(runs[[col]])
    }
    comb <- do.call(CJ, c(comb_list, list(unique = TRUE)))
    spans <- g0[comb, on = "seqnames", nomatch = 0L]
  } else {
    spans <- runs[, .(g_start = min(start), g_end = max(end)), by = by]
  }
  setkeyv(spans, by)
  
  # helper: weighted median for a single window
  wmedian <- function(v, w) {
    if (!length(v)) return(0)
    o <- order(v)
    v <- v[o]; w <- w[o]
    cw <- cumsum(w)
    half <- sum(w) / 2
    v[which(cw >= half)[1L]]
  }
  
  # worker per group
  worker <- function(x) {
    # group key
    k <- x[1L, ..by]
    # contig span for this group
    sp <- spans[k, on = by]
    if (!nrow(sp)) return(NULL)
    gs <- sp$g_start[1L]; ge <- sp$g_end[1L]
    if (!is.finite(gs) || !is.finite(ge) || ge < gs) return(NULL)
    
    # windows
    if (tail == "drop") {
      s_max <- ge - window_size + 1L
      if (s_max < gs) return(NULL)
      w_st <- seq.int(gs, s_max, by = step)
      w_en <- w_st + window_size - 1L
    } else { # keep tail
      w_st <- seq.int(gs, ge, by = step)
      w_en <- pmin(w_st + window_size - 1L, ge)
    }
    if (!length(w_st)) return(NULL)
    win <- data.table(w_id = seq_along(w_st), w_start = w_st, w_end = w_en)
    win[, width := (w_end - w_start + 1L)]
    
    # overlaps
    xx <- x[, .(start, end, count)]
    setkey(xx, start, end)
    setkey(win, w_start, w_end)
    ov <- foverlaps(xx, win, by.x = c("start","end"), by.y = c("w_start","w_end"), type = "any", nomatch = 0L)
    if (!nrow(ov)) {
      # no overlaps -> zeros
      out <- cbind(k[rep(1L, nrow(win))], win[, .(w_start, w_end)])
      out[, value := 0]
      return(out)
    }
    ol <- pmin(ov$end, ov$w_end) - pmax(ov$start, ov$w_start) + 1L
    ov <- ov[ol > 0]
    if (!nrow(ov)) {
      out <- cbind(k[rep(1L, nrow(win))], win[, .(w_start, w_end)])
      out[, value := 0]
      return(out)
    }
    ov[, ol := pmin(end, w_end) - pmax(start, w_start) + 1L]
    
    if (stat == "sum" || stat == "mean") {
      agg <- ov[, .(sumcov = sum(count * ol)), by = w_id]
      out <- merge(win[, .(w_id, w_start, w_end, width)], agg, by = "w_id", all.x = TRUE)
      out[is.na(sumcov), sumcov := 0]
      out[, value := if (stat == "sum") sumcov else sumcov / width]
      out <- cbind(k[rep(1L, nrow(out))], out[, .(w_start, w_end, value)])
      return(out[])
    } else { # median
      med <- ov[, .(value = wmedian(count, ol)), by = w_id]
      out <- merge(win[, .(w_id, w_start, w_end)], med, by = "w_id", all.x = TRUE)
      out[is.na(value), value := 0]
      out <- cbind(k[rep(1L, nrow(out))], out[, .(w_start, w_end, value)])
      return(out[])
    }
  }
  
  # split by groups
  parts <- split(runs, by = by, drop = TRUE, keep.by = TRUE)
  if (nworkers > 1L) {
    requireNamespace("future.apply", quietly = TRUE)
    requireNamespace("future", quietly = TRUE)
    oplan <- NULL
    if (!inherits(future::plan(), "multiprocess")) {
      oplan <- future::plan(future::multisession, workers = nworkers)
      on.exit({ if (!is.null(oplan)) future::plan(oplan) }, add = TRUE)
    }
    out_list <- future.apply::future_lapply(parts, worker, future.seed = TRUE)
  } else {
    out_list <- lapply(parts, worker)
  }
  ans <- rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(ans)) return(ans)
  setcolorder(ans, c(by, setdiff(names(ans), by)))
  ans[]
}