

tax_glom_fast2 <- function (ps, rank = NULL, rank_level = NULL, ignore_lineage = T,
          ignore_lineage_below_rank = T)
{
  require(data.table)
  ranks <- rank_names(ps)
  if (is.null(rank) && is.null(rank_level)) {
    stop("Either 'rank' or 'rank_level' must be provided.")
  }
  if (!is.null(rank) && !is.null(rank_level)) {
    message("Both 'rank' and 'rank_level' are provided. Using 'rank' to determine 'rank_level'.")
  }
  if (!is.null(rank)) {
    rank_level <- which(ranks == rank)
    if (length(rank_level) == 0) {
      stop("The provided 'rank' does not exist in 'ranks'.")
    }
  }
  else {
    if (rank_level > length(ranks) || rank_level < 1) {
      stop("The provided 'rank_level' is out of range.")
    }
    rank <- ranks[rank_level]
  }
  samp_names <- sample_names(ps)
  otutab <- data.frame(otu_table(ps))
  otutab <- data.table(rownames = rownames(otutab), otutab)
  taxtab <- data.frame(tax_table(ps))
  taxtab <- data.table(rownames = rownames(taxtab), taxtab)
  sampdat <- data.frame(sample_data(ps))
  samples <- rownames(sampdat)
  colnames(otutab)[-1] <- samp_names
  psdata <- merge(taxtab, otutab, by = "rownames", all = T)
  stopifnot(all(!is.na(psdata$rownames)))
  if (ignore_lineage) {
    ranks <- rank
    message("Ignoring lineages when glomerating.")
    psdata[, `:=`(lineage, get(rank))]
    stopifnot(nrow(psdata[get(rank) != lineage]) == 0)
    psdata[is.na(lineage), `:=`(lineage, "NA-lineage")]
  }
  else {
    if (ignore_lineage_below_rank) {
      ranks <- ranks[1:rank_level]
      message("Considering lineages, but ignoring rank levels below \"",
              rank, "\" when glomerating.")
      psdata[, `:=`(lineage, apply(.SD, 1, function(row) paste(na.omit(row),
                                                               collapse = ";"))), .SDcols = ranks]
    }
    else {
      ranks <- ranks
      message("Considering lineages on all rank levels when glomerating.")
      psdata[, `:=`(lineage, apply(.SD, 1, function(row) paste(na.omit(row),
                                                               collapse = ";"))), .SDcols = ranks]
    }
    psdata[is.na(get(rank)), `:=`(lineage, paste("unknown",
                                                 lineage, rank, sep = "_"))]
  }
  cols_to_group_by <- c("lineage", ranks)
  psdata.m <- melt.data.table(psdata, id.vars = cols_to_group_by,
                              measure.vars = samples, variable.name = "sample", value.name = "count")
  cols_to_sum <- c(cols_to_group_by, "sample")
  psdata.m.glom <- psdata.m[, .(glom_count = sum(count)), by = cols_to_sum]
  formula_str <- paste(paste(cols_to_group_by, collapse = "+"),
                       " ~ sample")
  formula <- as.formula(formula_str)
  psdata.glom <- dcast.data.table(psdata.m.glom, formula, value.var = "glom_count")
  taxtab <- data.frame(psdata.glom[, ..ranks], row.names = psdata.glom$lineage)
  otutab <- data.frame(psdata.glom[, ..samples], row.names = psdata.glom$lineage)
  colnames(otutab)[] <- samp_names
  ps.glom <- NULL
  try({
    ps.glom <- phyloseq(tax_table(as.matrix(taxtab)), otu_table(as.matrix(otutab),
                                                                taxa_are_rows = T), sample_data(sampdat))
  })
  if (is.null(ps.glom)) {
    rownames(taxtab) <- NULL
    rownames(otutab) <- NULL
    ps.glom <- phyloseq(tax_table(as.matrix(taxtab)), otu_table(as.matrix(otutab),
                                                                taxa_are_rows = T), sample_data(sampdat))
    taxtab <- data.frame(tax_table(ps.glom))
    otutab <- data.frame(otu_table(ps.glom))
    taxa_names(ps.glom) <- taxtab[, rank]
  }
  return(ps.glom)
}
