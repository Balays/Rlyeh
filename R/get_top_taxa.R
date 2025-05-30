
get_top_taxa <- function(ps, top = 20, method = 'sum', overlap = 'union', threshold_prevalence = 0.75) {

  require(phyloseq)
  require(data.table)

  # Validate input
  if (!method %in% c('sum', 'individual')) stop("method must be 'sum' or 'individual'")
  if (!overlap %in% c('union', 'intersect', 'core', 'mean')) stop("overlap must be 'union', 'intersect', 'core', or 'mean'")

  if (method == 'sum') {
    message("Method: summing taxa abundances across all samples.")
    top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:top]

  } else if (method == 'individual') {

    otutab <- data.table(data.frame(otu_table(ps)), keep.rownames = 'taxon')
    otutab_melt <- melt.data.table(otutab, id.vars = 'taxon', variable.name = 'sample', value.name = 'abundance')
    otutab_melt[, taxon_rank := frank(-abundance, ties.method = 'dense'), by = sample]
    otutab_top <- otutab_melt[taxon_rank <= top]

    if (overlap == 'union') {

      message("Taking the UNION of the top-", top, " taxa from each sample.")
      top_taxa <- unique(otutab_top$taxon)
      message("Number of taxa selected: ", length(top_taxa))

    } else if (overlap == 'intersect') {

      message("Taking the INTERSECT of the top-", top, " taxa from each sample.")
      taxa_per_sample <- otutab_top[, .(list(taxon)), by = sample]$V1
      top_taxa <- Reduce(intersect, taxa_per_sample)
      message("Number of taxa selected: ", length(top_taxa))

    } else if (overlap == 'core') {

      if (!is.null(prevalence_threshold)) {

        prevalence_threshold <- ceiling(length(unique(otutab_top$sample)) * threshold_prevalence)
        message("Identifying CORE taxa present in at least ", threshold_prevalence * 100, "% of samples, [ ",
                prevalence_threshold, "(N) samples ], ranked by mean abundance.")
      } else {
        message("Identifying CORE taxa present in at least ", prevalence_threshold, "(N) samples, ranked by mean abundance.")
      }

      core_candidates <- otutab_top[, .(
        n_samples = uniqueN(sample),
        mean_abundance = mean(abundance)
      ), by = taxon][n_samples >= prevalence_threshold]

      core_candidates <- core_candidates[order(-mean_abundance)]
      top_taxa <- head(core_candidates$taxon, top)

      message("Number of core taxa selected: ", length(top_taxa))

    } else if (overlap == 'mean') {

      message("Taking top-", top, " taxa based on mean abundance across samples.")

      mean_taxa_abundance <- otutab_top[, .(mean_abundance = mean(abundance)), by = taxon][order(-mean_abundance)]
      top_taxa <- head(mean_taxa_abundance$taxon, top)

      message("Number of taxa selected: ", length(top_taxa))
    }
  }

  return(top_taxa)
}
