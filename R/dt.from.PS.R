


dt.from.ps <- function(ps,
                       glom_taxa=F, i=NA, t=NA, ignore_lineage = T, keep.all=T,
                       melt = T,
                       glom_samples=F, summarise.samples = NULL, by.var = 'sample'
                       ) {


  if(glom_taxa) {
    if(is.na(i))   { i  <- 2      }
    if(is.na(t))   { t  <- rank_names(ps)[i]      }
    ps <- tax_glom_fast2(ps, rank = t, ignore_lineage = ignore_lineage)
    message('Taxa were combined at ', t, ' level.')
  }

  taxtable  <- data.frame(tax_table(ps) ) #[,c(1:n)])
  otutable  <- data.frame(otu_table(ps, taxa_are_rows = T))
  samptable <- data.frame(sample_data(ps))

  ranks     <- rank_names(ps)

  psdf <- merge(
    taxtable,
    otutable,
    by=0, all=keep.all
  )[,-1]

  setDT(psdf)

  if (melt) {

    cols <- intersect(colnames(psdf), c('lineage', ranks))

    psdt <- data.table(melt(psdf, id.vars = cols, value.name = 'count', variable.name = 'sample'))
    #'sample',

    psdt <- merge(psdt, samptable, by.x = 'sample', by.y = by.var, all.x=T)

    if(glom_samples) {

      psdt <- psdt[, .(sum_count  = sum( count),
                       mean_count = mean(count),
                       sd_count   = sd(  count)
      ), by = c(summarise.samples, cols)]

    } else if(!is.null(summarise.samples)) {


      psdt[, ':=' (sum_count  = sum( count),
                   mean_count = mean(count),
                   sd_count   = sd(  count)
      ), by = c(summarise.samples, cols)]

    }


  } else {

    psdt <- data.table(psdf)

  }

  return(psdt)

}

