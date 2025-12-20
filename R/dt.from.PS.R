


dt.from.ps <- function(ps,
                       glom.samples='Vregion',
                       keep.all=T,
                       glom_taxa=T, i=NA, t=NA, ignore_lineage = T) {
  
  
  if(!is.na(glom_taxa)) {
    if(is.na(i))   { i  <- 2      }
    if(is.na(t))   { t  <- rank_names(ps)[i]      }
    ps <- tax_glom_fast2(ps, rank = t, ignore_lineage = ignore_lineage)
    message('Taxa were combined at ', t, ' level.')
  }
  
  taxtable  <- as.data.frame(tax_table(ps) ) #[,c(1:n)])
  otutable  <- as.data.frame(otu_table(ps, taxa_are_rows = T))
  samptable <- data.table(data.frame(sample_data(ps))
  
  psdf <- merge(
    taxtable,
    otutable,
    by=0, all=keep.all
  )[,-1]
  
  if(!is.na(glom.samples)) {
    psdt <- data.table(melt(psdf, id.vars = 1, value.name = 'count', variable.name = 'sample'))
    #'sample',
    psdt <- merge(psdt, samptable, by = 'sample', all.x=T)
    
    psdt <- psdt[, .(count = sum(count)), by = .(get(glom.samples), get(t))]
    colnames(psdt)[1:2] <- c(glom.samples, t)
    
  } else {
    psdt <- data.table(psdf)
    
  }
  
  return(psdt)
  
}
