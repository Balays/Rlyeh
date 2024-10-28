#' Make plot data for the genome plot of plot.genome.regions function
#'
#' @export

gggenome.from.dataframe <- function (plot.data, visfrom =  NA, visto = NA, genename.toupper=T
                                     ) {

  require(gggenes)
  require(tidyverse)
  require(tidygenomics)
  require(fuzzyjoin)

  plot.data  <- plot.data[order(plot.data$strand, plot.data$start, plot.data$end), ]
  plot.data  <- data.frame(plot.data, xmin=plot.data$start, xmax = plot.data$end, ymin=0, ymax=0)

  plot.data$prime5[plot.data$strand == '-'] <- plot.data$end[plot.data$strand   == '-']
  plot.data$prime3[plot.data$strand == '-'] <- plot.data$start[plot.data$strand == '-']
  plot.data$prime5[plot.data$strand == '+'] <- plot.data$start[plot.data$strand == '+']
  plot.data$prime3[plot.data$strand == '+'] <- plot.data$end[plot.data$strand   == '+']

  plot.data  <-  genome_cluster(plot.data, by=c('seqnames', 'start', 'end'), 1, 'cluster')

  plot.data$ymin <- 0

  for (i in unique(plot.data$cluster)) {
    plot.data$ymin[plot.data$cluster == i] <- seq(1, nrow(plot.data[plot.data$cluster == i, ]))
  }


  plot.data$xmean <- apply(plot.data[,c('xmin', 'xmax')], 1, mean)
  plot.data$ymean <- apply(plot.data[,c('ymin', 'ymax')], 1, mean)

  if(!is.na(visfrom) & !is.na(visto)) {
    plot.data  <- plot.data[plot.data$start >= visfrom & plot.data$end <= visto, ]
  }

  plot.data$strand <- factor(plot.data$strand, levels = c('+', '-'))

  plot.data$orientation <- NA
  plot.data$orientation[plot.data$strand == '-'] <- 0
  plot.data$orientation[plot.data$strand == '+'] <- 1

  if(genename.toupper) {
    plot.data$gene <- toupper(plot.data$gene)
  }

  gene.plotdata <- plot.data
}
