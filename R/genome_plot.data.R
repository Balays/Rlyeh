require(gggenes)
require(tidyverse)
require(tidygenomics)
require(fuzzyjoin)

make.genome.plot.data <- function(feature.df, feature.col=5, feature.name='gene', visfrom=NA, visto=NA) {

  genome.plotdata <- feature.df
  colnames(genome.plotdata)[feature.col] <- feature.name

  if (all(!is.na(c(visfrom, visto)))) {
    genome.plotdata  <- genome.plotdata[genome.plotdata$start >= visfrom & genome.plotdata$end <= visto, ]
  }

  genome.plotdata  <- genome.plotdata[order(genome.plotdata$strand, genome.plotdata$start, genome.plotdata$end), ]
  
  #genome.plotdata <- genome.plotdata %>% group_by(seqnames, strand, gene) %>% reframe(seqnames, strand, gene, ID, start, end, gene.start=min(start), gene.end=max(end))
  
  genome.plotdata  <- data.frame(genome.plotdata, xmin=genome.plotdata$start, xmax = genome.plotdata$end, ymin=0, ymax=0)

  genome.plotdata$prime5 <- genome.plotdata$start; genome.plotdata$prime3 <- genome.plotdata$end 
  genome.plotdata  <- add.primes(genome.plotdata)
   
  genome.plotdata$strand <- factor(genome.plotdata$strand, levels = c('+', '-', '*'))

  genome.plotdata  <- genome_cluster(genome.plotdata, by=c('seqnames', 'start', 'end'), 1, 'cluster')

  #genome.plotdata$ymin <- 0

  for (i in unique(genome.plotdata$cluster)) {
    genome.plotdata$ymin[genome.plotdata$cluster == i] <- seq(1, nrow(genome.plotdata[genome.plotdata$cluster == i, ]))
  }

  genome.plotdata$xmean <- apply(genome.plotdata[,c('xmin', 'xmax')], 1, mean)
  genome.plotdata$ymean <- apply(genome.plotdata[,c('ymin', 'ymax')], 1, mean)

  genome.plotdata$orientation <- NA
  genome.plotdata$orientation[genome.plotdata$strand == '-'] <- 0
  genome.plotdata$orientation[genome.plotdata$strand == '+'] <- 1

  genome.plotdata
}


