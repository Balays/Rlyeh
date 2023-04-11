#' Visualize 5 or 3-primes on a genome
#'
#' @export

plot.genome.region <- function(visfrom, visto, plot.data, gene.plotdata, geom = geom_bar(), geom2 = NULL,
                               add.genome.plot=T, genome.only=F, add.cageTSS, cagefr.clust=NULL, cagefr.clust.dist=6,
                               gene.label=T, add.unstranded=T, add.feature=F, gene.feature.width.thresh=700, flip.gene.y=T, force.gene.y=T, force.all.gene.down=F, angle=-60,
                               prime='prime5', crop.FALSE=F, samples=NA, comb.all.samples=NA,
                               sizes=16, gene.sizes, gene.label.col='white', palette=pal_npg()(10), alpha=1, plot.title=NA,
                               tolower=F, vline=NULL, vline2=NULL, breakseq=5000, margins=unit(c(0,0,0,0), "cm"), #margins=margin(10,1,1,5)
                               scales='fixed', genomplot.scale=5.5, ylim=NULL, y.log10=F, y.multip=1.5,
                               strip.text =  element_text(angle = 45, size = sizes, hjust = 0.5),
                               gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene),
                               facet_cropF = facet_grid(rows = vars(hpi), cols=NULL, scales ),
                               ...) {


  #### Genome and annotation plot ####
  if(add.genome.plot) {
  visregion      <- data.frame(seqnames=genome, start=visfrom, end=visto)
  gene.plotsub   <- dplyr::select(genome_join(gene.plotdata, visregion, by=c('seqnames', 'start', 'end')), -c(seqnames.y, start.y, end.y) )
  colnames(gene.plotsub)[is.element(colnames(gene.plotsub), c('seqnames.x', 'start.x', 'end.x'))] <- c('seqnames', 'start', 'end')

  gene.plotsub$strand <- factor(gene.plotsub$strand, levels=c('+', '-', '*'))

  for (i in unique(gene.plotsub$cluster)) {
    gene.plotsub$ymin[gene.plotsub$cluster == i] <- seq(1, nrow(gene.plotsub[gene.plotsub$cluster == i, ]))
  }

  ### Add sub-gene positions
  gene.plotsub$subgene.start <- gene.plotsub$start
  gene.plotsub$subgene.end   <- gene.plotsub$end
  for (i in unique(gene.plotsub$gene)) {
    gene.plotsub$start[gene.plotsub$gene == i] <- min(gene.plotsub$start[gene.plotsub$gene == i])
    gene.plotsub$end  [gene.plotsub$gene == i] <- max(gene.plotsub$end  [gene.plotsub$gene == i])
  }
  gene.plotsub$width <- abs(gene.plotsub$end - gene.plotsub$start)+1

  ## make gene names lowercase
  if(tolower) {gene.plotsub$gene <- tolower(gene.plotsub$gene)}

  ## ORi regions strand !!!
  gene.plotsub$strand[grep('Ori', gene.plotsub$gene)] <- '*' #### !!!!

  ## flip
  if (flip.gene.y)  { gene.plotsub$ymin <- gene.plotsub$ymin * -1 }

  ## force into one row
  gene.plotsub$ylab <- gene.plotsub$ymin
  if (force.gene.y) { gene.plotsub$ymin <- 1 }

  ## stranded and unstranded features for geom_gene arrow and geom_rect
  gene.unstranded <- gene.plotsub[gene.plotsub$strand == '*', ]#; gene.unstranded$strand <- factor(gene.unstranded$strand, levels=  '*')
  gene.stranded   <- gene.plotsub[gene.plotsub$strand != '*', ]#; gene.stranded$strand   <- factor(gene.stranded$strand,   levels=c('+', '-'))

  ## gene labels pointing up and down
  if (force.gene.y) {
    gene.feature.upwards   <- as.data.frame(unique.data.frame(gene.stranded[gene.stranded$ylab  != Mode(gene.plotsub$ylab) |
                                                                            gene.stranded$width < gene.feature.width.thresh, c('seqnames', "gene")])) [c(T,F), 'gene']
    gene.feature.downwards <- as.data.frame(unique.data.frame(gene.stranded[gene.stranded$ylab  != Mode(gene.plotsub$ylab) |
                                                                            gene.stranded$width < gene.feature.width.thresh, c('seqnames', "gene")])) [c(F,T), 'gene']
  } else {
    gene.feature.upwards   <- as.data.frame(unique.data.frame(gene.stranded[gene.stranded$width < gene.feature.width.thresh, c('seqnames', "gene")])) [c(T,F), 'gene']
    gene.feature.downwards <- as.data.frame(unique.data.frame(gene.stranded[gene.stranded$width < gene.feature.width.thresh, c('seqnames', "gene")])) [c(F,T), 'gene']
  }

  if (force.all.gene.down) {
    gene.feature.upwards   <- NULL
    gene.feature.downwards <- gene.stranded$gene
  }

  ## add unstranded !!!
  gene.unstranded$start[gene.unstranded$width < 100] <- gene.unstranded$start - 20
  gene.unstranded$end  [gene.unstranded$width < 100] <- gene.unstranded$end   + 20

  ## force unstranded y position
  gene.unstranded$ymin <- max(gene.plotsub$ymin)

  ## recalculate xlim of plots to include gene annotations coverage
  visfrom    <- min(c(visfrom, min(gene.plotsub$start)))
  visto      <- max(c(visto,   max(gene.plotsub$end  )))
  visregion  <- data.frame(seqnames=genome, start=visfrom, end=visto)

  ### CAGE TSS clusters as genome feature
  if(add.cageTSS ) {
    cagefr.clust <- dplyr::select(genome_join(cagefr.clust, visregion, by=c('seqnames', 'start', 'end')), -c(seqnames.y, start.y, end.y) )
    colnames(cagefr.clust)[is.element(colnames(cagefr.clust), c('seqnames.x', 'start.x', 'end.x'))] <- c('seqnames', 'start', 'end')
    cagefr.clust$strand <- factor(cagefr.clust$strand, levels=c('+', '-', '*'))
    cagefr.clust$ymin  <-  max(gene.plotsub$ymin) + cagefr.clust.dist
    cagefr.clust$ymax  <-  max(gene.plotsub$ymin) + cagefr.clust.dist + 1
    cagefr.clust$y[cagefr.clust$strand == '-'] <- cagefr.clust$ymin
    cagefr.clust$y[cagefr.clust$strand == '+'] <- cagefr.clust$ymax
    cagefr.clust$region_start <- min(cagefr.clust$start)
    cagefr.clust$region_end   <- max(cagefr.clust$end)
    ## fix
    #cagefr.clust$strand[cagefr.clust$strand == '+'] <- '*'

  } else {cagefr.clust <- data.frame(ymin=max(gene.plotsub$ymin), ymax=max(gene.plotsub$ymin))}

  ### Set y-axis limits
  ylims <- c(min(c(gene.plotsub$ymin, cagefr.clust$ymin))-y.multip,
             max(c(gene.plotsub$ymin, cagefr.clust$ymax))+y.multip)


  ### Build genome annotation subplot
  gggenome <- ggplot(#gene.plotsub, #[gene.plotsub$strand != '*', ],
    mapping = gene.aes) +

    ## Gene arrows
    geom_gene_arrow(
        data = gene.stranded,
        arrowhead_height  = grid::unit(gene.sizes$gene_arrowhead_height,  "mm"),
        arrow_body_height = grid::unit(gene.sizes$gene_arrow_body_height, "mm"),
        alpha=alpha$gene_geom,
        fill='white') +

    ## Sub-gene arrows
    geom_subgene_arrow(
        aes(xsubmin = subgene.start, xsubmax = subgene.end, ),
        data = gene.stranded,
        arrowhead_height  = grid::unit(gene.sizes$gene_arrowhead_height,  "mm"),
        arrow_body_height = grid::unit(gene.sizes$gene_arrow_body_height, "mm"),
        alpha=alpha$gene_geom) +

    ## Gene arrow labels
    { if(gene.label) geom_gene_label(
        data  = gene.stranded[!is.element(gene.stranded$gene, c(gene.feature.upwards, gene.feature.downwards)),],
        height = grid::unit(gene.sizes$gene_label_height, "mm"), grow = F, colour=gene.label.col, align = "centre")
      } +

    ## Gene feature labels --> pointing upwards
    { if(gene.label & add.feature ) geom_feature(
      data = gene.stranded[is.element(gene.stranded$gene, gene.feature.upwards),  ],
      aes(x = xmean, y = ymin), forward = 1,
      feature_height = grid::unit(gene.sizes$gene_feature_height, "mm"),
      feature_width  = grid::unit(gene.sizes$gene_feature_width,  "mm"))
    } +
    { if(gene.label & add.feature ) geom_feature_label(
      data = gene.stranded[is.element(gene.stranded$gene, gene.feature.upwards), ],
      aes(x = xmean, y = ymin, label=gene), forward = 1,
      feature_height = grid::unit(gene.sizes$gene_feature_label_height, "mm"),
      label_height   = grid::unit(gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
    )
    } +
    ## Gene feature labels  --> pointing downwards
    { if(gene.label & add.feature ) geom_feature(
      data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      aes(x = xmean, y = ymin), forward = 1,
      feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
      feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm"))
    } +
    { if(gene.label & add.feature ) geom_feature_label(
      data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      aes(x = xmean, y = ymin, label=gene), forward = 1,
      feature_height = grid::unit(-gene.sizes$gene_feature_label_height, "mm"),
      label_height   = grid::unit(gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
    )
    } +
    ## Force all gene features pointing downwards
    { if(gene.label & force.all.gene.down ) geom_feature(
      data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      aes(x = xmean, y = ymin), forward = 1,
      feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
      feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm"))
    } +
    { if(gene.label & force.all.gene.down ) geom_feature_label(
      data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      aes(x = xmean, y = ymin, label=gene, angle=angle), forward = 1,
      feature_height = grid::unit(-gene.sizes$gene_feature_label_height, "mm"),
      label_height   = grid::unit(gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
    )
    } +

    ## Unstranded features with feature label
    { if(add.unstranded ) geom_rect(
        data=gene.unstranded,
        aes(xmin=start, xmax=end,
            ymin=ymin-gene.sizes$unstranded_rect_height/2,
            ymax=ymin+gene.sizes$unstranded_rect_height/2,
            fill=strand), color='black', fill='white')
      } +
    { if(add.unstranded ) geom_rect(
        data=gene.unstranded,
        aes(xmin=subgene.start, xmax=subgene.end,
            ymin=ymin-gene.sizes$unstranded_rect_height/2,
            ymax=ymin+gene.sizes$unstranded_rect_height/2,
            fill=strand), color='black', alpha=alpha$unstranded_geom)
    } +
    { if(add.unstranded & add.feature ) geom_feature(
        data = gene.unstranded,
        aes(x = xmean, y = ymin), forward = 1,
        feature_height = grid::unit(gene.sizes$unstranded_feature_height, "mm"),
        feature_width  = grid::unit(gene.sizes$unstranded_feature_width,  "mm"))
      } +
    { if(add.unstranded & add.feature ) geom_feature_label(
        data = gene.unstranded,
        aes(x = xmean, y = ymin, label=gene), forward = 1,
        feature_height = grid::unit(gene.sizes$unstranded_feature_label_height, "mm"),
        label_height   = grid::unit(gene.sizes$unstranded_feature_label_text,   "mm") #,size=gene.sizes[5]
          )
      } +


    ## CageFighteR results
    { if(add.cageTSS ) geom_gene_arrow(
      arrowhead_height  = grid::unit(gene.sizes$genome_feature_arrowhead_height,  "mm"),
      arrow_body_height = grid::unit(gene.sizes$genome_feature_arrow_body_height, "mm"),
      data=cagefr.clust,#[cagefr.clust$strand == '-', ],
      aes(xmin = region_start, xmax = region_end, y = y), fill="white")
    } +
    { if(add.cageTSS ) geom_subgene_arrow(
      arrowhead_height  = grid::unit(gene.sizes$genome_feature_arrowhead_height,  "mm"),
      arrow_body_height = grid::unit(gene.sizes$genome_feature_arrow_body_height, "mm"),
      data = cagefr.clust,#[cagefr.clust$strand == '-', ],
      aes(xmin = region_start, xmax = region_end, y = y, fill = strand, xsubmin = start, xsubmax = end, color=strand))
    } +


    ## Filling colors
    #scale_fill_manual(values=alpha(palette[c(2,3,1)],   alpha=alpha$gene_geom)) + #
    #scale_color_manual(values=alpha(palette[c(2,3,1)],   alpha)) +
    { if(length(unique(as.character(gene.plotsub$strand))) == 3 ) scale_fill_manual(values=palette[c(2,3,1)]) } +
    { if(length(unique(as.character(gene.plotsub$strand))) == 2 ) scale_fill_manual(values=palette[c(1,2,3)]) } +

    ## x-axis scale
    scale_x_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = seq(0, l_genome, breakseq) ) +
    ## y-axis scale
    scale_y_continuous(breaks = unique(gene.plotsub$ymin)) +

    ## theme components
    #theme_ipsum() +
    theme(axis.ticks.length.y = unit(0, "cm"),
          axis.text.y = element_blank(), axis.title.y = element_blank(), #element_text(size = 10)
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          panel.grid.minor.x = element_blank(), #panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
          legend.text  = element_text(size = sizes),
          legend.title = element_text(size = sizes),
          legend.position = 'none',
          plot.margin = margins
    )


    ## set axes' limits
    gggenome <- gggenome + coord_cartesian(ylim = c(ylims[1], ylims[2])
                                         , xlim = c(visfrom, visto) )
  }


  #### Genome coverage plot ####
  if(!genome.only) {

    plot.sub <- plot.data[plot.data$start >= visfrom & plot.data$end <= visto , ]
    plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-', '*'))

    if(!is.na(samples)) {
      plot.sub <- plot.sub[plot.sub$hpi == samples, ]
    }

    if(!is.na(comb.all.samples)) {
      plot.sub$hpi <- comb.all.samples
    }

    if (crop.FALSE) {
      if (prime == 'prime3') {
        plot.sub <- plot.sub[plot.sub$correct_tes ==T, ]
        gg.prime  <- ggplot(plot.sub, aes(x=prime3, fill=strand,  color=strand ))

      } else if (prime == 'prime5') {
        plot.sub <- plot.sub[plot.sub$correct_tss ==T, ]
        gg.prime  <- ggplot(plot.sub, aes(x=prime5, fill=strand,  color=strand ))

      } else if (prime == 'both')   {
        plot.sub <- plot.sub[plot.sub$correct_tss ==T
                             & plot.sub$correct_tes ==T, ]
        gg.prime  <- ggplot(plot.sub, aes(x=pos,    fill=endtype, color=endtype))

      }
      gg.prime <- gg.prime + facet_cropF

    } else {
      if (is.na(prime))   {
        plot.sub <- plot.sub[,]
        gg.prime  <- ggplot(plot.sub, aes(x=pos, fill=strand, color=strand)) + #, color=strand
          facet_wrap(~ hpi, scales = scales, strip.position = 'right', ncol=1)  #### facet_grid(rows = vars(hpi, correct_tss), scales = scales)
      } else if (prime == 'prime3') {
        plot.sub <- plot.sub[, ]
        gg.prime  <- ggplot(plot.sub, aes(x=prime3, fill=strand, color=strand)) +
          facet_wrap(~ hpi + correct_tes, scales = scales, strip.position = 'right', ncol=1)  #### facet_grid(rows = vars(hpi, correct_tes), scales = scales)
      } else if (prime == 'prime5') {
        plot.sub <- plot.sub[,]
        gg.prime  <- ggplot(plot.sub, aes(x=prime5, fill=strand, color=strand)) + #, color=strand
          facet_wrap(~ hpi + correct_tss, scales = scales, strip.position = 'right', ncol=1)  #### facet_grid(rows = vars(hpi, correct_tss), scales = scales)
      } else if (prime == 'both')   {
        plot.sub <- plot.sub[,]
        gg.prime  <- ggplot(plot.sub, aes(x=pos, fill=endtype, color=endtype)) + #, color=strand
          facet_wrap(~ hpi, scales = scales, strip.position = 'right', ncol=1)  #### facet_grid(rows = vars(hpi, correct_tss), scales = scales)
      }
    }


    gg.prime <- gg.prime +
      geom +
      { if(!is.null(geom2)) geom2 } +
      scale_fill_manual(values=palette) +
      #scale_color_manual(values=alpha(palette, alpha)) +
      scale_x_continuous(breaks = seq(0, l_genome, breakseq),  # c(0, 10000, 50000, 75000, 100000, 125000
                         labels=function(x) format(x, big.mark = " ", scientific = FALSE) ) +
      geom_vline(xintercept = vline,  linetype='dashed', color='blue', size=1) +
      geom_vline(xintercept = vline2, linetype='dashed', color='blue', size=1) +
      #theme_ipsum() +
      theme(strip.text =  strip.text, #axis.text.y = element_blank(),
            strip.background = element_rect(fill='lightgrey'),
            strip.placement = 'outside' ,
            plot.subtitle = element_text(size=sizes),
            axis.text.x = element_text(size=sizes), axis.ticks.length.y = unit(0, "cm"),
            axis.text.y = element_text(size=sizes, margin = margin(r = 0, l=0)),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.minor.x = element_blank(), #panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(), #panel.grid.major.y = element_blank(),
            legend.text  = element_text(size = sizes),
            legend.title = element_text(size = sizes + 2),
            plot.margin = margins,
            ...)

    if(!is.na(plot.title)) {
      gg.prime <- gg.prime + labs(title=plot.title)
    }

    if (!is.null(ylim)) { gg.prime <- gg.prime + coord_cartesian(xlim = c(visfrom, visto), ylim = ylim )
    } else {
      gg.prime <- gg.prime + coord_cartesian(xlim = c(visfrom, visto)) # xlim(c(visfrom, visto))
    }

    if(y.log10)       {gg.prime <- gg.prime + scale_y_log10() }
  }

  #### Assemble plot ####
  if(genome.only) {
    plot.exp <- gggenome
  } else {
    if(add.genome.plot) {
      plot.exp <- cowplot::plot_grid(gg.prime, gggenome, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(genomplot.scale, 1))
    } else {
      plot.exp <- gg.prime
    }
  }

  plot.exp
}
