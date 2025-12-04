#' Visualize 5 or 3-primes on a genome
#'
#' @export

plot.genome.region <-
  function(visfrom, visto, genome=NA, plot.data, gene.plotdata=NA, TR.merged.data=NA,
           geom = geom_bar(), geom2 = NULL,
           prime='prime5', crop.FALSE=F, samples=NA, comb.all.samples=NA,
           sum.counts.in.window=T, bin_width = 200, sum.fun='sum', add.all.pos=F, y.thresh=0, minus.strand.down=F,

           #what.to.plot = c('coverage', 'gene annotation', 'transcripts'),

           add.coverage=T, add.transcripts.plot=F, add.genome.plot=T, flip.panels=F, genome.only=F, transcripts.only=F,
           genome.and.transcripts=F, transcript.plot.scale=3,
           facet_genes = facet_nested(rows=vars(seqnames), scales=scales, drop=T),

           add.TR.labels=T, asterisk.size=4, tr.size.multip = 1, add.CAGE_significance.to.TRs =T, TR.y.add=1,
           facet_TR = facet_nested(rows=vars(seqnames), scales=scales, drop=T),

           ## CAGE (or other) clusters
           add.cageTSS=F, cagefr.clust=NULL, cagefr.clust.dist=6,
           cluster.aes = aes(xmin = region_start, xmax = region_end, y = y, forward = orientation, fill = strand, xsubmin = start, xsubmax = end),
           add.cageTSS.features = F,

           gene.label=T, add.unstranded=T, add.feature=F, gene.feature.width.thresh=700, flip.gene.y=T, force.gene.y=T, force.gene.y.all = F, 
           force.all.gene.down=F, angle=-60, force.one.lab.per.gene=T, force.all.gene.up.and.down=F,
           transcript.sizes=NULL,
           annot_fill_column='strand', tr.palette = NULL,

           sizes=16, gene.sizes=NULL, gene.label.col='white', palette=rep("#3CBC75D9", 10), alpha=0.8, plot.title=NA, # pal_npg()(10)
           gene_name_col = 'gene', tolower=F, vline=NULL, vline2=NULL, breakseq=5000, n_breaks = 5, margins=unit(c(1,1,1,1), "cm"),
           scales='fixed', genomplot.scale=5.5, ylim=c(0, NA), y.log10=F, y.log2=F, ylims.gene=NULL, y.multip=1.5,
           ybreaks=NULL, ybreak_n=10,
           strip.text =  element_text(angle = 45, size = sizes*1.5, hjust = 0.5), facet_space=1,
           gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene),
           facet_cropF = facet_nested(rows = vars(hpi), cols=NULL, scales ),
           labels=NULL, legend.position.prime = 'top', theme_general = theme_ipsum(),
           xlab_name=NULL, ylab_name=NULL,
           return.plot.data=F,
           ...) {


    
  ## Adaptive breakseq
  if(is.null(breakseq) ) {
    
    d <- (visto - visfrom) / n_breaks
    
    s <- round( log10(d), 0)
    
    f <- cumsum(rep(10^(s)/2, 20))
    
    breakseq <- f[abs(f -d) == min(abs(f -d))]
    
  } else { breakseq <- breakseq }
  
  message('OK')
  
  #### Vlines
    
  vline2 <- if (!is.null(vline2)) {
    list(
      new_scale_colour(),                             # must come before the vline
      vline2,
      scale_colour_manual(values = palette[-(1:3)])
    )
  } else {
    list()                                           # empty list if nothing to add
  }
  
  vline  <- if (!is.null(vline)) {
    list(
      new_scale_colour(),                            # must come before the vline
      vline,
      scale_colour_manual(values = palette[-(1:3)])
    )
  } else {
    list()                                           # empty list if nothing to add
  }
  
  
    
  #### Genome and annotation plot ####

  ### subset the features
  visregion      <- data.frame(seqnames=genome, start=visfrom, end=visto)

  #### Transcript annotation

  if (add.transcripts.plot) {

    ## Starting from TR.merged.data

    TR.subdata <- make_read_plot_data(DT = TR.merged.data, #[transcript_id %in% tr_toplot,],
                                      prime3_cluster_window = 25, cluster_neighbours = T, optimal_overlap = T,
                                      add.introns = T, is.transcript.pos=T,
                                      genome = genome, visfrom = visfrom, visto=visto, strand.toplot = NULL,
                                      add.class_code = '', flip_minus = T)



    TR.subdata$strand <- factor(TR.subdata$strand, levels = c('+', '-', '*'))
    tr.size.multip
    ## Filling colors
    if(is.null(tr.palette)) { tr.palette <- palette }
    
    ## Annotation sizes
    if (is.null(transcript.sizes)) {
      transcript.sizes <- data.frame(
        transcript.label.size = 1.25,
        exon.rect_height = 0.1, exon.rect.linewidth=0.01,
        intron.linetype = 'dashed', intron.linewidth = 0.1)

    }
    exon.rect_height      <- transcript.sizes$exon.rect_height # 0.75
    exon.rect.linewidth   <- transcript.sizes$exon.rect.linewidth # 0.2
    intron.linetype       <- transcript.sizes$intron.linetype # 'dashed'
    intron.linewidth      <- transcript.sizes$intron.linewidth # 0.1
    transcript.label.size <- transcript.sizes$transcript.label.size

    TR.subdata$size <- exon.rect.linewidth /  10000

    ## Axis limit
    ylim.tr    <- c(TR.subdata[, min(ypos) ] - TR.y.add,
                    TR.subdata[, max(ypos) ] + TR.y.add   )

    if(is.null(visfrom) & is.null(visto)) {
      xlim       <- c(TR.subdata[, min(start.exon) ],
                      TR.subdata[, max(end.exon) ]   )
    } else {
      xlim       <- c(visfrom, visto)
    }

    CAGE.nudge_x <- (visregion$end - visregion$start) / 100

    ggtr <- ggplot(TR.subdata) +
      
      ## Gene arrows for the first exons
      geom_gene_arrow(data=TR.subdata[last_exon == T,],
                      aes(fill = as.factor(!!sym((annot_fill_column))),
                          xmin = start.exon, xmax = end.exon,
                          y    = ypos, forward = orientation,
                          #ymin=ypos-0.2, ymax=ypos+0.2
                      ),
                      arrowhead_height  = grid::unit(gene.sizes$gene_arrowhead_height  * tr.size.multip,  "mm"),
                      arrow_body_height = grid::unit(gene.sizes$gene_arrow_body_height * tr.size.multip, "mm"),
                      alpha=alpha$gene_geom) +

      ## Geom rect for other exons
      geom_gene_rect(data=TR.subdata[last_exon == F,],
                     aes(fill = as.factor(!!sym((annot_fill_column))),
                         xmin = start.exon, xmax = end.exon,
                         #ymin = ypos - (exon.rect_height ), ymax = ypos + (exon.rect_height )
                         y = ypos
                         #, forward = orientation
                         #,size=size
                         #,color=genome
                     ),
                     #linewidth=exon.rect.linewidth,
                     #arrowhead_height  = grid::unit(gene.sizes$gene_arrowhead_height,  "mm"),
                     #arrow_body_height = grid::unit(gene.sizes$gene_arrow_body_height, "mm"),
                     rect_height = grid::unit(gene.sizes$gene_arrow_body_height * tr.size.multip, "mm"),
                     color='black',
                     alpha=alpha$gene_geom) +
      
      
      ## Introns
      geom_segment(data=TR.subdata[exon_number != 1,],
                   aes(x=intron_start-1, xend=intron_end+1, y=ypos, yend=ypos),
                   color='black', linetype=intron.linetype, linewidth=intron.linewidth) +
      
      ## Transcript labels
      { if (add.TR.labels)
        geom_gene_label(data = TR.subdata[last_exon == T,],
                        aes(xmin=start.exon, xmax=end.exon, y = ypos, label=transcript_id)
                        # , forward = 1
                        ,height = grid::unit(transcript.label.size, "mm"), grow = F, colour=gene.label.col, align = "centre"
                        #feature_height = grid::unit(linewidth, "mm"),
                        #feature_width  = grid::unit(linewidth, "mm")
        ) } +
      
      ## Label CAGE supported transcripts with stars
      { if (add.CAGE_significance.to.TRs) 
        geom_text(aes(x=ifelse(strand == '+',
                               prime5.TR - CAGE.nudge_x,
                               prime5.TR + CAGE.nudge_x),
                      y=ypos,
                      label = ifelse(is.na(CAGE_significance), '', CAGE_significance)),
                  size = asterisk.size,
                  fontface = "bold"
        ) } +
      
      
      ## Vertical lines
      { if(!is.null(vline ))  vline   } +
      { if(!is.null(vline2))  vline2  } +
      
      #if(!is.null(vline )) { ggtr <- ggtr + new_scale_colour() + vline  + scale_colour_manual(values = palette[-c(1:3)]) }
      
      
      ## Axes and coordinates
      scale_y_continuous(breaks = c(-5, 0, 5, 10, 15)) +
      coord_cartesian(
        xlim=xlim,    # c(visfrom, visto), #
        ylim=ylim.tr  # c(1530, 1590)
      ) +
      
      ## X-axis scale
      scale_x_continuous(labels  = function(x) format(x, big.mark = " ", scientific = FALSE),
                         #limits  = c(xlim[1], xlim[2]),
                         breaks  = seq(0, l_genome, breakseq) ) +
      
      ## Theme elements
      #labs(x='Genomic Position', y='Transcripts (grouped by 3-prime end)', title='Transcript Visualization') +
      theme_general +
      theme(
        strip.text       = strip.text,
        strip.background = element_rect(fill='lightgrey'),
        strip.placement  = 'outside' ,
        
        axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.x = element_blank(),
        
        ## Remove Y-axis ticks and lines
        axis.ticks.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        
        legend.position = 'none',
        ...
        
      ) +
      
      ## Plot margins
      
      # If margins is set explicitly, it will apply to each part of the plot
      { if (!is.null(margins))   theme(plot.margin = margins) } +
      
      # Default is when the "margins" parameter is NULL.
      # In this case the spaces between the genome and transcriptome and "coverage" plots are minimized.
      { if (is.null(margins) & !transcripts.only) theme(plot.margin = unit(c(5,5,1,5), 'mm')) } +
      
      ## Faceting -->> ???
      facet_TR
    #NULL
    
    ## Filling colors for strand or other category
    
    if(annot_fill_column == 'strand') {
      strands <- factor(unique(as.character(TR.subdata$strand)), levels = c('+', '-', '*'))
      strands <- strands[order(strands)]
      
      if        ( all(is.element(c('+', '-'), strands)) )     {
        ggtr <- ggtr + scale_fill_manual(values=alpha(tr.palette[c(1,2)],   alpha=alpha$gene_geom))
        
      } else if ( all(is.element(c('+', '*'), strands)) )     {
        ggtr <- ggtr + scale_fill_manual(values=alpha(tr.palette[c(1,3)],   alpha=alpha$gene_geom))
        
      } else if ( all(is.element(c('-', '*'), strands))  )    {
        ggtr <- ggtr + scale_fill_manual(values=alpha(tr.palette[c(2,3)],   alpha=alpha$gene_geom))
        
      } else if ( all(is.element(c('+', '-', '*'), strands)) ){
        ggtr <- ggtr + scale_fill_manual(values=alpha(tr.palette[c(1,2,3)], alpha=alpha$gene_geom))
        
      } else if ( all(is.element(c('+'), strands)) )     {
        ggtr <- ggtr + scale_fill_manual(values=alpha(tr.palette[c(1, 2)],   alpha=alpha$gene_geom))
        
      } else if ( all(is.element(c('-'), strands))  )    {
        ggtr <- ggtr + scale_fill_manual(values=alpha(tr.palette[c(2, 1)],   alpha=alpha$gene_geom))
        
      }
    } else {
      ggtr <- ggtr + scale_fill_manual(values=alpha(tr.palette[],   alpha=alpha$gene_geom))
    }
    
    # ggtr
    # ggsave('test.jpg', ggtr, width = 9, height=15)
    
    
    
    ## set axes' limits fo Gene annotation also
    visfrom  <- xlim[1]
    visto    <- xlim[2]
    
    
    ### Update visregion
    visregion      <- data.frame(seqnames=genome, start=visfrom, end=visto)
    
  }
  
  #### Gene annotation
  #if(add.genome.plot) {
  
  gene.plotsub   <- dplyr::select(
    genome_join(gene.plotdata, visregion, by=c('seqnames', 'start', 'end')),
    -c(seqnames.y, start.y, end.y) )
  colnames(gene.plotsub)[is.element(colnames(gene.plotsub), c('seqnames.x', 'start.x', 'end.x'))] <- c('seqnames', 'start', 'end')
  
  ## make DT
  setDT(gene.plotsub)
  
  ## factorize strand
  gene.plotsub[, strand := factor(strand, levels=c('+', '-', '*'))]
  
  ## facet
  gene.plotsub[, hpi  := genome]
  
  
  ## put 'clusters' into one y-position
  for (i in unique(gene.plotsub$cluster)) {
    gene.plotsub$ymin[gene.plotsub$cluster == i] <- seq(1, nrow(gene.plotsub[gene.plotsub$cluster == i, ]))
    
  }
  #gene.plotsub[, ymin := seq(1, n_distinct(start)), by=cluster]
  
  ## back to DF
  setDF(gene.plotsub)
  
  ### Add sub-gene positions
  gene.plotsub$subgene.start <- gene.plotsub$start
  gene.plotsub$subgene.end   <- gene.plotsub$end
  for (i in unique(gene.plotsub$gene)) {
    gene.plotsub$start[gene.plotsub$gene == i] <- min(gene.plotsub$start[gene.plotsub$gene == i])
    gene.plotsub$end  [gene.plotsub$gene == i] <- max(gene.plotsub$end  [gene.plotsub$gene == i])
  }
  gene.plotsub$width <- abs(gene.plotsub$end - gene.plotsub$start)+1
  
  ### Gene names
  gene.plotsub[,'gene_name'] <- gene.plotsub[,gene_name_col]
  ## make gene names lowercase
  if(tolower) { gene.plotsub$gene_name <- tolower(gene.plotsub$gene_name) }
  
  ## ORi regions strand !!!
  #gene.plotsub$strand[grep('Ori', gene.plotsub$gene)] <- '*' #### !!!!
  
  ## flip
  if (flip.gene.y)  { gene.plotsub$ymin <- gene.plotsub$ymin * -1 }
  
  
  gene.plotsub$ylab <- gene.plotsub$ymin
  
  
  ## force into one row BY STRAND
  if (force.gene.y) {
    gene.plotsub$ymin[gene.plotsub$strand=='+'] <- 1
    gene.plotsub$ymin[gene.plotsub$strand=='-'] <- 0
    gene.plotsub$ymin[gene.plotsub$strand=='*'] <- 0.5
  }
  
  ## stranded and unstranded features for geom_gene arrow and geom_rect
  gene.unstranded <- gene.plotsub[gene.plotsub$strand == '*', ] #; gene.unstranded$strand <- factor(gene.unstranded$strand, levels=  '*')
  gene.stranded   <- gene.plotsub[gene.plotsub$strand != '*', ] #; gene.stranded$strand   <- factor(gene.stranded$strand,   levels=c('+', '-'))
  
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
  
  if (force.all.gene.up.and.down) {
    gene.feature.upwards   <- gene.stranded$gene[gene.stranded$strand == '+']
    gene.feature.downwards <- gene.stranded$gene[gene.stranded$strand == '-']
  }
  
  
  ## add unstranded !!!
  if (add.unstranded) {
    #gene.unstranded$start[gene.unstranded$width < 100] <- gene.unstranded$start - 20
    #gene.unstranded$end  [gene.unstranded$width < 100] <- gene.unstranded$end   + 20
    
    ## force unstranded y position
    try({ gene.unstranded$ymin <- max(gene.plotsub$ymin) })
  }
  
  ## force giving only one label for each gene (in case of multi-exoned genes)
  if (force.one.lab.per.gene) {
    multi.genes <- dup(gene.stranded$gene)
    for (mgene in multi.genes) {
      gene.stranded[gene.stranded$gene == mgene & !is.na(gene.stranded$gene),'gene_name'][1]  <- mgene
      gene.stranded[gene.stranded$gene == mgene & !is.na(gene.stranded$gene),'gene_name'][-1] <- NA
    }
  }
  
  ## recalculate xlim of plots to include gene annotations coverage
  if(!exists('xlim')) {
    visfrom    <- min(c(visfrom, min(gene.plotsub$start)))
    visto      <- max(c(visto,   max(gene.plotsub$end  )))
    visregion  <- data.frame(seqnames=genome, start=visfrom, end=visto)
  }
  
  gene.plotsub$strand    <- factor(gene.plotsub$strand,    levels=c('+', '-', '*'))
  gene.stranded$strand   <- factor(gene.stranded$strand,   levels=c('+', '-', '*'))
  gene.unstranded$strand <- factor(gene.unstranded$strand, levels=c('+', '-', '*'))
   
  ### Check if there are CAGE clusters in the region
  if( add.cageTSS ) {
    cagefr.clust <- dplyr::select(genome_join(cagefr.clust, visregion, by=c('seqnames', 'start', 'end')), -c(seqnames.y, start.y, end.y) )
    colnames(cagefr.clust)[is.element(colnames(cagefr.clust), c('seqnames.x', 'start.x', 'end.x'))] <- c('seqnames', 'start', 'end')
    cagefr.clust$strand <- factor(cagefr.clust$strand, levels=c('+', '-', '*'))
    
    if (nrow(cagefr.clust) > 0) { add.cageTSS <- T } else { add.cageTSS <- F }
    
  }
  
  ### CAGE TSS clusters as genome feature
  if(add.cageTSS ) {
    
    ## faceting
    cagefr.clust$hpi <- genome
    
    ## y-positions
    #cagefr.clust$ymin  <-  max(gene.plotsub$ymin) + cagefr.clust.dist
    #cagefr.clust$ymax  <-  max(gene.plotsub$ymin) + cagefr.clust.dist + 1
    #cagefr.clust$y[cagefr.clust$strand == '-'] <- cagefr.clust$ymin
    #cagefr.clust$y[cagefr.clust$strand == '+'] <- cagefr.clust$ymax
    
    # Adjust 'y' values for '+' strand based on 'dist'
    cagefr.clust$y[cagefr.clust$strand == '+'] <-
      ifelse(cagefr.clust.dist > 0,
             max(gene.plotsub$ymin) + cagefr.clust.dist + 1,
             min(gene.plotsub$ymin) - cagefr.clust.dist)
    
    # Adjust 'y' values for '-' strand based on 'dist'
    cagefr.clust$y[cagefr.clust$strand == '-'] <-
      ifelse(cagefr.clust.dist > 0,
             max(gene.plotsub$ymin) + cagefr.clust.dist,
             min(gene.plotsub$ymin) - cagefr.clust.dist - 1)
    
    cagefr.clust$ymin  <-  min(cagefr.clust$y)
    cagefr.clust$ymax  <-  max(cagefr.clust$y)
    
    # crop other clusters
    #cagefr.clust$region_start <- visfrom # min(cagefr.clust$start)
    #cagefr.clust$region_end   <- visto   # max(cagefr.clust$end)
    cagefr.clust$visfrom <- visfrom # min(cagefr.clust$start)
    cagefr.clust$visto   <- visto   # max(cagefr.clust$end)
    
    
    ## factorize strands
    cagefr.clust$strand <-  factor(cagefr.clust$strand, levels = c('+' ,'-', '*'))
    
    ## y-breaks
    ybreaks.genome <- unique(c(gene.plotsub$ymin, cagefr.clust$y))
    
  } else {
    cagefr.clust   <- data.frame(ymin=min(gene.plotsub$ymin), ymax=max(gene.plotsub$ymin))
    
    ybreaks.genome <- unique( c(gene.plotsub$ymin  ))
    
  }
  
  
  ## force EVERY ANNOTATION into one row
  gene.plotsub$ylab <- gene.plotsub$ymin
  if (force.gene.y.all) {
    gene.stranded$ymin    <- 0
    gene.unstranded$ymean <- 0
    
    message('all annotations were forced into one row!')
  }
  
  
  ### Set y-axis limits
  if (is.null(ylims.gene)) {
    ylims.gene <- c(min(c(gene.plotsub$ymin, cagefr.clust$ymin)) - y.multip,
                    max(c(gene.plotsub$ymin, cagefr.clust$ymax)) + y.multip)
  } else {
    ylims.gene <- ylims.gene
  }
  
  
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
    
    ## Gene Arrow labels or Gene Feature Labels?
    { if(gene.label & !add.feature & !force.all.gene.down & !force.all.gene.up.and.down)
      message('Adding Gene Arrow Labels')
      geom_gene_label(
        data  = gene.stranded[!is.element(gene.stranded$gene, c(gene.feature.upwards, gene.feature.downwards)),],
        height = grid::unit(gene.sizes$gene_label_height, "mm"), grow = F, colour=gene.label.col, align = "centre")
    } +
    
    ## Gene feature labels --> pointing upwards
    { if(gene.label & add.feature  & !force.all.gene.down & !force.all.gene.up.and.down)
      message('Adding Gene Feature Labels')
      geom_feature(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.upwards),  ],
        aes(x = xmean, y = ymin), forward = 1,
        feature_height = grid::unit(gene.sizes$gene_feature_height, "mm"),
        feature_width  = grid::unit(gene.sizes$gene_feature_width,  "mm")
      ) } +
    
    { if(gene.label & add.feature  & !force.all.gene.down & !force.all.gene.up.and.down)
      geom_feature_label(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.upwards), ],
        aes(x = xmean, y = ymin, label=gene_name), forward = 1,
        feature_height = grid::unit(gene.sizes$gene_feature_label_height, "mm"),
        label_height   = grid::unit(gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
      ) } +
    
    ## Gene feature labels  --> pointing downwards
    { if(gene.label & add.feature  & !force.all.gene.down & !force.all.gene.up.and.down)
      geom_feature(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
        aes(x = xmean, y = ymin), forward = 1,
        feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
        feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm")
      ) } +
    
    { if(gene.label & add.feature  & !force.all.gene.down & !force.all.gene.up.and.down)
      geom_feature_label(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
        aes(x = xmean, y = ymin, label=gene_name), forward = 1,
        feature_height = grid::unit(-gene.sizes$gene_feature_label_height, "mm"),
        label_height   = grid::unit(gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
      ) } +
    
    ## Force all gene features pointing downwards
    { if(gene.label & add.feature & force.all.gene.down ) geom_feature(
      data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      aes(x = xmean, y = ymin), forward = 1,
      #size=gene.sizes$gene_label_height,
      feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
      feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm"))
    } +
    { if(gene.label & add.feature & force.all.gene.down )
      #geom_feature_label(
      #data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      #aes(x = xmean, y = ymin, label=gene, angle=angle), forward = 1,
      ##size=gene.sizes$gene_label_height,
      #feature_height = grid::unit(-gene.sizes$gene_feature_label_height, "mm"),
      #label_height   = grid::unit( gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
      #)
      geom_text(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
        aes(x = xmean, y = ymin-gene.sizes$gene_feature_label_height,
            label=gene_name, angle=angle),
        size=gene.sizes$gene_label_height)
    } +
    
    ## Force (+) gene features pointing upwards and (-) downwards
    # Down
    { if(gene.label & add.feature & force.all.gene.up.and.down )
      geom_feature(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
        aes(x = xmean, y = ymin), forward = 1,
        #size=gene.sizes$gene_label_height,
        feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
        feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm"))
    } +
    { if(gene.label & add.feature & force.all.gene.up.and.down )
      #geom_feature_label(
      #data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      #aes(x = xmean, y = ymin, label=gene, angle=angle), forward = 1,
      ##size=gene.sizes$gene_label_height,
      #feature_height = grid::unit(-gene.sizes$gene_feature_label_height, "mm"),
      #label_height   = grid::unit( gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
      #)
      geom_text(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
        aes(x = xmean, y = ymin-gene.sizes$gene_feature_label_height,
            label=gene_name, angle=angle),
        size=gene.sizes$gene_label_height)
    } +
    # Up
    { if(gene.label & add.feature & force.all.gene.up.and.down )
      geom_feature(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.upwards),],
        aes(x = xmean, y = ymin), forward = 1,
        #size=gene.sizes$gene_label_height,
        feature_height = grid::unit(gene.sizes$gene_feature_height, "mm"),
        feature_width  = grid::unit(gene.sizes$gene_feature_width,  "mm"))
    } +
    { if(gene.label & add.feature & force.all.gene.up.and.down )
      #geom_feature_label(
      #data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
      #aes(x = xmean, y = ymin, label=gene, angle=angle), forward = 1,
      ##size=gene.sizes$gene_label_height,
      #feature_height = grid::unit(-gene.sizes$gene_feature_label_height, "mm"),
      #label_height   = grid::unit( gene.sizes$gene_feature_label_text,   "mm") #,size=gene.sizes[5]
      #)
      geom_text(
        data = gene.stranded[is.element(gene.stranded$gene, gene.feature.upwards),],
        aes(x = xmean, y = ymin+gene.sizes$gene_feature_label_height,
            label=gene_name, angle=angle),
        size=gene.sizes$gene_label_height)
    } +
    
    #gggenome +
    ## Unstranded features with feature label
    { if(add.unstranded ) geom_rect(
      data=gene.unstranded,
      aes(xmin=start, xmax=end,
          ymin=ymean - gene.sizes$unstranded_rect_height/2,
          ymax=ymean + gene.sizes$unstranded_rect_height/2,
          
          #ymin=ymin-gene.sizes$unstranded_rect_height/2,
          #ymax=ymin+gene.sizes$unstranded_rect_height/2,
          fill=strand), color='black', fill='white')
    } +
    { if(add.unstranded ) geom_rect(
      data=gene.unstranded,
      aes(xmin=subgene.start, xmax=subgene.end,
          ymin=ymean - gene.sizes$unstranded_rect_height/2,
          ymax=ymean + gene.sizes$unstranded_rect_height/2,
          
          #ymin=ymin-gene.sizes$unstranded_rect_height/2,
          #ymax=ymin+gene.sizes$unstranded_rect_height/2,
          fill=strand), color='black', alpha=alpha$unstranded_geom)
    } +
    { if(add.unstranded & add.feature ) geom_feature(
      data = gene.unstranded,
      aes(x = xmean, y = ymin), forward = 1,
      feature_height = grid::unit(-gene.sizes$unstranded_feature_height, "mm"),
      feature_width  = grid::unit(-gene.sizes$unstranded_feature_width,  "mm"))
    } +
    { if(add.unstranded & add.feature ) geom_feature_label(
      data = gene.unstranded,
      aes(x = xmean, y = ymin, label=gene_name), forward = 1,
      size = gene.sizes$unstranded_feature_label_text,
      feature_height = grid::unit(gene.sizes$unstranded_feature_label_height, "mm"),
      label_height   = grid::unit(gene.sizes$unstranded_feature_label_text,   "mm") #,size=gene.sizes[5]
    )
    } +
    #{ if(add.unstranded ) geom_text(
    #  data = gene.unstranded,
    #  aes(x = xmean, y = ymin-gene.sizes$gene_feature_label_height,
    #      label=gene_name, angle=angle),
    #  size=gene.sizes$gene_label_height)
    #} +
    
    
    ## CAGE results
    { if( add.cageTSS ) geom_gene_arrow(
      arrowhead_height  = grid::unit(gene.sizes$genome_feature_arrowhead_height,  "mm"),
      arrow_body_height = grid::unit(gene.sizes$genome_feature_arrow_body_height, "mm"),
      data = cagefr.clust,      #[cagefr.clust$strand == '-', ],
      aes(xmin = visfrom, xmax = visto, y = y, forward = orientation), fill="white"
    ) } + 
    { if( add.cageTSS ) geom_gene_arrow(
      arrowhead_height  = grid::unit(gene.sizes$genome_feature_arrowhead_height,  "mm"),
      arrow_body_height = grid::unit(gene.sizes$genome_feature_arrow_body_height, "mm"),
      data = cagefr.clust,    #[cagefr.clust$strand == '-', ],
      mapping = cluster.aes 
    ) } + 
    { if( add.cageTSS ) geom_feature(
      feature_height  = unit(3, "mm"),
      feature_width   = unit(3, "mm"),
      arrowhead_width = unit(2, "mm"),
      data  = cagefr.clust[cagefr.clust$strand == '+', ],
      aes(x = start, y = y, forward = orientation)
    ) } +
    { if( add.cageTSS & add.cageTSS.features) geom_feature_label(
      feature_height  = unit(3, "mm"),
      feature_width   = unit(3, "mm"),
      arrowhead_width = unit(2, "mm"),
      data  = cagefr.clust[cagefr.clust$strand == '+', ],
      aes(x = start, y = y, label = gene, forward = orientation)
    ) } +
    { if( add.cageTSS ) geom_feature(
      feature_height  = unit(-3, "mm"),
      feature_width   = unit(-3, "mm"),
      arrowhead_width = unit(-2, "mm"),
      data  = cagefr.clust[cagefr.clust$strand == '-', ],
      aes(x = start, y = y, forward = orientation)
    ) } +
    { if( add.cageTSS & add.cageTSS.features) geom_feature_label(
      feature_height  = unit(-3, "mm"),
      feature_width   = unit(-3, "mm"),
      arrowhead_width = unit(-2, "mm"),
      data  = cagefr.clust[cagefr.clust$strand == '-', ],
      aes(x = start, y = y, label = gene, forward = orientation)
    ) } +
    
    
    ## X-axis scale
    scale_x_continuous(labels  = function(x) format(x, big.mark = " ", scientific = FALSE),
                       #limits  = c(visfrom, visto),
                       breaks  = seq(0, l_genome, breakseq) ) +
    ## Y-axis scale
    scale_y_continuous(labels   = function(x) format(x, big.mark = " ", scientific = FALSE),
                       #limits   = ylims.gene,
                       breaks   = ybreaks.genome,
                       sec.axis = sec_axis(~.,
                                           labels = function(x) format(x, big.mark = " ", scientific = FALSE),
                                           breaks = ybreaks.genome
                       )) +
    
    #if(add.cageTSS) unique(c(gene.plotsub$ymin, cagefr.clust$y)) else unique(c(gene.plotsub$ymin  ))
    #{ if( add.cageTSS ) scale_y_continuous(breaks = unique(c(gene.plotsub$ymin, cagefr.clust$y))) } +
    #{ if(!add.cageTSS ) scale_y_continuous(breaks = unique(c(gene.plotsub$ymin  ))) } +
    
    
    ## Vertical lines
    { if(!is.null(vline ))  vline  } +
    { if(!is.null(vline2))  vline2 } +
    
    
    ## X and Y axis limits
    coord_cartesian(ylim = c(ylims.gene[1], ylims.gene[2]),
                    xlim = c(visfrom, visto) ) +
    
    #xlim(c(visfrom, visto)) +
    
    ## theme components
    theme_general +
    theme(strip.text       = strip.text,
          strip.background = element_rect(fill='lightgrey'),
          strip.placement  = 'outside' ,
          
          axis.text.x  = if (genome.and.transcripts | genome.only | transcripts.only) { element_text(size = sizes) } else { element_blank() },
          axis.text.y  = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), #element_text(size = 10)
          
          axis.ticks.length.y = unit(0, "cm"),
          
          panel.grid.minor.x = element_blank(),
          #panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          
          legend.position = 'none',
          legend.text     = element_text(size = sizes),
          legend.title    = element_text(size = sizes + 5),
          
          panel.spacing = grid::unit(facet_space, "mm"),
          
          plot.title    = element_text(size=sizes*2),
          plot.subtitle = element_text(size=sizes)
          
          #,...
    ) +
    
    
    ## Plot margins
    
    # If margins is set explicitly, it will apply to each part of the plot
    { if (!is.null(margins))   theme(plot.margin = margins) } +
    
    # Default is when the "margins" parameter is NULL.
    # In this case the spaces between the genome and transcriptome and "coverage" plots are minimized.
    { if (is.null(margins) & !genome.only)
      theme(plot.margin = unit(c(1,5,5,5), 'mm')) }
  
    ## --- Filling colors for gene annotations -----------------------------
    # Determine which strand labels appear
    present_strands <- unique(c(
      as.character(gene.stranded$strand),
      as.character(gene.unstranded$strand),
      if (exists("cagefr.clust") && "strand" %in% names(cagefr.clust))
        as.character(cagefr.clust$strand)
    ))
    present_strands <- intersect(c("+","-","*"), present_strands)
    
    # Build a named vector of fill colors with alpha transparency
    fill_vals <- c(
      "+" = scales::alpha(palette[1], alpha = alpha$gene_geom),
      "-" = scales::alpha(palette[2], alpha = alpha$gene_geom),
      "*" = scales::alpha(palette[3], alpha = alpha$gene_geom)
    )
    
    # Apply the manual fill scale
    gggenome <- gggenome +
      scale_fill_manual(
        values = fill_vals[present_strands],
        limits = present_strands
      )
    

    ## Faceting --> this
    gggenome <- gggenome + facet_genes #facet_cropF


    ## Finished annotation plot
  #}



  #### Genome coverage plot ####
  if(transcripts.only | genome.only) { add.coverage <- F }
  if(add.coverage) {

    ### Manipulate data

    ## Filter for region
    plot.sub <- plot.data[plot.data$start >= visfrom & plot.data$end <= visto , ]
    strands  <- as.character(unique(plot.sub$strand))

    ## Filter for samples?
    if(!all(is.na(samples))) {
      plot.sub     <- plot.sub[is.element(plot.sub$hpi, samples), ]
    }

    ## Combine all samples to one?
    if(!is.na(comb.all.samples)) {
      plot.sub$hpi <- comb.all.samples
    }

    ## Add all X-positions (zeros) ?
    if(add.all.pos) {
      fullvisreg <- data.frame(strand=c('+', '-', '*'), start=rep(c(visfrom:visto), 3), end=rep(c(visfrom:visto), 3))

      # Use map to create a list of data frames
      fullvisreg.hpi  <- rbindlist(purrr::map(unique(plot.sub$hpi),
                                       ~{df_temp <- fullvisreg; df_temp$hpi <- .; df_temp}))

      fullvisreg.hpi.correct_prime <- data.frame(rbindlist(purrr::map(as.character(unique(plot.sub[,correct_prime])),
                                                    ~{df_temp <- fullvisreg.hpi; df_temp[,correct_prime] <- .; df_temp})))
      fullvisreg.hpi.correct_prime[,correct_prime] <- as.logical(fullvisreg.hpi.correct_prime[,correct_prime])

      plot.sub <- data.frame(merge(data.table(plot.sub),
                                   data.table(fullvisreg.hpi.correct_prime),
                                   by=c('strand', 'start', 'end', 'hpi', correct_prime), all=T) )
      plot.sub[,prime] <- plot.sub[,'start']
      plot.sub$count[is.na(plot.sub$count)] <- 0


    }

    ## filter for minimum y-value
    if (!is.na(y.thresh)) {
      plot.sub <- plot.sub[plot.sub[,'count'] >= y.thresh,]
    }

    ### Summarise counts in a window of size: bin_width for density-like area plotting
    if (sum.counts.in.window) {
      df <- plot.sub
      df$count[is.na(df$count)] <- 0
      df$position <- df$start
      df$coverage <- df$count
      df$sample <- df$hpi

      # Define bin width
      if(prime=='prime5') {correct_prime <- 'correct_tss'} else if(prime == 'prime3') {correct_prime <- 'correct_tes'}
      # Create bins based on position
      df$bin <- cut(df$position, breaks = seq(from = min(df$position), to = max(df$position), by = bin_width), include.lowest = TRUE)

      # Calculate sum of coverage for each bin for each strand and each sample
      if (sum.fun == 'sum') {
        sum_coverage <- tapply(df$coverage, INDEX = list(df$sample, df$strand, df$bin, data.frame(df)[,correct_prime]), FUN = sum)
      } else if (sum.fun== 'mean') {
        sum_coverage <- tapply(df$coverage, INDEX = list(df$sample, df$strand, df$bin, data.frame(df)[,correct_prime]), FUN = mean)
      }

      # Convert the result to a data frame
      df_sum <- as.data.frame.table(sum_coverage, responseName = "coverage")
      names(df_sum)[1:4] <- c("sample", "strand", "bin", correct_prime)

      # Convert bin to numeric for plotting
      df_sum$bin <- gsub('\\[', '(', df_sum$bin)
      df_sum$bin_start <- as.numeric(gsub("\\((.*),.*\\]", "\\1", df_sum$bin))
      df_sum$bin_end   <- as.numeric(gsub("\\]", "", gsub(".*,", "", df_sum$bin)))
      df_sum$position  <- apply(df_sum[,c('bin_start', 'bin_end')], 1, mean)
      #as.numeric(gsub("\\((.*),.*\\]", "\\1", df_sum$bin)) + bin_width / 2

      colnames(df_sum)[1] <- 'hpi'
      df_sum[,prime] <- df_sum[,'position']

      df_sum$count <- df_sum$coverage
      plot.sub <- df_sum
    }

    ### Log?
    if(y.log10)   {
      plot.sub$count <- log10(plot.sub$count)
      plot.sub$count[plot.sub$count == '-Inf'] <- 0
    } else if(y.log2)   {
      plot.sub$count <- log2(plot.sub$count)
      plot.sub$count[plot.sub$count == '-Inf'] <- 0
    }


    ### Minus strand downwards
    if(minus.strand.down)   {
      plot.sub$count[plot.sub$strand == '-'] <- plot.sub$count[plot.sub$strand == '-'] * -1
    }

    ## Keep only original strands
    plot.sub <- plot.sub[is.element(plot.sub$strand, strands), ]

    ## Add unstranded elements and factorize strands
    if (add.unstranded) {
        plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-', '*'))
      } else {
        plot.sub[is.element(plot.sub$strand, c('+', '-')), ]
        plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-'))
      }

    ## factorise samples
    if(all(is.na(samples))) {
      samples <- unique(plot.sub$hpi)
    }
    plot.sub$hpi <- factor(plot.sub$hpi, levels = samples)
    samples      <- levels(plot.sub$hpi)


    ## Build plot
    if (crop.FALSE == T) {
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

      } else if (prime == 'pos')   {
        gg.prime  <- ggplot(plot.sub, aes(x=pos,  fill=strand, color=strand))
      }
      gg.prime <- gg.prime + facet_cropF

    } else if (crop.FALSE == F) {
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
    } else {
      plot.sub  <- plot.sub[,]
      gg.prime  <- ggplot(plot.sub) + facet_cropF
    }


    gg.prime <- gg.prime +
      geom +
      { if(!is.null(geom2)) geom2 }
      
      
    ## Filling colors
    strands <- factor(unique(as.character(plot.sub$strand)), levels = c('+', '-', '*'))
    strands <- strands[order(strands)]
    
    if        ( all(is.element(c('+', '-'), strands)) )     {
      gg.prime <- gg.prime + scale_fill_manual(values=alpha(palette[c(1,2)],   alpha=alpha$gene_geom))
      
    } else if ( all(is.element(c('+', '*'), strands)) )     {
      gg.prime <- gg.prime + scale_fill_manual(values=alpha(palette[c(1,3)],   alpha=alpha$gene_geom))
      
    } else if ( all(is.element(c('-', '*'), strands))  )    {
      gg.prime <- gg.prime + scale_fill_manual(values=alpha(palette[c(2,3)],   alpha=alpha$gene_geom))
      
    } else if ( all(is.element(c('+', '-', '*'), strands)) ){
      gg.prime <- gg.prime + scale_fill_manual(values=alpha(palette[c(1,2,3)], alpha=alpha$gene_geom))
      
    } else if ( all(is.element(c('+'), strands)) )     {
      gg.prime <- gg.prime + scale_fill_manual(values=alpha(palette[c(1, 2)],  alpha=alpha$gene_geom))
      
    } else if ( all(is.element(c('-'), strands))  )    {
      gg.prime <- gg.prime + scale_fill_manual(values=alpha(palette[c(2, 1)],  alpha=alpha$gene_geom))
      
    }
    
    gg.prime <- gg.prime +
      
      #scale_fill_manual(values =alpha(palette, alpha$cov_geom)) +
      #scale_color_manual(values=alpha(palette, alpha$cov_geom)) +

      ## x-axis scale
      scale_x_continuous(name=xlab_name,
                         labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                         #limits =  c(visfrom, visto),
                         breaks = seq(0, l_genome, breakseq) ) +

      ## Y-axis scale
      scale_y_continuous(name=ylab_name,
                         labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                         { if (!is.null(ybreaks)) breaks = ybreaks },
                         sec.axis = sec_axis(~.,
                                             labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                                             { if (!is.null(ybreaks)) breaks = ybreaks }
                         )) +

      #{ if(!is.null(ybreaks)) scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), breaks = ybreaks) } +

      ## X and Y axis limits
      #xlim(c(visfrom, visto)) +
      coord_cartesian(ylim = ylim,
                      xlim = c(visfrom, visto) ) +


      ## vertical lines
      #{ if(!is.null(vline))  geom_vline(xintercept = vline,  linetype='dashed', color='blue', size=1) } +
      #{ if(!is.null(vline2)) geom_vline(xintercept = vline2, linetype='dashed', color='blue', size=1) } +
      { if(!is.null(vline ))  vline  } +
      { if(!is.null(vline2))  vline2 } +

      ## Theme elements
      theme_general +
      theme(strip.text       = strip.text,
            strip.background = element_rect(fill='lightgrey'),
            strip.placement  = 'outside' ,

            axis.text.x  = element_text(size=sizes),
            axis.text.y  = element_text(size=sizes, margin = margin(r = 0, l=0)),

            # Axis titles
            axis.title.x  = element_text(size=sizes),
            axis.title.y  = element_text(size=sizes, margin = margin(r = 0, l=0), angle=90),
            #axis.title.y = element_blank(),
            #axis.title.x = element_blank(),

            axis.ticks.length.y = unit(0, "cm"),

            panel.grid.minor.x = element_blank(),
            #panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            #panel.grid.major.y = element_blank(),

            legend.position = legend.position.prime,
            legend.text     = element_text(size = sizes),
            legend.title    = element_text(size = sizes + 5),

            panel.spacing = grid::unit(facet_space, "mm"),

            plot.title    = element_text(size=sizes*2),
            plot.subtitle = element_text(size=sizes),

            ...

      ) +

    ## Plot margins

    # If margins is set explicitly, it will apply to each part of the plot
    { if (!is.null(margins))   theme(plot.margin = margins) } +

    # Default is when the "margins" parameter is NULL.
    # In this case the spaces between the genome and transcriptome and "coverage" plots are minimized.
    { if (is.null(margins) & add.coverage)
      theme(plot.margin = unit(c(5,5,1,5), 'mm')) }



    ## Legend labels
    if(!is.null(labels)) { gg.prime <- gg.prime + labels }


    ## PLOT TITLE
    if(!is.na(plot.title)) {
      gg.prime <- gg.prime + labs(title=plot.title)

      if (genome.and.transcripts) {
        ggtr <- ggtr + labs(title=plot.title)
      }
    }

  }

  #### Assemble plot ####
  if(genome.only) {
    plot.exp <- gggenome
  } else if (transcripts.only) {
    plot.exp <- ggtr
  } else if (genome.and.transcripts) {
    plot.exp <- cowplot::plot_grid(ggtr, gggenome, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(genomplot.scale, transcript.plot.scale))
  } else {
          if (add.genome.plot & !add.transcripts.plot) {
      if (!flip.panels) {
        plot.exp <- cowplot::plot_grid(gg.prime, gggenome, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(genomplot.scale, 1))
      } else {
        plot.exp <- cowplot::plot_grid(gggenome, gg.prime, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(1, genomplot.scale))
      }
    } else if (add.genome.plot & add.transcripts.plot) {
      if (!flip.panels) {
        plot.exp <- cowplot::plot_grid(gg.prime, ggtr, gggenome, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(genomplot.scale, transcript.plot.scale, 1))
      } else {
        plot.exp <- cowplot::plot_grid(ggtr, gggenome, gg.prime, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(transcript.plot.scale, 1, genomplot.scale))
      }
    }
    else {
      plot.exp <- gg.prime
    }
  }

  if (return.plot.data) {
    if (genome.only) {
      return(gene.plotsub)
      } else {
      return(plot.sub)     }
    } else {
    return(plot.exp)
    }
}
