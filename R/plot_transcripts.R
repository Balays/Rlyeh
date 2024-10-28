#' Visualize 5 or 3-primes on a genome
#'
#' @export

plot_transcripts <-
  function(visfrom, visto, plot.data, gene.plotdata, geom = geom_bar(), geom2 = NULL,
           sum.counts.in.window=T, bin_width = 200, sum.fun='sum', add.all.pos=F, y.thresh=0, minus.strand.down=T,
           genome=NA, add.genome.plot=T, flip.panels=F, genome.only=F, add.cageTSS=F, cagefr.clust=NULL, cagefr.clust.dist=6,
           gene.label=T, add.unstranded=T, add.feature=F, gene.feature.width.thresh=700, flip.gene.y=T, force.gene.y=T,
           force.all.gene.down=F, angle=-60, force.one.lab.per.gene=T, force.all.gene.up.and.down=F,
           prime='prime5', crop.FALSE=F, samples=NA, comb.all.samples=NA,
           sizes=16, gene.sizes, gene.label.col='white', palette=pal_npg()(10), alpha=0.8, plot.title=NA,
           tolower=F, vline=NULL, vline2=NULL, breakseq=5000, margins=unit(c(0,0,0,0), "cm"), #margins=margin(10,1,1,5)
           scales='fixed', genomplot.scale=5.5, ylim=c(0, NA), y.log10=F, y.log2=F, ylims.gene=NULL, y.multip=1.5,
           strip.text =  element_text(angle = 45, size = sizes*1.5, hjust = 0.5), facet_space=1,
           gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene),
           facet_cropF = facet_nested(rows = vars(hpi), cols=NULL, scales ),
           labels=NULL, return.plot.data=F, 
           add.introns = T, add.ymin = F,
           ...) {
    
    require(fuzzyjoin)
    require(tidygenomics)
    require(ggh4x)
    require(hrbrthemes)
    require(gggenes)
    require(ggsci)
    
    #### Genome and annotation plot ####
    
    ### subset the features
    visregion      <- data.frame(seqnames=genome, start=visfrom, end=visto)
    
    if(add.genome.plot) {
      
      gene.plotsub   <- dplyr::select(
        genome_join(gene.plotdata, visregion, by=c('seqnames', 'start', 'end')), 
        -c(seqnames.y, start.y, end.y) )
      colnames(gene.plotsub)[is.element(colnames(gene.plotsub), c('seqnames.x', 'start.x', 'end.x'))] <- c('seqnames', 'start', 'end')
      
      gene.plotsub$strand <- factor(gene.plotsub$strand, levels=c('+', '-', '*'))
      
      ### (Re-)assign ymin
      if (add.ymin) {
        for (i in unique(gene.plotsub$cluster)) {
          gene.plotsub$ymin[gene.plotsub$cluster == i] <- seq(1, nrow(gene.plotsub[gene.plotsub$cluster == i, ]))
        }
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
      
      ## flip
      if (flip.gene.y)  { gene.plotsub$ymin <- gene.plotsub$ymin * -1 }
      
      ## force into one row
      gene.plotsub$ylab <- gene.plotsub$ymin
      if (force.gene.y) { 
        gene.plotsub$ymin[gene.plotsub$strand=='+'] <- 0 
        gene.plotsub$ymin[gene.plotsub$strand=='-'] <- 1 
        gene.plotsub$ymin[gene.plotsub$strand=='*'] <- 0.5
      }
      
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
      
      if (force.all.gene.up.and.down) {
        gene.feature.upwards   <- gene.stranded$gene[gene.stranded$strand == '+']
        gene.feature.downwards <- gene.stranded$gene[gene.stranded$strand == '-']
      }
      
      
      ## add unstranded !!!
      #gene.unstranded$start[gene.unstranded$width < 100] <- gene.unstranded$start - 20
      #gene.unstranded$end  [gene.unstranded$width < 100] <- gene.unstranded$end   + 20
      
      ## force unstranded y position
      try({ gene.unstranded$ymin <- max(gene.plotsub$ymin) })
      
      ## force giving only one label for each gene (in case of multi-exoned genes)
      if (force.one.lab.per.gene) {
        multi.genes <- dup(gene.stranded$gene)
        for (mgene in multi.genes) {
          gene.stranded$gene[gene.stranded$gene == mgene & !is.na(gene.stranded$gene)][1]  <- mgene
          gene.stranded$gene[gene.stranded$gene == mgene & !is.na(gene.stranded$gene)][-1] <- NA
        }
      }
      
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
      if (is.null(ylims.gene)) {
        ylims <- c(min(c(gene.plotsub$ymin, cagefr.clust$ymin))-y.multip,
                   max(c(gene.plotsub$ymin, cagefr.clust$ymax))+y.multip)
      } else {
        ylims <- ylims.gene
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
        { if(gene.label & add.feature ) 
          geom_feature(
            data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
            aes(x = xmean, y = ymin), forward = 1,
            feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
            feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm"))
        } +
        { if(gene.label & add.feature ) 
          geom_feature_label(
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
          #size=gene.sizes$gene_label_height,
          feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
          feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm"))
        } +
        { if(gene.label & force.all.gene.down )
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
                label=gene, angle=angle),
            size=gene.sizes$gene_label_height)
        } +
        
        ## Force (+) gene features pointing upwards and (-) downwards
        # Down
        { if(gene.label & force.all.gene.up.and.down ) geom_feature(
          data = gene.stranded[is.element(gene.stranded$gene, gene.feature.downwards),],
          aes(x = xmean, y = ymin), forward = 1,
          #size=gene.sizes$gene_label_height,
          feature_height = grid::unit(-gene.sizes$gene_feature_height, "mm"),
          feature_width  = grid::unit(-gene.sizes$gene_feature_width,  "mm"))
        } +
        { if(gene.label & force.all.gene.up.and.down )
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
                label=gene, angle=angle),
            size=gene.sizes$gene_label_height)
        } +
        # Up
        { if(gene.label & force.all.gene.up.and.down ) geom_feature(
          data = gene.stranded[is.element(gene.stranded$gene, gene.feature.upwards),],
          aes(x = xmean, y = ymin), forward = 1,
          #size=gene.sizes$gene_label_height,
          feature_height = grid::unit(gene.sizes$gene_feature_height, "mm"),
          feature_width  = grid::unit(gene.sizes$gene_feature_width,  "mm"))
        } +
        { if(gene.label & force.all.gene.up.and.down )
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
                label=gene, angle=angle),
            size=gene.sizes$gene_label_height)
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
        { if(add.unstranded & force.all.gene.down ) geom_feature(
          data = gene.unstranded,
          aes(x = xmean, y = ymin), forward = 1,
          feature_height = grid::unit(-gene.sizes$unstranded_feature_height, "mm"),
          feature_width  = grid::unit(-gene.sizes$unstranded_feature_width,  "mm"))
        } +
        { if(add.unstranded & add.feature ) geom_feature_label(
          data = gene.unstranded,
          aes(x = xmean, y = ymin, label=gene), forward = 1,
          feature_height = grid::unit(gene.sizes$unstranded_feature_label_height, "mm"),
          label_height   = grid::unit(gene.sizes$unstranded_feature_label_text,   "mm") #,size=gene.sizes[5]
        ) 
        } +
        { if(add.unstranded ) geom_text(
          data = gene.unstranded,
          aes(x = xmean, y = ymin-gene.sizes$gene_feature_label_height,
              label=gene, angle=angle),
          size=gene.sizes$gene_label_height)  
        } +
        
        ## Introns
        { if(add.introns)
          geom_segment(data=gene.stranded[gene.stranded$exon_number != 1,], 
                       aes(x=intron_start, xend=intron_end, y=ypos, yend=ypos), 
                       color='black', linetype=gene.sizes$linetype, linewidth=gene.sizes$linewidth) 
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
        
        
        ## x-axis scale
        scale_x_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), 
                           breaks = seq(0, l_genome, breakseq) ) +
        ## y-axis scale
        scale_y_continuous(breaks = unique(gene.plotsub$ymin)) +
        
        ## theme components
        theme_ipsum() +
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
      
      
      ## Filling colors
      strands <- unique(as.character(gene.plotsub$strand)) #; strands <- strands[order(strands)
      if ( all(is.element(strands, c('+', '-'))) )     {
        gggenome <- gggenome + scale_fill_manual(values=alpha(palette[c(1,2)],   alpha=alpha$gene_geom))
      } else 
        if ( all(is.element(strands, c('+', '*'))) )     {
          gggenome <- gggenome + scale_fill_manual(values=alpha(palette[c(3,1)],   alpha=alpha$gene_geom))
        } else 
          if ( all(is.element(strands, c('-', '*')))  )    {
            gggenome <- gggenome + scale_fill_manual(values=alpha(palette[c(2,3)],   alpha=alpha$gene_geom))
          } else
            if ( all(is.element(strands, c('+', '-', '*'))) ){
              gggenome <- gggenome + scale_fill_manual(values=alpha(palette[c(2,3,1)], alpha=alpha$gene_geom))
            }
      
      
    }
    
    
    
    #### Genome coverage plot ####
    if(!genome.only) {
      
      plot.sub <- plot.data[plot.data$start >= visfrom & plot.data$end <= visto , ]
      strands  <- as.character(unique(plot.sub$strand))
      #plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-', '*'))
      
      if(!all(is.na(samples))) {
        plot.sub     <- plot.sub[is.element(plot.sub$hpi, samples), ]
        plot.sub$hpi <- factor(plot.sub$hpi, levels = samples)
      } else { 
        samples <- levels(plot.sub$hpi) 
        
      }
      
      if(!is.na(comb.all.samples)) {
        plot.sub$hpi <- comb.all.samples
      }
      
      if(add.all.pos) {
        fullvisreg <- data.frame(strand=c('+', '-', '*'), start=rep(c(visfrom:visto), 3), end=rep(c(visfrom:visto), 3))
        
        # Use map to create a list of data frames
        fullvisreg.hpi  <- rbindlist(map(unique(plot.sub$hpi),
                                         ~{df_temp <- fullvisreg; df_temp$hpi <- .; df_temp}))
        
        fullvisreg.hpi.correct_prime <- data.frame(rbindlist(map(as.character(unique(plot.sub[,correct_prime])),
                                                                 ~{df_temp <- fullvisreg.hpi; df_temp[,correct_prime] <- .; df_temp})))
        fullvisreg.hpi.correct_prime[,correct_prime] <- as.logical(fullvisreg.hpi.correct_prime[,correct_prime])
        
        plot.sub <- data.frame(merge(data.table(plot.sub),
                                     data.table(fullvisreg.hpi.correct_prime),
                                     by=c('strand', 'start', 'end', 'hpi', correct_prime), all=T) )
        plot.sub[,prime] <- plot.sub[,'start']
        plot.sub$count[is.na(plot.sub$count)] <- 0
        
        ## Factorize
        if (add.unstranded) {
          plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-', '*'))
        } else {
          plot.sub[is.element(plot.sub$strand, c('+', '-')), ]
          plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-'))
        }
        
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
      
      ### Log10
      if(y.log10)   {
        plot.sub$count <- log10(plot.sub$count)
        plot.sub$count[plot.sub$count == '-Inf'] <- 0
      }
      
      if(y.log2)   {
        plot.sub$count <- log2(plot.sub$count)
        plot.sub$count[plot.sub$count == '-Inf'] <- 0
      }
      
      
      ### Minus strand downwards
      if(minus.strand.down)   {
        plot.sub$count[plot.sub$strand == '-'] <- plot.sub$count[plot.sub$strand == '-'] * -1
      }
      
      ## Keep only original strands
      plot.sub <- plot.sub[is.element(plot.sub$strand, strands), ]
      
      ## Factorize
      if (add.unstranded) {
        plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-', '*'))
      } else {
        plot.sub[is.element(plot.sub$strand, c('+', '-')), ]
        plot.sub$strand <- factor(plot.sub$strand, levels=c('+', '-'))
      }
      
      plot.sub$hpi <- factor(plot.sub$hpi, levels=samples)
      
      
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
        { if(!is.null(geom2)) geom2 } +
        scale_fill_manual(values =alpha(palette, alpha$cov_geom)) +
        scale_color_manual(values=alpha(palette, alpha$cov_geom)) +
        
        ## x-axis scale
        scale_x_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE), 
                           breaks = seq(0, l_genome, breakseq) ) +
        
        { if(!is.null(vline))  geom_vline(xintercept = vline,  linetype='dashed', color='blue', size=1) } +
        { if(!is.null(vline2)) geom_vline(xintercept = vline2, linetype='dashed', color='blue', size=1) } +
        
        theme_ipsum() +
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
              panel.spacing=grid::unit(facet_space, "mm"),
              plot.title = element_text(size=sizes*2),
              #legend.position=legend.position.prime,
              ...)
      
      
      ## Legend labels
      if(!is.null(labels)) { gg.prime <- gg.prime + labels }
      
      ## X-axis limits
      gg.prime <- gg.prime + coord_cartesian(xlim = c(visfrom, visto), ylim = ylim)
      
      ## Y-axis limits
      if (!is.null(ylim)) {
        gg.prime <- gg.prime #+
        #coord_cartesian(ylim = ylim)
        #ylim(ylim)
      }
      
      if(y.log10)   {
        gg.prime <- gg.prime + scale_y_log10() } else {
          gg.prime <- gg.prime + scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE))
        }
      
      if(!is.na(plot.title)) {
        gg.prime <- gg.prime + labs(title=plot.title)
      }
      
    }
    
    #### Assemble plot ####
    if(genome.only) {
      plot.exp <- gggenome
    } else {
      if(add.genome.plot) {
        if (!flip.panels) {
          plot.exp <- cowplot::plot_grid(gg.prime, gggenome, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(genomplot.scale, 1))
        } else {
          plot.exp <- cowplot::plot_grid(gggenome, gg.prime, ncol = 1, align = 'v', axis=c('tblr'), rel_heights = c(1, genomplot.scale))
        }
      } else {
        plot.exp <- gg.prime
      }
    }
    
    if (return.plot.data) {
      if (genome.only) {
        return(gene.plotsub) } else {
          return(plot.sub)     }
    }
    else return(plot.exp)
  }
