
require(igraph)
require(tidygenomics)
require(fuzzyjoin)
require(data.table)

make_read_plot_data <- function(DT, 
                          prime3_cluster_window = 25, cluster.TR.overlaps=T, cluster_neighbours = T, optimal_overlap = T, 
                          is.transcript.pos=T,
                          add.introns=T,
                          genome=NULL, visfrom=NULL, visto=NULL, strand.toplot=NULL,
                          flip_minus=T,
                          add.class_code = ''
                          ) {
  
  plot.data <- DT
  
  ### transcript positions
  if(is.transcript.pos) {
    plot.data[,start := start.TR][,end:=end.TR]
  }
  
  ### Add class_code column if missing
  if(!'class_code' %in% colnames(plot.data)) { plot.data[,class_code := add.class_code]}
  
  
  ### Subset to the region to be visualized
  if (!is.null(visfrom) & !is.null(visto) ) {
    
    visregion    <- data.table(seqnames=genome, start=visfrom, end=visto); setkey(visregion, seqnames, start, end)
    plot.subdata <- foverlaps(plot.data, 
                              visregion,
                              type='any',
                              by.x=c('seqnames', 'start', 'end'), 
                              minoverlap = 1)[!is.na(seqnames) & !is.na(start) & !is.na(end)]
    plot.subdata <- plot.subdata[,start:=i.start][,end:=i.end][,i.start:=NULL][,i.end:=NULL]
    
  } else {
    plot.subdata <- plot.data
  }
  
  
  ### Add introns
  if (add.introns) {
    plot.subdata <- plot.subdata[order(transcript_id, start)]
    
    # Calculate intron start and end
    plot.subdata[, intron_start := c(NA, head(end.exon, -1) + 1), by=transcript_id]
    plot.subdata[, intron_end   := start.exon - 1]
    
    # Filter out non-intron rows (first exon of each transcript will have NA for intron start)
    #intron_dt <- data[!is.na(intron_start), .(TR_ID, intron_start, intron_end)]
  }
  
  
  if (nrow(plot.subdata) > 1) {
  
    ### Cluster only the 3prime ends based on a window
    if(!is.null(prime3_cluster_window)) {
      cluster.prime3 <- plot.subdata[,.(seqnames, start, end, strand, transcript_id)]
      cluster.prime3[strand == '+', prime3:=end]
      cluster.prime3[strand == '-', prime3:=start]
      cluster.prime3[,prime3.win.start:=prime3 - prime3_cluster_window]
      cluster.prime3[,prime3.win.end  :=prime3 + prime3_cluster_window]
      
      cluster.prime3.uni <- 
        setnames(cluster.prime3[,.(seqnames, strand, prime3.win.start, prime3.win.end )], 
                 c('prime3.win.start', 'prime3.win.end'),
                 c('start', 'end'))
      
      cluster.prime3.uni  <- unique(cluster.prime3.uni)
      cluster.prime3.uni  <- data.table(genome_cluster(cluster.prime3.uni, by=c('seqnames', 'start', 'end'), 1, 'cluster'))
      cluster.prime3.uni[,cluster := cluster + 1]
      
      cluster.prime3  <- merge(cluster.prime3, cluster.prime3.uni, 
                               by.x=c('seqnames', 'strand', 'prime3.win.start', 'prime3.win.end'), 
                               by.y=c('seqnames', 'strand', 'start', 'end'))
      
      cluster.prime3  <- unique(cluster.prime3[,.(seqnames, strand, prime3.win.start, prime3.win.end, cluster, transcript_id)])
      
      ## merge clusters with trancsripts
      plot.subdata <- merge(plot.subdata, cluster.prime3,
                         by.x=c('seqnames', 'strand', 'transcript_id'),
                         by.y=c('seqnames', 'strand', 'transcript_id') )
      
      cluster.prime3.stats <- unique(cluster.prime3[,.(cluster.prime3.start = min(prime3.win.start),
                                                       cluster.prime3.end   = max(prime3.win.end)),
                                                    by=.(seqnames, strand, cluster)])
      cluster.prime3.stats[,cluster_width:=abs(cluster.prime3.end - cluster.prime3.start)]
      
    } else {
      ### each prime3 end will be a separate cluster
      plot.subdata[,cluster:=.GRP, by=.(seqnames, prime3)]
    }
    #plot.subdata <- plot.subdata[order(seqnames, cluster, transcript_id, start, end)]
    
    ### cluster based on the overlaps among the transcripts themselves
    if(cluster.TR.overlaps) {
    
      ## cluster starts and ends
      cluster.tr_stats <- unique(plot.subdata[,.(seqnames, strand, start=start.TR, end=end.TR) ])
      cluster.tr_stats[, cluster:=.GRP, by=.(seqnames, strand, start, end)]
      
      ## overlaps between each cluster
      cluster.tr_ov     <- foverlaps2(cluster.tr_stats, cluster.tr_stats, by = c('seqnames', 'strand', 'start', 'end'), minoverlap = prime3_cluster_window)
      cluster.tr_ov_sum <- cluster.tr_ov[, .(cluster, i.cluster)][,ov:=1] 
      
      ## filter out self-overlaps -->> check if neccessary !
      cluster.tr_ov_sum <- cluster.tr_ov_sum[cluster != i.cluster,]
      
      ## include non-overlapping clusters
      cluster_ov_all <- unique(cluster.tr_stats[,.(cluster)])
      cluster_ov_all[,i.cluster:=cluster]
      cluster_ov_all <- merge(cluster.tr_ov_sum, cluster_ov_all, by=c('cluster', 'i.cluster'), all=T)
      
      ## check if the cluster overlaps with its neighbour !
      # Convert cluster_ov_sum table to an adjacency matrix
      cluster.tr_ov_sum <- unique(cluster_ov_all) %>% spread(i.cluster, ov)
      rownames(cluster.tr_ov_sum) <- cluster.tr_ov_sum$cluster
      cluster.tr_ov_sum <- as.matrix(cluster.tr_ov_sum[,-1]) # Exclude the first column which is 'cluster'
      
      # Replace NA with 0
      cluster.tr_ov_sum[is.na(cluster.tr_ov_sum)] <- 0
      
   
      ### Cluster (again) based on overlaps
      if (nrow(cluster.tr_ov) > 1 ) {
      
        if (cluster_neighbours) {
          
          ## cluster starts and ends
          cluster_stats <- plot.subdata[,.(start = min(start), end=(max(end))), 
                                        by=.(seqnames, strand, cluster) ]
          ## overlaps between each cluster
          cluster_ov     <- foverlaps2(cluster_stats, cluster_stats, by = c('seqnames', 'strand', 'start', 'end'), minoverlap = prime3_cluster_window)
          cluster_ov_sum <- cluster_ov[cluster != i.cluster, .(cluster, i.cluster)][,ov:=1] 
          
          ## include non-overlapping clusters
          cluster_ov_all <- unique(cluster_stats[,.(cluster)])
          cluster_ov_all[,i.cluster:=cluster]
          cluster_ov_sum <- merge(cluster_ov_sum, cluster_ov_all, by=c('cluster', 'i.cluster'), all=T)
          
          ## check if the cluster overlaps with its neighbour !
          # Convert cluster_ov_sum table to an adjacency matrix
          cluster_ov_sum <- unique(cluster_ov_sum) %>% spread(i.cluster, ov)
          rownames(cluster_ov_sum) <- cluster_ov_sum$cluster
          cluster_ov_sum <- as.matrix(cluster_ov_sum[,-1]) # Exclude the first column which is 'cluster'
          
          # Replace NA with 0
          cluster_ov_sum[is.na(cluster_ov_sum)] <- 0
          
          if (nrow(cluster_ov_sum) > 1) {
          
            if (!optimal_overlap) {
              # Create a graph from the adjacency matrix
              g <- graph_from_adjacency_matrix(cluster_ov_sum, mode = "undirected", diag = FALSE)
              
              # Find the connected components
              components <- components(g)
              clusters   <- data.frame(cluster = as.numeric(names(components$membership)), cluster2=components$membership)
              
              # Add the new cluster information to original data
              cluster_ov <- merge(cluster_ov, clusters, by='cluster', all=T)
              
            } else {  ### NOT FINISHED THIS PART -->> TOO BIG/COMPLEX GRAPH
              # Create a graph from the adjacency matrix
              g <- graph_from_adjacency_matrix(cluster_ov_sum, mode = "undirected", weighted = TRUE)
              
              # Find the lowest total overlap
              communities <- cluster_optimal(g, weights = E(g)$weight)
              
              # Get the membership vector representing the new cluster assignments
              clusters    <- data.table(cluster=as.numeric(communities$names), cluster2=communities$membership)
              
              # Add the new cluster information to original data
              cluster_ov <- merge(cluster_ov, clusters, by='cluster', all=T)
              
            }
          } else {cluster_ov[,cluster2:= cluster]}
        } else {
          
          ## combining transcript clusters
          cluster_ov_all[,cluster2 := paste(i.cluster, collapse='_'), by=.(cluster)]
          cluster_ov_all[,cluster3 := .GRP, by=.(cluster2)]
          cluster_ov_all[,cluster2 := cluster3]
          
          cluster_ov <- unique(cluster_ov_all[,.(cluster, cluster2)])
          
        }
      
      
        # merge the cluster2 info back to the read data
        plot.subdata <- merge(plot.subdata,
                              unique(cluster_ov[,.(cluster, cluster2)]),
                              by='cluster')
      
      } else {
        plot.subdata[,cluster2:= cluster]
      }
    
  } else {
    plot.subdata[,cluster2:= cluster]
    }
  } else {
    plot.subdata[,cluster := 1]
    plot.subdata[,cluster2:= cluster]
  }
  
  ### Subset to strand
  if (!is.null(strand.toplot)) {
    plot.subdata <- plot.subdata[strand == strand.toplot]  
  } 
    
  
  ### Order the data and add the y-positions
  
  ### Flip negative strand reads to minus ?
  if (flip_minus) {
    
    plot.sub.plus  <- data.table(NULL)
    plot.sub.minus <- data.table(NULL)
    #setorder(plot.sub.plus,  cluster2, prime5)
    #setorder(plot.sub.minus, cluster2, -prime5)
    
    try({ plot.sub.plus  <- plot.subdata[strand == '+', ][order(cluster2, cluster, start, start.exon, end)][, ypos := match(transcript_id, unique(transcript_id)), by = cluster2] })
    try({ plot.sub.minus <- plot.subdata[strand == '-', ][order(cluster2, cluster, start, start.exon, end)][, ypos := match(transcript_id, unique(transcript_id)), by = cluster2] })
    plot.subdata <- data.table(plyr::rbind.fill( plot.sub.plus, plot.sub.minus ))
    
    plot.subdata[strand == '-', ypos:= ypos * -1]
  } else {
    setorder(plot.subdata, cluster2, prime5)
    plot.subdata[, ypos := match(transcript_id, unique(transcript_id)), by = cluster2] 
  }
  
  ### Add 'orientation'
  plot.subdata[, orientation := ifelse(strand == "+", 1, 0)]
  
  ### Add 'xmean'
  plot.subdata$xmean <- apply(plot.subdata[,.(start, end)], 1, mean, na.rm=T)
  
  ### Return the DT
  return(plot.subdata)
}
