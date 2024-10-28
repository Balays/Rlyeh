require(ggplot2)
require(grid)

# Define the custom geom function
geom_gene_rect <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity",
                           na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, 
                           rect_height = grid::unit(3, "mm"), ...) {
  ggplot2::layer(
    stat = stat, data = data, mapping = mapping, geom = GeomGeneRect, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      rect_height = rect_height,
      na.rm = na.rm,
      #linewidth = linewidth,
      ...
    )
  )
}




GeomGeneRect <- ggproto('GeomGeneRect', Geom,
                        required_aes = c('xmin', 'xmax', 'y'),
                        default_aes = aes(
                          colour = "black", # Default outline color
                          fill = "grey",    # Default fill color
                          alpha = 1,        # Default opacity
                          size = 0.5        # Default line width
                        ),
                        
                        draw_panel = function(data, panel_scales, coord, rect_height) {
                          # Apply coordinate transformations if necessary
                          data <- coord$transform(data, panel_scales)
                          n <- nrow(data)
                          
                          grobs <- vector('list', n)
                          for (i in seq_len(n)) {
                            grobs[[i]] <- rectGrob(
                              x = unit(data$xmin[i], "native") + unit(data$xmax[i] - data$xmin[i], "native") / 2,
                              y = unit(data$y[i], "native"),
                              width = unit(data$xmax[i] - data$xmin[i], "native"),
                              height = rect_height,
                              just = "centre",
                              #lwd = linewidth,
                              gp = gpar(
                                fill = alpha(data$fill[i], data$alpha[i]), 
                                col = data$colour[i], 
                                lwd = 0.75, #linewidth, #data$size[i] * .pt,
                                lty = data$linetype[i]
                              )
                            )
                          }
                          
                          class(grobs) <- 'gList'
                          return(grobs)
                        }
)
