#########################################################################################
#
#       circular plot of the time points
#
#########################################################################################

circular_plot <- function(x,
                          y,
                          factors,
                          color,
                          leg_name,
                          leg_col){
              # x: The log2 fold-changes of the dots to be plotted
              # y: The adjusted p-values of the dots to be plotted
              # factors: The variable wihic separate the dots on the peaces
              #          coresponding to the each time point
              # color: The color of the each dot in the circular plot
              # leg_name: Names of the different colors for the legend
              # leg_col: Color names for the legend
              if (!is(x, "numeric")) {
                stop("Input (x_data) is of wrong class.")
              }
              if (!is(y, "numeric")) {
                stop("Input (y) is of wrong class.")
              }
              if (!is(y, "numeric")) {
                stop("Input (y) is of wrong class.")
              }
              if (!is(factors, "character")) {
                stop("Input (factors) is of wrong class.")
              }
              if (!is(leg_name, "character")) {
                stop("Input (leg_name) is of wrong class.")
              }
              if (!is(leg_col, "character")) {
                stop("Input (leg_col) is of wrong class.")
              }
              if ((length(x) != length(y)) & 
                  (length(x) != length(factors))) {
                stop("Dimension do not correspond")
              }
              if (length(leg_name) != length(leg_col)) {
                stop("Dimension do not correspond")
              }
              circos.par("track.height" = 0.6) 
              
              # initialise the chart giving factor and x-axis.
              circos.initialize(factors = factors, 
                                x = x
              )
              # build regions. 
              circos.trackPlotRegion(factors = factors,
                                     y = y,
                                     panel.fun = function(x, y) {
                                       circos.axis(h = "top",                   
                                                   labels = TRUE,               
                                                   major.tick = TRUE,           
                                                   labels.cex = 0.5,            
                                                   labels.font = 1,             
                                                   direction = "outside",       
                                                   minor.ticks = 4,             
                                                   major.tick.length = 0.1,     
                                                   lwd = 3                      
                                       )
                                     }
              )
              #  add points
              circos.trackPoints(sectors = factors,
                                 x = x, 
                                 y = y, 
                                 col = color,
                                 pch = 16,
                                 cex = 1.2
              )
              # add legend
              legend("bottomleft", 
                     legend = leg_name, 
                     col = leg_col, 
                     pch = 19, 
                     bty = "n", 
                     pt.cex = 2, 
                     cex = .9, 
                     text.col = "black", 
                     horiz = F , 
                     inset = c(0.4, 0.35)
              )
}
