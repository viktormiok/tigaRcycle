#########################################################################################
#
#       circadian heatmap plot 
#
#########################################################################################

circadian_heatmap <- function(data,
                              colinf,
                              sumdat,
                              title,
                              filename,
                              colors, 
                              distfun,
                              hclustfun,
                              legs, 
                              legc){
                if (!is(data, "list")) {
                  stop("Input (data) is of wrong class.")
                }
                if (!is(colinf, "character")) {
                  stop("Input (colinf) is of wrong class.")
                }
                if (!is(sumdat, "matrix")) {
                  stop("Input (sumdat) is of wrong class.")
                }  
                if (ncol(sumdat) != length(colinf)) {
                  stop("Dimension of the design matrices do not correspond")
                }
                if (!is(title, "character")) {
                  stop("Input (title) is of wrong class.")
                }
                if (!is(filename, "character")) {
                  stop("Input (filename) is of wrong class.")
                }
                if (!is(colors, "character")) {
                  stop("Input (colors) is of wrong class.")
                }
                if (!is(distfun, "function")) {
                  stop("Input (distfun) is of wrong class.")
                }
                if (!is(hclustfun, "function")) {
                  stop("Input (hclustfun) is of wrong class.")
                }
                if (!is(legs, "character")) {
                  stop("Input (legs) is of wrong class.")
                }
                if (!is(legc, "character")) {
                  stop("Input (legc) is of wrong class.")
                }
                if ((length(data) != length(legs)) & 
                    (length(data) != length(legc))) {
                  stop("Dimension of the design matrices do not correspond")
                }
                pdf(file = filename)
                
                x <- as.matrix(sumdat[rownames(sumdat) %in% setdiff(data[[1]], data[[2]]), colinf %in% names(data[1])])
                Rowv <- rowMeans(x)
                distr <- distfun(x)
                hcr <- hclustfun(distr)
                ddr <- as.dendrogram(hcr)
                ddr <- reorder(ddr, Rowv)
                rowInd <- order.dendrogram(ddr)
                # select all the samples of interest for heatmap
                samp <- colinf %in% c(names(data[1]),names(data[2]))
                mat <- as.matrix(cbind(sumdat[rownames(sumdat) %in% setdiff(data[[1]], data[[2]]), samp][rev(rowInd), ]))
                hm <- heatmap.2(x = mat,
                                col = colors,
                                scale = "row",
                                trace = "none",
                                main = title,
                                Rowv = F,
                                Colv = F,
                                keysize = 1.1,
                                dendrogram = "row",
                                distfun = distfun,
                                hclustfun = hclustfun, 
                                ColSideColors = rep(legc, each = 6)
                )
                legend("topright",       
                       legend = legs,
                       col = legc, 
                       pch = 19 
                )
                dev.off() 
                eval(hm$call)
                legend("topright",       
                       legend = legs,
                       col = legc, 
                       pch = 19 
                )
}
