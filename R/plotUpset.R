#########################################################################################
#
#       function for plotting upset plot
#
#########################################################################################

plotUpset <- function(analyses,  
                      cut_off, 
                      col_interest,  
                      text.scale,
                      matrix.color,
                      main.bar.color,
                      sets.bar.color, ...) {
          # analyses: Data set
          # cut_off: The number of intersections from each set (to cut off at) when aggregating by sets
          # col_interest: The columum number to be used
          # text.scale: Numeric, value to scale the text sizes, applies to all axis labels,
          #             tick labels, and numbers above bar plot. 
          # matrix.color: Color of the intersection points
          # main.bar.color: Color of the main bar plot
          # sets.bar.color: Color of set size bar plot
          if (!is(analyses, "list")) {
            stop("Input (analyses) is of wrong class.")
          }
          if (!is(cut_off, "numeric")) {
            stop("Input (cut_off) is of wrong class.")
          }
          if (length(cut_off) != 1) {
            stop("Input (cut_off) is of wrong length.")
          }
          if (is.na(cut_off)) {
            stop("Input (cut_off) is not a positive integer.")
          }
          if ((cut_off < 0)&(cut_off > 1)) {
            stop("Input (cut_off) is not a positive integer.")
          }
          if (!is(col_interest, "numeric")) {
            stop("Input (col_interest) is of wrong class.")
          }
          if (length(col_interest) != 1) {
            stop("Input (col_interest) is of wrong length.")
          }
          if (is.na(col_interest)) {
            stop("Input (col_interest) is not a positive integer.")
          }
          if (col_interest < 0) {
            stop("Input (cut_off) is not a positive integer.")
          }
          if (!is(text.scale, "numeric")) {
            stop("Input (text.scale) is of wrong class.")
          }
          if (length(text.scale) != 1) {
            stop("Input (text.scale) is of wrong length.")
          }
          if (is.na(text.scale)) {
            stop("Input (text.scale) is not a positive integer.")
          }
          if (text.scale < 0) {
            stop("Input (text.scale) is not a positive integer.")
          }
          if (!is(matrix.color, "character")) {
            stop("Input (matrix.color) is of wrong class.")
          }
          if (!is(main.bar.color, "character")) {
            stop("Input (main.bar.color) is of wrong class.")
          }
          if (!is(sets.bar.color, "character")) {
            stop("Input (sets.bar.color) is of wrong class.")
          }
          sig <- lapply(analyses, 
                        function(x) 
                          rownames(x)[which(x[,col_interest] < cut_off)]
          )
          
          if (sum(sapply(sig, length) > 0) > 1) {
            aux <- UpSetR::fromList(sig)
            UpSetR::upset(aux, 
                          nsets=length(aux),
                          order="freq",
                          text.scale=text.scale,
                          matrix.color=matrix.color,
                          main.bar.color=main.bar.color,
                          sets.bar.color=sets.bar.color
            )
          } else {
            stop("Two non-empty sets are needed for an upset plot. There is ",
                 sum(sapply(sig, length) > 0), 
                 "."
            )
          }
}
