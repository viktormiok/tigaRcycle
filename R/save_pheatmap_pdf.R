#########################################################################################
#
#       function for saving the pheatmap plot as pdf file
#
#########################################################################################

save_pheatmap_pdf <- function(obj, 
                              filename, 
                              width=10, 
                              height=10) {
                    # obj: pheatmap plot object
                    # filename: location path and file name
                    # width: width of the plot
                    # hight: hight of the plot
                    if (!is(obj, "pheatmap")) {
                      stop("Input (obj) is of wrong class.")
                    }
                    if (!is(filename, "character")) {
                      stop("Input (filename) is of wrong class.")
                    }
                    if (!is(width, "numeric")) {
                      stop("Input (width) is of wrong class.")
                    }
                    if (length(width) != 1) {
                      stop("Input (width) is of wrong length.")
                    }
                    if (is.na(width)) {
                      stop("Input (width) is not a positive integer.")
                    }
                    if (width < 0) {
                      stop("Input (width) is not a positive integer.")
                    }
                    if (!is(height, "numeric")) {
                      stop("Input (height) is of wrong class.")
                    }
                    if (length(height) != 1) {
                      stop("Input (height) is of wrong length.")
                    }
                    if (is.na(height)) {
                      stop("Input (height) is not a positive integer.")
                    }
                    if (height < 0) {
                      stop("Input (height) is not a positive integer.")
                    }
                    stopifnot(!missing(obj))
                    stopifnot(!missing(filename))
                    pdf(filename,
                        width=width,
                        height=height
                    )
                    grid::grid.newpage()
                    grid::grid.draw(obj$gtable)
                    dev.off()
}
