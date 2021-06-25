#####################################################################################################
#
#       visualization of the a global test results in terms of the contributions of the covariates
#
####################################################################################################

association <- function (object,
                         what=c("p-value", "statistic", "z-score", "weighted"),
                         cluster="average",
                         alpha=0.05, 
                         sort=TRUE, 
                         zoom=FALSE, 
                         legend=TRUE, 
                         plot=TRUE,
                         plot_name="global_test",
                         colors, 
                         alias, 
                         help.lines=FALSE,
                         cex.labels=0.6,
                         ylim, 
                         pdf, 
                         trace,
                         mirror) {
          # object: A gt.object, usually created by a call to gt function
          # what: Gives a choice between various presentations of the same plot
          # cluster: Type of the hierarchical clustering performed for the dendrogram
          # alpha: Parameter between 0 and 1. Set the level of the family-wise error control
          #        in the multiple testing procedure performed on the dendrogram
          # sort: If TRUE, the plot sorts the bars with the most significant covariates and subjects to the left, 
          #       as far as is possible within the constraints of the dendrogram (if present).
          # zoom: If TRUE, discards non-significant branches from the dendrogram with the corresponding covariates. 
          # legend: If TRUE, draws a legend in the plot
          # plot: If FALSE, suppress all plotting.
          # plot_name: Title of the plot
          # colors: The colors to be used for the bars
          # alias: Optional alternative labels for the bars in the plots.
          # help.lines: If TRUE, prints grey dotted lines that help connect the dendrogram to the bars
          # cex.labels: Magnification factor for the x-axis labels.
          # ylim: Override for the y axis limits of the barplot.
          # pdf: Optional filename (character) of the pdf file to which the plots are to be written.
          # trace: If TRUE, prints progress information.
          # mirror:	If TRUE, plots the reverse of the scores for the subjects with negative residual response,
          #         so that "good" scores are positive for all subjects.
          if (!is(object, "gt.object")) {
            stop("Input (object) is of wrong class.")
          }
          if (!is(what, "character")) {
            stop("Input (what) is of wrong class.")
          }
          if (!is(what, "character")) {
            if (!(what %in% c("p-value",
                              "statistic", 
                              "z-score", 
                              "weighted"))) {
              stop("Input (what) ill-specified.")
            }
          } 
          if (!is(alpha, "numeric")) {
            stop("Input (alpha) is of wrong class.")
          }
          if (length(alpha) != 1) {
            stop("Input (alpha) is of wrong length.")
          }
          if (is.na(alpha)) {
            stop("Input (alpha) is not a positive integer.")
          }
          if ((alpha < 0)&(alpha > 1)) {
            stop("Input (cut_off) is not a positive integer.")
          }
          if (!is(sort, "logical")) {
            stop("Input (sort) is of wrong class.")
          }
          if (!is(zoom, "logical")) {
            stop("Input (zoom) is of wrong class.")
          }
          if (!is(legend, "logical")) {
            stop("Input (legend) is of wrong class.")
          }
          if (!is(plot, "logical")) {
            stop("Input (plot) is of wrong class.")
          }
          if (!is(plot_name, "character")) {
            stop("Input (plot_name) is of wrong class.")
          }                
          if (!is(help.lines, "logical")) {
            stop("Input (help.lines) is of wrong class.")
          }
          if (!is(cex.labels, "numeric")) {
            stop("Input (cex.labels) is of wrong class.")
          }
          if (length(cex.labels) != 1) {
            stop("Input (cex.labels) is of wrong length.")
          }
          if (is.na(cex.labels)) {
            stop("Input (cex.labels) is not a positive integer.")
          }
          if (cex.labels < 0) {
            stop("Input (cex.labels) is not a positive integer.")
          }
          if ((length(object) > 1) && missing(pdf)) 
            stop("length(object) > 1. Please reduce to a single test result or specify an output file.")
          if (missing(trace)) 
            trace <- gt.options()$trace
          if (missing(alias)) 
            alias <- NULL
          if (!missing(pdf)) {
            if (tolower(substr(pdf, 
                               nchar(pdf) - 3, 
                               nchar(pdf))) != ".pdf") 
              pdf <- paste(pdf, 
                           ".pdf", 
                           sep=""
              )
            pdf(pdf)
          }
          what <- substr(match.arg(tolower(what), 
                                   c("p-value", 
                                     "statistic", 
                                     "z-score",
                                     "weighted")), 1, 1)
          dendrogram <- plot && (((!is.logical(cluster)) || (cluster)) && (substr(cluster, 1, 1) != "n"))
          if (is.character(legend)) {
            object@legend$cov <- legend
            legend <- TRUE
          }
          if (!is.null(alias)) {
            if (is.environment(alias) || is(alias, "AnnDbBimap")) 
              alias <- as.list(alias)
            if (is.list(alias)) 
              alias <- unlist(alias)
            if (length(alias) == object@functions$df()[3] && is.null(names(alias))) 
              names(alias) <- object@functions$cov.names()
          }
          for (jj in (1:length(object))[size(object) > 0]) {
            obj <- object[jj]
            if (is.null(obj@weights)) 
              weights <- rep(1, size(obj))
            else weights <- obj@weights[[1]]
            if (is.null(obj@subsets)) {
              subset <- seq_len(size(obj))
              obj@subsets <- list(subset)
            }
            else subset <- obj@subsets[[1]]
            ttl <- names(obj)
            if (!is.null(alias(obj))) 
              ttl <- paste(ttl, "-", alias(obj))
            test <- function(set) {
              obj@functions$test(subset[set], weights[set])
            }
            leaves <- matrix(0, size(obj), 5)
            colnames(leaves) <- c("p", "S", "ES", "sdS", "ncov")
            rownames(leaves) <- obj@functions$cov.names(subset)
            progress <- 1
            if (dendrogram) {
              progress <- 1
              K <- 2 * size(obj) - 1
            }
            else {
              K <- size(obj)
              progress <- 0
            }
            digitsK <- trunc(log10(K) + 1)
            for (i in 1:size(obj)) {
              if (trace) {
                if (progress > 0) 
                  cat(rep("\b", digitsK + trunc(log10(progress)) + 
                            4), sep="")
                progress <- progress + 1
                cat(progress, "/", K)
                flush.console()
              }
              leaves[i, ] <- test(i)
            }
            if (is.null(rownames(leaves))) 
              rownames(leaves) <- subset
            if (what == "w") {
              leaves[, c("S", "ES", "sdS")] <- leaves[, c("S", "ES", "sdS")] * matrix(weights(obj),
                                                                                      size(obj),
                                                                                      3
              )
            }
            leaves[leaves[, "sdS"] == 0, "sdS"] <- 1
            pps <- -log10(leaves[, "p"])
            maxlogp <- max(pps[pps != Inf],
                           0, 
                           na.rm=TRUE
            )
            pps[pps == Inf] <- maxlogp + 3
            bars <- switch(what, 
                           p=pps,
                           z=(leaves[, "S"] - leaves[, "ES"])/leaves[, "sdS"],
                           w=, 
                           s=leaves[, "S"]
            )
            names(bars) <- rownames(leaves)
            if (!is.null(alias)) {
              if (length(alias) == length(bars) && (is.null(names(alias)) || 
                                                    is.null(names(bars)))) 
                names(bars) <- alias
              else {
                names(bars) <- alias[names(bars)]
              }
            }
            if (sort) 
              order.bars <- -pps
            else order.bars <- 1:length(bars)
            margins <- par("mai")
            if (is.logical(cluster) && dendrogram) 
              cluster <- "average"
            if (dendrogram) {
              if (dim(leaves)[1] == 1) {
                hc=list(1)
                attr(hc, "label") <- row.names(leaves)[1]
                attr(hc, "members") <- 1
                attr(hc, "leaf") <- TRUE
                attr(hc, "height") <- 0
                attr(hc, "class") <- "dendrogram"
                d2s <- dendro2sets(hc)
                obj@extra <- data.frame(inheritance=p.value(obj))
                if (!is.null(alias)) 
                  alias(obj) <- names(bars)
                obj@subsets <- list(row.names(leaves)[1])
                names(obj) <- d2s$names
                sorter <- unlist(hc)
              }
              else {
                cors <- obj@functions$dist.cor(subset)
                if (obj@directional) 
                  dd <- as.dist(1 - cors)
                else dd <- as.dist(1 - abs(cors))
                hc <- as.dendrogram(hclust(dd, 
                                           method=cluster)
                )
                hc <- reorder(hc, 
                              wts=order.bars,
                              agglo.FUN=min
                )
                sorter <- unlist(hc)
                obj@result=rbind(obj@result,
                                   leaves
                )
                obj@subsets=c(list(obj@functions$cov.names(obj@subsets[[1]])), 
                                unlist(obj@functions$cov.names(subset)))
                obj@extra <- NULL
                obj@structure <- NULL
                if (!is.null(alias)) 
                  alias(obj) <- c("", names(bars))
                d2s <- dendro2sets(hc)
                obj <- inheritance(obj, 
                                   sets=d2s$sets,
                                   ancestors=d2s$ancestors, 
                                   trace=trace, 
                                   Shaffer=TRUE,
                                   homogeneous=TRUE
                )
                obj <- obj[sort.list(names(obj))]
              }
              if (zoom) {
                if (dendrogram) {
                  select <- unique(do.call(c, 
                                           subsets(leafNodes(obj, 
                                                             alpha=alpha))
                                           
                  )
                  )
                  select <- which(rownames(leaves) %in% select)
                }
                else {
                  select <- obj@extra[, "holm"] <= alpha
                }
                sorter <- sorter[sorter %in% select]
              }
              leaves <- leaves[sorter, , drop=FALSE]
              bars <- bars[sorter]
              branch2name <- names(d2s$branch)
              names(branch2name) <- unlist(lapply(d2s$branch, 
                                                  function(br) paste(".", 
                                                                     br,
                                                                     collapse="", 
                                                                     sep="")
              )
              )
              sigcol <- function(tree, branch, sig, top) {
                branch.name <- branch2name[paste(".", 
                                                 branch, 
                                                 collapse="", 
                                                 sep=""
                )]
                sig <- sig && (obj@extra[branch.name, "inheritance"] <= 
                                 alpha)
                uit <- tree
                attr(uit, "edgePar") <- list(col=ifelse(sig, 
                                                        1, 
                                                        gray(0.8)), 
                                                        lwd=ifelse(sig, 2, 1))
                if (sig && top) 
                  attr(uit, "nodePar") <- list(pch=20)
                if (!is.leaf(tree)) {
                  select.branch <- 1:length(tree)
                  if (zoom) {
                    sigs <- lapply(1:length(tree), function(i) {
                      subbranch <- branch2name[paste(".", 
                                                     c(branch, i), 
                                                     collapse="", 
                                                     sep="")]
                      obj@extra[subbranch, "inheritance"] <= 
                        alpha
                    })
                    if (do.call(any, sigs) && (!do.call(all, 
                                                        sigs))) {
                      select.branch <- which(do.call(c, sigs))
                      attrs <- attributes(uit)
                      uit <- list()
                      attributes(uit) <- attrs
                    }
                    attr(uit, "members") <- sum(select %in% unlist(tree))
                  }
                  for (i in 1:length(select.branch)) {
                    uit[[i]] <- Recall(tree[[select.branch[i]]], 
                                       c(branch, select.branch[i]), sig, FALSE)
                  }
                  if (length(uit) == 1) 
                    attr(uit, "midpoint") <- attr(uit[[1]], "midpoint")
                  else {
                    leftmid <- attr(uit[[1]], "midpoint")
                    rightmid <- attr(uit[[2]], "midpoint") + 
                      attr(uit[[1]], "members")
                    attr(uit, "midpoint") <- (leftmid + rightmid)/2
                  }
                }
                else {
                  attr(uit, "midpoint") <- 0
                }
                return(uit)
              }
              hc <- sigcol(hc,
                           numeric(0), 
                           sig=TRUE,
                           top=TRUE
              )
              par(mai=c(0, max(1, margins[2]), margins[3:4]))
              layout(as.matrix(1:2), heights=c(1, 2))
              ylab <- ifelse(obj@directional, "correlation", "abs. correlation")
              plot(hc,
                   leaflab="none", 
                   yaxt="n",
                   ylab=ylab, 
                   mgp=c(4, 1, 0),
                   main=plot_name)
              axis(2,
                   at=seq(0, 2, by=0.2), 
                   labels=1 - seq(0, 2, by=0.2), las=2)
            }
            else {
              obj@subsets <- as.list(row.names(leaves))
              obj@result <- leaves
              obj@extra <- NULL
              if (!is.null(alias)) 
                alias(obj) <- names(bars)
              obj@weights <- NULL
              colnames(obj@result) <- c("p-value", 
                                        "Statistic", 
                                        "Expected", 
                                        "Std.dev",
                                        "#Cov")
              names(obj) <- row.names(leaves)
              sorter <- sort.list(order.bars)
              if (zoom) 
                obj <- p.adjust(obj, "holm")
              obj <- obj[sorter]
              bars <- bars[sorter]
              if (trace) 
                cat("\n")
            }
            labwidth <- max(strwidth(names(bars), "inches", cex.labels)) + 
              0.2
            par(mai=c(max(margins[1], labwidth * 1.3),
                        max(1, margins[2]), 
                        if (dendrogram) 0 else margins[3], margins[4])
            )
            positive <- obj@functions$positive(subset)[sorter]
            if (missing(colors)) 
              if (max(positive) <= 2) 
                colors <- 3:2
            else colors <- rainbow(length(obj@legend$cov), start=0, 
                                   end=1/2)
            if (all(positive %in% 0:1)) 
              cols <- ifelse(positive, colors[1], colors[2])
            else cols <- colors[positive]
            ylab <- switch(what,
                           p="p-value",
                           z="z-score", 
                           s="test statistic", 
                           w="weighted test statistic"
            )
            if (!missing(ylim)) {
              ylims <- ylim
            }
            else {
              ylims <- switch(what, z=, p=range(bars), w=, 
                              s=range(c(bars, leaves[, "ES"])))
              if (ylims[1] > 0) 
                ylims[1] <- 0
            }
            if (legend) {
              nbars <- length(bars)
              room <- (ylims[2] - max(bars[trunc(nbars * 0.6):nbars]))/diff(ylims)
              ylims[2] <- ylims[2] + diff(ylims) * max(0, 0.1 * 
                                                         length(colors) - room)
            }
            if(cluster == "average"){
              mids <- drop(barplot(bars, 
                                   yaxt="n",
                                   border=NA,
                                   las=2, 
                                   ylab=ylab, 
                                   ylim=ylims,
                                   mgp=c(4, 1, 0), 
                                   col=cols, 
                                   cex.names=cex.labels,
                                   plot=plot)
              )
            } else{
              mids <- drop(barplot(bars, 
                                   yaxt="n",
                                   border=NA,
                                   las=2, 
                                   ylab=ylab, 
                                   ylim=ylims,
                                   mgp=c(4, 1, 0), 
                                   col=cols, 
                                   cex.names=cex.labels,
                                   plot=plot,
                                   main=plot_name)
              )
            }
            
            if (dendrogram && help.lines) {
              mb <- max(bars)
              for (i in 1:length(bars)) lines(c(mids[i], mids[i]), 
                                              c(max(0, bars[i]) + 0.01 * mb, mb * 1.2), 
                                              col=gray(0.8), 
                                              lty=3)
            }
            if (plot) {
              if (what == "p") {
                labs <- seq(0, ylims[2], by=max(1, ylims[2]%/%5))
                if (length(labs) == 1) {
                  minp <- 10^-ylims[2]
                  if (minp < 0.5) {
                    labs <- log10(c(1, 2, 10/3, 5))
                    labs <- labs[labs < ylims[2]]
                  }
                  else {
                    fc <- 10^-floor(log10(1 - minp))
                    labs <- -log10(c(1, ceiling(minp * fc)/fc))
                  }
                }
                else if (length(labs) <= 2) 
                  labs <- outer(log10(c(1, 2, 5)), labs, "+")
                else if (length(labs) <= 4) 
                  labs <- outer(log10(c(1, 10/3)), labs, "+")
                if (max(bars) > ylims[2]) 
                  axis(2, at=c(labs, max(bars)), labels=c(10^-labs, 
                                                              0), las=2)
                else axis(2, at=labs, labels=10^-labs, las=2)
              }
              else axis(2, las=2)
              if (what %in% c("s", "w")) {
                sapply(seq_along(mids), function(i) {
                  lines(c(mids[i] - 0.5, mids[i] + 0.5), rep(leaves[i, 
                                                                    "ES"], 2), lwd=3)
                  sapply(seq_len(max(0, (bars[i] - leaves[i,"ES"])/leaves[i, "sdS"])), 
                         function(k) lines(c(mids[i] - 0.5, mids[i] + 0.5), 
                                           rep(leaves[i, "ES"] + k * leaves[i, "sdS"], 2)
                         )
                  )
                })
              }
              abline(0, 0)
              if (legend) 
                legend("topright", obj@legend$cov, fill=colors, 
                       bg="white")
              layout(1)
              par(mai=margins)
              if (!missing(pdf)) {
                title(ttl)
              }
            }
          }
          if (!missing(pdf)) 
            dev.off()
          if (length(object) == 1) {
            out <- obj
          }
          else out <- NULL
          return(invisible(out))
}
