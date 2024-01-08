#Creator: Viktorian Miok
#Day of creation: 05.01.2022
#Description: This file contains functions of interest for the package tigaRcycle

#'association
#'@description This function creates visualization of the a global test results in terms of the contributions of the covariates.
#'@param object A gt.object, usually created by a call to gt function.
#'@param what Gives a choice between various presentations of the same plot.
#'@param cluster Type of the hierarchical clustering performed for the dendrogram.
#'@param alpha Parameter between 0 and 1. Set the level of the family-wise error control in the multiple testing procedure performed on the dendrogram.
#'@param sort If TRUE, the plot sorts the bars with the most significant covariates and subjects to the left, as far as is possible within the constraints of the dendrogram (if present).
#'@param zoom If TRUE, discards non-significant branches from the dendrogram with the corresponding covariates. 
#'@param legend If TRUE, draws a legend in the plot.
#'@param plot If FALSE, suppress all plotting.
#'@param plot_name Title of the plot.
#'@param colors The colors to be used for the bars.
#'@param alias Optional alternative labels for the bars in the plots.
#'@param help.lines If TRUE, prints grey dotted lines that help connect the dendrogram to the bars.
#'@param cex.labels Magnification factor for the x-axis labels.
#'@param ylim Override for the y axis limits of the barplot.
#'@param pdf Optional filename (character) of the pdf file to which the plots are to be written.
#'@param trace If TRUE, prints progress information.
#'@param mirror If TRUE, plots the reverse of the scores for the subjects with negative residual response, so that "good" scores are positive for all subjects.
#'@examples
#'@author Viktorian Miok
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
  #####################################################################################################
  #
  #       visualization of the a global test results in terms of the contributions of the covariates
  #
  ####################################################################################################
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
  ####################################################################################################
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


#'circadian_gene_plot
#'@description This function performs plotting individual gene over time.
#'@param gene_name: Name of the gene of interest, which exist in the row names of the data.
#'@param data1 First data set.
#'@param tp1 Time points of the first data set.
#'@param group1 Group information of the first data set.
#'@param data2 Second data set.
#'@param tp2 Time points of the second data set. 
#'@param group2 Group information of the second data set.
#'@examples
#'@author Viktorian Miok
circadian_gene_plot <- function(gene_name,
                                data1,
                                tp1,
                                group1,
                                data2,
                                tp2,
                                group2){
  #########################################################################################
  #
  #       function for plotting individual gene over time
  #
  ########################################################################################@
  # gene_name: Name of the gene of interest, 
  #            which exist in the row names of the data
  # data1: First data set
  # tp1: Time points of the first data set
  # group1: Group information of the first data set
  # data2: Secon data set
  # tp1: Time points of the first dat set
  # group2: Group information of the second data set
  if (!is(gene_name, "character")) {
    stop("Input (gene_name) is of wrong class.")
  }
  if (!is(data1, "matrix")) {
    stop("Input (data1) is of wrong class.")
  }
  if (!is(tp1, "factor")) {
    stop("Input (tp1) is of wrong class.")
  }
  if (!is(group1, "character")) {
    stop("Input (group1) is of wrong class.")
  }
  if ((ncol(data1) != length(tp1)) & 
      (ncol(data1) != length(group1))) {
    stop("Dimension do not correspond")
  }
  if (!is(data2, "matrix")) {
    stop("Input (data2) is of wrong class.")
  }
  if (!is(tp2, "factor")) {
    stop("Input (tp2) is of wrong class.")
  }
  if (!is(group2, "character")) {
    stop("Input (group2) is of wrong class.")
  }
  if ((ncol(data2) != length(tp2)) & 
      (ncol(data2) != length(group2))) {
    stop("Dimension of the design matrices do not correspond")
  }
  plot_panel <- function(data,
                         tp,
                         diet,
                         gene){
    
    d <- tibble(time=tp, 
                value=as.numeric(data[rownames(data) == gene_name, ]),
                diet=diet
    ) %>% arrange(time)
    
    df <- data_summary(d,
                       varname="value", 
                       groupnames=c("diet", "time")
    )
    
    p1 <- ggplot(df, aes(x=time, 
                         y=mean,
                         group=diet, 
                         color=diet)) + 
      geom_line() + 
      geom_point() + 
      ggtitle("hippocampus") + 
      xlab("time points") + 
      ylab("expression") +
      theme(plot.title=element_text(hjust=0.5)) +
      geom_errorbar(aes(ymin=mean - sd, 
                        ymax=mean + sd), 
                    width=.2, 
                    position=position_dodge(0.05)) + 
      scale_color_manual(values=c("lightskyblue", 
                                  "indianred1")) +
      annotate("rect", 
               xmin=3.5, 
               xmax=Inf, 
               ymin=-Inf,
               ymax=Inf, 
               alpha=.3
      )
  }
  
  p1 <- plot_panel(data=data1,
                   tp=tp1,
                   diet=group1,
                   gene=gene_name
  )
  p2 <- plot_panel(data=data2,
                   tp=tp2,
                   diet=group2,
                   gene=gene_name
  )
  
  grid.arrange(p1,# + ylim(min(hip_data[rownames(hip_data)==gene,])-0.2, max(hip_data[rownames(hip_data)==gene,])+0.2), 
               p2,# + ylim(min(hpt_data[rownames(hpt_data)==gene,])-0.2, max(hpt_data[rownames(hpt_data)==gene,])+0.2),
               nrow=1, 
               top=gene_name
  )
}


#'circadian_heatmap
#'@description This function performs circadian heatmap plot.
#'@param data numeric matrix of the values to be plotted.
#'@param colinf column information.
#'@param sumdat sum data information.
#'@param title main title of the heatmap.
#'@param filename name of the file.
#'@param colors colors used for the image. Defaults to heat colors (heat.colors). 
#'@param distfun function used to compute the distance (dissimilarity) between both rows and columns.
#'@param hclustfun function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms.
#'@param legs information about the legends.
#'@param legc colors for the legends.
#'@examples
#'@author Viktorian Miok
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
  #########################################################################################
  #
  #       circadian heatmap plot 
  #
  #########################################################################################
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
  pdf(file=filename)
  
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
  hm <- heatmap.2(x=mat,
                  col=colors,
                  scale="row",
                  trace="none",
                  main=title,
                  Rowv=F,
                  Colv=F,
                  keysize=1.1,
                  dendrogram="row",
                  distfun=distfun,
                  hclustfun=hclustfun, 
                  ColSideColors=rep(legc, each=6)
  )
  legend("topright",       
         legend=legs,
         col=legc, 
         pch=19 
  )
  dev.off() 
  eval(hm$call)
  legend("topright",       
         legend=legs,
         col=legc, 
         pch=19 
  )
}

#'circular_plot
#'@description This function performs circular plot of the time points.
#'@param x The log2 fold-changes of the dots to be plotted.
#'@param y The adjusted p-values of the dots to be plotted.
#'@param factors The variable which separate the dots on the peaces corresponding to the each time point.
#'@param color The color of the each dot in the circular plot.
#'@param leg_name Names of the different colors for the legend.
#'@param leg_col Color names for the legend.
#'@examples
#'@author Viktorian Miok
circular_plot <- function(x,
                          y,
                          factors,
                          color,
                          leg_name,
                          leg_col){
  #########################################################################################
  #
  #       circular plot of the time points
  #
  #########################################################################################
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
  circos.par("track.height"=0.6) 
  
  # initialise the chart giving factor and x-axis.
  circos.initialize(factors=factors, 
                    x=x
  )
  # build regions. 
  circos.trackPlotRegion(factors=factors,
                         y=y,
                         panel.fun=function(x, y) {
                           circos.axis(h="top",                   
                                       labels=TRUE,               
                                       major.tick=TRUE,           
                                       labels.cex=0.5,            
                                       labels.font=1,             
                                       direction="outside",       
                                       minor.ticks=4,             
                                       major.tick.length=0.1,     
                                       lwd=3                      
                           )
                         }
  )
  #  add points
  circos.trackPoints(sectors=factors,
                     x=x, 
                     y=y, 
                     col=color,
                     pch=16,
                     cex=1.2
  )
  # add legend
  legend("bottomleft", 
         legend=leg_name, 
         col=leg_col, 
         pch=19, 
         bty="n", 
         pt.cex=2, 
         cex=.9, 
         text.col="black", 
         horiz=F , 
         inset=c(0.4, 0.35)
  )
}

#'data_summary
#'@description This function performs data summary per time point for each diet and tissue - function for plotting.
#'@param data Data set, numeric matrix of the values.
#'@param varname Name of the variable.
#'@param groupnames Name of the variable, for grouping.
#'@examples
#'@author Viktorian Miok
data_summary <- function(data,
                         varname,
                         groupnames){
  #########################################################################################
  #
  #       data summary per time point for each diet and tissue - function for plotting 
  #
  #########################################################################################
  # data: Data set
  # varname: Name of the variable
  # groupnames: Name of the variable, for grouping
  if (!is(data, "data.frame")) {
    stop("Input (data) is of wrong class.")
  }
  if (!is(varname, "character")) {
    stop("Input (varname) is of wrong class.")
  }
  if (!is(groupnames, "character")) {
    stop("Input (groupnames) is of wrong class.")
  }
  summary_func <- function(x, col){
    c(mean=mean(x[[col]], 
                na.rm=TRUE),
      sd=sqrt(var(x[[col]],
                  na.rm=TRUE)/length(x[[col]])))/2
  }
  data_sum<-ddply(data,
                  groupnames,
                  .fun=summary_func,
                  varname
  )
  return(data_sum)
}


#'dendro2sets
#'@description format converting function - create set names from the dendrogram.
#'@param hc Data set for conversion.
#'@examples
#'@author Viktorian Miok
dendro2sets <- function(hc){
  ###############################################################################
  #
  #         format converting function - create set names from the dendrogram
  #
  ###############################################################################
  do.sets <- function(x,
                      parent.label, 
                      parent.branch, 
                      struct) {
    for (i in 1:length(x)) {
      parentnum <- as.character(which(struct$names == parent.label))
      currentnum <- as.character(length(struct$names) + 1)
      newset <- labels(x[[i]])
      newname <- paste(parent.label,
                       "[",
                       i,
                       sep=""
      )
      if (is.leaf(x[[i]]))
        newname <- paste(newname, 
                         labels(x[[i]]),
                         sep=":"
        )
      newancestors <- c(struct$ancestors[[parentnum]], 
                        parent.label
      )
      newbranch <- c(parent.branch, i)
      struct$sets[[currentnum]] <- newset
      struct$ancestors[[currentnum]] <- newancestors
      struct$branch[[currentnum]] <- newbranch
      struct$names <- c(struct$names, 
                        newname
      )
      if (!is.leaf(x[[i]])) {
        struct <- do.sets(x=x[[i]], 
                          parent.label=newname, 
                          parent.branch=newbranch, 
                          struct
        )
      }
    }
    return(struct)
  }
  
  struct <- list()
  struct$sets=new.env(hash=TRUE)
  struct$sets$"1" <- labels(hc)
  struct$ancestors <- new.env(hash=TRUE)
  struct$branch <- new.env(hash=TRUE)
  struct$branch$"1" <- numeric(0)
  struct$names <- "O"
  
  if (!is.leaf(hc))
    struct <- do.sets(hc, 
                      "O",
                      numeric(0),
                      struct
    )
  else
    struct$names <- paste(struct$names,
                          labels(hc),
                          sep=":")
  
  #convert to lists
  struct$sets <- as.list(struct$sets)
  names(struct$sets) <-
    struct$names[as.numeric(names(struct$sets))]
  struct$ancestors <- as.list(struct$ancestors)
  names(struct$ancestors) <- 
    struct$names[as.numeric(names(struct$ancestors))]
  struct$branch <- as.list(struct$branch)
  names(struct$branch) <- 
    struct$names[as.numeric(names(struct$branch))]
  
  return(struct)
}

#'dist.pear
#'@description This function performs distance calculation in heatmap.
#'@param x Data set for distance calculation in heatmap
#'@examples
#'@author Viktorian Miok
dist.pear <- function(x) as.dist(1 - cor(t(x)))


#'enrichSelectPath
#'@description This function performs KEGG enrichment analysis of selected pathways.
#'@param genes genes for enrichment analysis.
#'@param selectPathway KEGG names, if NULL all pathway are included.
#'@param cut_off p-value cutoff on enrichment tests to report.
#'@param name name of the comparison.
#'@param orgdb organism data base.
#'@param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'.
#'@examples
#'@author Viktorian Miok
enrichSelectPath <- function(genes, 
                             selectPathway, 
                             cut_off, 
                             name, 
                             orgdb,
                             organism){
  #########################################################################################
  #
  #       function for KEGG enrichment analysis of selected pathways
  #
  #########################################################################################
  # genes: genes for enrichment analysis
  # selectPathway: KEGG names, if NULL all pathway are included
  # cut_off: p-value cutoff on enrichment tests to report
  # name: name of the comparison
  # orgdb: organism data base
  # organism: supported organism listed in 
  #           'http://www.genome.jp/kegg/catalog/org_list.html'
  if (!is(genes, "character")) {
    stop("Input (genes) is of wrong class.")
  }
  if (!is(selectPathway, "character")) {
    stop("Input (selectPathway) is of wrong class.")
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
  if (!is(name, "character")) {
    stop("Input (name) is of wrong class.")
  }
  if (!is(orgdb, "OrgDb")) {
    stop("Input (orgdb) is of wrong class.")
  }
  if (!is(organism, "character")) {
    stop("Input (organism) is of wrong class.")
  }
  sig.gene <- bitr(geneID=genes, 
                   fromType="SYMBOL",
                   toType="ENTREZID",
                   OrgDb=orgdb
  )
  obj <- enrichKEGG(gene=sig.gene[,2],
                    organism=organism,
                    pvalueCutoff=cut_off
  )
  obj <- obj@result
  empt <- data.frame(t(rep('',5)))
  colnames(empt) <- c("Description","p.adjust","geneID","Count","Comparison")
  rownames(empt) <- c("genes")
  
  if(nrow(obj)==0){
    return(empt)
  } else {
    if(nrow(obj[obj$p.adjust < cut_off, ]) > 0){
      for(i in 1:nrow(obj)){
        obj[i,8] <- paste(suppressMessages(bitr(geneID=unlist(strsplit(obj[i, 8], 
                                                                       "/")),
                                                fromType="ENTREZID",
                                                toType="SYMBOL",
                                                OrgDb=orgdb
        )[,2]
        ), collapse="/"
        )
      } 
      if(is.null(selectPathway)){
        res <- obj[obj$p.adjust < cut_off, c(2,6,8,9)]
      } else{
        res <- obj[(obj$ID %in% selectPathway)&(obj$p.adjust < cut_off), c(2,6,8,9)]
      }
      if(nrow(res) > 0){
        res$Comparison <- name
        return(res)
      } else {
        return(empt)
      }
    } else {
      return(empt)
    }  
  }
}



#'getFileSignF
#'@description Extract the field separator character(FILE_SEP), the set of quoting characters(FILE_QUOTE), the character used for decimal points(FILE_DEC)
#'@param filestyle Data set for conversion.
#'@examples
#'@author Viktorian Miok
getFileSignF <- function(filestyle){
  #########################################################################################
  #
  #         extract the field separator character(FILE_SEP), the set of quoting
  #         characters(FILE_QUOTE), the character used for decimal points(FILE_DEC)
  #
  #########################################################################################
  # filestyle: a character vector. 
  #            The data format of input file, must be "txt", or "csv".
  if (!is(filestyle, "character")) {
    stop("Input (filestyle) is of wrong class.")
  }
  if (!is(filestyle, "character")) {
    if (!(filestyle %in% c("csv",
                           "txt"))) {
      stop("Input (filestyle) ill-specified.")
    }
  }
  file_quote2 <- FALSE;
  if (length(filestyle) == 1) {
    if (filestyle == "csv") {
      file_sep <- ",";
      file_quote <- "\"";
      file_quote2 <- TRUE;
      file_dec <- ".";
    } else if (filestyle == "txt") {
      file_sep <- "\t";
      file_quote <- "";
      file_dec <- ".";
    } else {
      stop(c("Please set 'filestyle' before running this function ",
             "(the 'filesyle' could be set as 'txt' or 'csv').\n") );
    }
  } else if (length(filestyle) > 1) {
    if (length(filestyle) == 3) {
      file_sep <- filestyle[1];
      file_quote <- filestyle[2];
      file_quote2 <- TRUE;
      file_dec <- filestyle[3];
    } else {
      stop(c("Please set 'filestyle' before running this function ",
             "(the 'filesyle' should be assigned a vector containing ",
             "three characters, which corresponding to symbols used to ",
             "separate columns, quote values and used for decimal points).\n") );
    }
  } else {
    stop(c("Please set 'filestyle' before running this function ",
           "(the 'filesyle' could be set as 'txt' or 'csv', or could be ",
           "assigned a vector containing three characters, ",
           "which corresponding to symbols used to separate columns, ",
           "quote values and used for decimal points).\n") );
  }
  return(list("sep"=file_sep,
              "quote"=file_quote,
              "quote2"=file_quote2,
              "dec"=file_dec)
  )
}

#'hclust.ave
#'@description This function performs hierahical clustering.
#'@param x Data set for distance calculation in heatmap.
#'@param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#'@examples
#'@author Viktorian Miok
hclust.ave <- function(x) hclust(x, 
                                 method="average"
)



#'pathway_association
#'@description This function performs association of selected pathways with peripheral data.
#'@param x_data The response vector of the regression model.
#'@param y_data The part of the design matrix corresponding to the alternative hypothesis.
#'@param pathways The character vector of pathways names.
#'@param pathways_genes  List object with all the gene names corresponding to the pathway names.
#'@examples
#'@author Viktorian Miok
pathway_association <- function(x_data, 
                                y_data, 
                                pathways, 
                                pathways_genes){
  #########################################################################################
  #
  #       Association of selected pathways with peripheral data
  #
  #########################################################################################
  # x_data: The response vector of the regression model.
  # y_data: The part of the design matrix corresponding to the alternative hypothesis.
  # pathways: The character vector of pathways names
  # pathways_genes:  List object with all the gene names corresponding to the pathway names
  if (!is(x_data, "matrix")) {
    stop("Input (x_data) is of wrong class.")
  }
  if (!is(y_data, "numaric")) {
    stop("Input (y_data) is of wrong class.")
  }
  if (nrow(x_data) != length(y_data)) {
    stop("Dimension of the design matrices do not correspond")
  }
  if (!is(pathways, "charater")) {
    stop("Input (pathways) is of wrong class.")
  }
  if (!is(pathways_genes, "list")) {
    stop("Input (pathways_genes) is of wrong class.")
  }
  if (length(pathways) != length(pathways_genes)) {
    stop("Dimension of the pathways do not correspond")
  }
  pf_tissue_diet <- numeric()
  r_pf_tissue_diet <- data.frame()
  for(i in 1:length(pathways_genes)){
    pf_tissue_diet <- gt(y_data, 
                         x_data[,colnames(x_data)%in%pathways_genes[[i]]]
    )
    if(pf_tissue_diet@result[1] < 0.05){
      pfTissueDiet <- extract(association(object=pf_tissue_diet,
                                          col=c("darkslategray4",
                                                "goldenrod4"
                                          ),
                                          plot=TRUE,
                                          cluster=FALSE,
                                          plot_name=pathways[i]
      )      
      )
      res <- data.frame(pfTissueDiet@result[,1], 
                        pfTissueDiet@extra$direction
      )
      colnames(res) <- c("association", "direction")
      
      res[,2] <- ifelse(pfTissueDiet@extra$direction == 2, 1, -1) 
      res$pathway <- pathways[i]
      r_pf_tissue_diet <- rbind(r_pf_tissue_diet,res)
    }
  }
  return(r_pf_tissue_diet)
}



#'plotUpset
#'@description This function for plotting upset plot.
#'@param analyses Data set.
#'@param cut_off The number of intersections from each set (to cut off at) when aggregating by sets.
#'@param col_interest The columum number to be used.
#'@param text.scale Numeric, value to scale the text sizes, applies to all axis labels, tick labels, and numbers above bar plot.
#'@param matrix.color Color of the intersection points.
#'@param main.bar.color Color of the main bar plot.
#'@param sets.bar.color Color of set size bar plot.
#'@examples
#'@author Viktorian Miok
plotUpset <- function(analyses,  
                      cut_off, 
                      col_interest,  
                      text.scale,
                      matrix.color,
                      main.bar.color,
                      sets.bar.color, ...) {
  #########################################################################################
  #
  #       function for plotting upset plot
  #
  #########################################################################################
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


#'rhythm_gene
#'@description This function for detecting rhythmic gene expression over time.
#'@param infile a character string. The name of input file containing time-series data.
#'@param filestyle a character vector. The data format of input file, must be "txt", or "csv".
#'@param method specify the method for detection of rhythmic signals, must be "tigaR", "JTK", "ARS", "LS".
#'@param family Character string. Either equal to "poisson", "zip" (zero-inflated Poisson), "nb" (negative binomial), or "zinb" (zero-inflated negative binomial): likelihood to be used..
#'@param knot number of knots to used for splines.
#'@param deg degree of splines to be used for splines.
#'@param ncpus Integer. The number of cpus to use for parallel computations.
#'@examples
#'@author Viktorian Miok
rhythm_gene <- function(infile,
                        filestyle,
                        method,
                        family,
                        timepoints,
                        knot,
                        deg,
                        ncpus){
  #########################################################################################
  #
  #       detect rhythmic gene expression over time
  #
  #########################################################################################
  # infile: a character string. The name of input file containing time-series data.
  # filestyle: a character vector. The data format of input file, must be "txt", or "csv"
  # method: specify the method for detection of rhythmic signals,
  #         must be "tigaR", "JTK", "ARS", "LS"
  # family: Character string. Either equal to "poisson", "zip" (zero-inflated Poisson),
  #         "nb" (negative binomial), or "zinb" (zero-inflated negative binomial): 
  #         likelihood to be used.
  # knot: number of knots to used for splines
  # deg: degree of splines to be used for splines
  # ncpus: Integer. The number of cpus to use for parallel computations.
  
  if (!is(infile, "character")) {
    stop("Input (infile) is of wrong class.")
  }
  if (!is(filestyle, "character")) {
    stop("Input (filestyle) is of wrong class.")
  }
  if (!is(filestyle, "character")) {
    if (!(filestyle %in% c("csv",
                           "txt"))) {
      stop("Input (filestyle) ill-specified.")
    }
  }
  if (!is(method, "character")) {
    stop("Input (method) is of wrong class.")
  }
  if (!is(method, "character")) {
    if (!(method %in% c("tigaR",
                        "JTK",
                        "ARS",
                        "LS"))) {
      stop("Input (method) ill-specified.")
    }
  }
  if (!is(family, "character")) {
    stop("Input (family) is of wrong class.")
  }
  if (!is(family, "character")) {
    if (!(family %in% c("nb",
                        "zinb", 
                        "poisson", 
                        "zip",
                        "gaussian"))) {
      stop("Input (family) ill-specified.")
    }
  }
  if (!is(knot, "numeric")) {
    stop("Input (knot) is of wrong class.")
  }
  if (length(knot) != 1) {
    stop("Input (knot) is of wrong length.")
  }
  if (is.na(knot)) {
    stop("Input (knot) is not a positive integer.")
  }
  if (knot < 0) {
    stop("Input (knot) is not a positive integer.")
  }
  if (!is(deg, "numeric")) {
    stop("Input (deg) is of wrong class.")
  }
  if (length(deg) != 1) {
    stop("Input (deg) is of wrong length.")
  }
  if (is.na(deg)) {
    stop("Input (deg) is not a positive integer.")
  }
  if (deg < 0) {
    stop("Input (deg) is not a positive integer.")
  }
  if (!is(ncpus, "numeric")) {
    stop("Input (ncpus) is of wrong class.")
  }
  if (length(ncpus) != 1) {
    stop("Input (ncpus) is of wrong length.")
  }
  if (is.na(ncpus)) {
    stop("Input (ncpus) is not a positive integer.")
  }
  if (ncpus < 0) {
    stop("Input (ncpus) is not a positive integer.")
  }
  ################################################################
  #           read data
  ################################################################
  file_sign <- getFileSignF(filestyle)
  FILE_SEP <- file_sign$sep
  FILE_QUOTE <- file_sign$quote
  FILE_DEC <- file_sign$dec
  
  data <- read.table(infile, 
                     header=TRUE, 
                     sep=FILE_SEP,
                     quote=FILE_QUOTE,
                     dec=FILE_DEC, 
                     stringsAsFactors=FALSE
  )
  ################################################################
  #           tigar method
  ################################################################
  if(method == "tigar"){
    dat <- data.frame(data[,-1],
                      row.names=data[,1]
    )
    p <- vector()
    for(i in colnames(dat)){
      p[i] <- as.numeric(strsplit(i, "[.]")[[1]][2])
    } 
    p[is.na(p)] <- 0
    dat <- as.matrix(dat[,order(p)])
    timefac <- timepoints
    groupfac <- as.factor(sort(p))
    dd <- data.frame(a=groupfac)
    design <- as.matrix(sparse.model.matrix(~ -1 + a, dd))
    
    # define knots used for the p-spline at equal spaced quantiles of the covariate
    knots <- quantile(unique(timefac),
                      seq(0, 1, length=(knot+2))[-c(1, (knot+2))]
    )
    # obtain the design matrix for random coefficients using a radial basis
    z_k <- (abs(outer(timefac, knots, "-")))^deg
    omega <- (abs(outer(knots, knots, "-")))^deg
    svd.omega <- svd(omega)
    sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
    # define variables
    ZSpline <- t(solve(sqrt.omega, t(z_k)))
    
    ################################
    #    shrink hypreparameters
    ################################
    if (family %in% c("zinb", "zip", "nb", "poisson")){
      shrinkA <- tigaRshrinkSeq(form=y ~ 1 + f(timefac,
                                               model="z",
                                               Z=ZSpline,
                                               initial=3,
                                               prior="loggamma",
                                               param=c(1, 0.00001)
      ), 
      dat=dat,
      maxiter=3,
      timefac=timefac,
      groupfac=groupfac,
      ZSpline=ZSpline,
      fams=family,
      shrinkrandom="timefac",
      ncpus=ncpus,
      addpackage=c("splines")
      )
      shrink0 <- tigaRshrinkSeq(form=y ~ 1,
                                dat=dat,
                                maxiter=3, 
                                timefac=timefac,
                                groupfac=groupfac,
                                ZSpline=ZSpline,
                                fams=family,
                                ncpus=ncpus, 
                                addpackage=c("splines")
      )   
    }
    else if (family == "gaussian"){
      shrinkA <- tigaRshrinkSeq(form=y ~ f(timefac,
                                           model="z",
                                           Z=ZSpline,
                                           initial=3,
                                           prior="loggamma",
                                           param=c(1, 0.00001)
      ), 
      dat=dat,
      maxiter=3, 
      timefac=timefac,
      groupfac=groupfac,
      ZSpline=ZSpline,
      shrinkrandom="timefac",
      ncpus=ncpus, 
      addpackage=c("splines")
      )
      
      shrink0 <- tigaRshrinkSeq(form=y ~ 1,
                                dat=dat,
                                maxiter=3, 
                                timefac=timefac,
                                groupfac=groupfac,
                                ZSpline=ZSpline,
                                ncpus=ncpus, 
                                addpackage=c("splines")
      )
    } else{
      stop("Choose correct likelihood to be used.")
    }
    ##############################
    #   fit models
    ##############################
    seqFitA <- tigaRshrinkFit(forms=y ~ f(timefac, 
                                          model="z",
                                          Z =ZSpline,
                                          initial=3,
                                          prior="loggamma", 
                                          param=c(1, 0.00001)
    ),
    dat=dat[1:100,],
    timefac=timefac, 
    groupfac=groupfac,
    ZSpline=ZSpline,
    fams=family,
    shrinksimul=shrinkA, 
    ncpus=ncpus
    )
    seqFit0 <- tigaRshrinkFit(forms=y ~ 1,
                              dat=dat[1:100,],
                              timefac=timefac, 
                              groupfac=groupfac,
                              ZSpline=ZSpline,
                              fams=family,
                              shrinksimul=shrink0, 
                              ncpus=ncpus
    )
    ######################
    #   testing
    ######################
    res <- tdge(dat=dat[1:100,],
                fitAlt=seqFitA,
                fitNull=seqFit0,
                design=design,
                ZSpline=ZSpline
    )
    return(list(result=res, 
                fitA=seqFitA, 
                fit0=seqFit0)
    )
  }
  ################################################################
  #           JTK method
  ################################################################
  else if(method == "JTK"){
    res <- meta2d(infile=infile,
                  filestyle=filestyle, 
                  timepoints=timepoints,
                  cycMethod=method, 
                  outIntegration="noIntegration",
                  outputFile=FALSE
    )$JTK
    return(res)
  }
  ################################################################
  #           ARS method
  ################################################################
  else if(method == "ARS"){
    res <- meta2d(infile=infile,
                  filestyle=filestyle, 
                  timepoints=timepoints,
                  cycMethod=method, 
                  outIntegration="noIntegration",
                  outputFile=FALSE
    )$ARS
    return(res)
  }
  ################################################################
  #           LS method
  ################################################################
  else if(method == "LS"){
    res <- meta2d(infile=infile,
                  filestyle=filestyle, 
                  timepoints=timepoints,
                  cycMethod=method, 
                  outIntegration="noIntegration",
                  outputFile=FALSE
    )$LS
    return(res)
  } else{
    stop("Choose correct method to be used.")
  }
}


#'save_pheatmap_pdf
#'@description This function for saving the pheatmap plot as pdf file.
#'@param obj pheatmap plot object.
#'@param filename location path and file name.
#'@param width width of the plot.
#'@param hight hight of the plot.
#'@examples
#'@author Viktorian Miok
save_pheatmap_pdf <- function(obj, 
                              filename, 
                              width=10, 
                              height=10) {
  #########################################################################################
  #
  #       function for saving the pheatmap plot as pdf file
  #
  #########################################################################################
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


#'shrink_dge
#'@description This function for differential expression analysis at each time point using schrinkage.
#'@param dds DESeqDataSet object. Result of DESeq function from DESeq2 package.
#'@param method Character. Indicated whether DESeq2 or ShrinkBayes package is employed for the analysis.
#'@param contrast Character string. Indicating the which groups to be compared.
#'@param family Character string. Either equal to "poisson", "zip" (zero-inflated Poisson), "nb" (negative binomial), or "zinb" (zero-inflated negative binomial): likelihood to be used.
#'@param num_gene Integer. The number of genes for ShrinkBayes method testing
#'@param maxiter Integer. Number of iteration of ShrinkBayes method.
#'@param ncpus Integer. The number of cpus to use for parallel computations.
#'@examples
#'@author Viktorian Miok
shrink_dge <- function(dds, 
                       method="shrinkbayes",
                       contrast,
                       family="nb",
                       num_gene=100,
                       maxiter=3,
                       ncpus=2){
  #########################################################################################
  #
  #       function for differential expression analysis at each time point using schrinkage
  #
  #########################################################################################
  # dds: DESeqDataSet object. Result of DESeq function from DESeq2 package
  # method: Character. Indicated whether DESeq2 or ShrinkBayes package is 
  #         employed for the analysis
  # contrast: Character string. Indicating the which groups to be compared
  # family: Character string. Either equal to "poisson", 
  #         "zip" (zero-inflated Poisson), "nb" (negative binomial), or
  #         "zinb" (zero-inflated negative binomial): likelihood to be used.
  # num_gene: Integer. The number of genes for shrinkbayes method testing
  # maxiter: Integer. Number of iteration of shrinkbayes meethod
  # ncpus: Integer. The number of cpus to use for parallel computations.
  if (!is(method, "character")) {
    stop("Input (method) is of wrong class.")
  }
  if (!is(method, "character")) {
    stop("Input (method) is of wrong class.")
  }
  if (!is(method, "character")) {
    if (!(method %in% c("shrinkbayes",
                        "deseq2"))) {
      stop("Input (method) ill-specified.")
    }
  }
  if (!is(family, "character")) {
    stop("Input (family) is of wrong class.")
  }
  if (!is(family, "character")) {
    if (!(family %in% c("nb",
                        "zinb", 
                        "poisson", 
                        "zip",
                        "gaussian"))) {
      stop("Input (family) ill-specified.")
    }
  }
  if (!is(num_gene, "numeric")) {
    stop("Input (num_gene) is of wrong class.")
  }
  if (length(num_gene) != 1) {
    stop("Input (num_gene) is of wrong length.")
  }
  if (is.na(num_gene)) {
    stop("Input (num_gene) is not a positive integer.")
  }
  if (num_gene < 0) {
    stop("Input (num_gene) is not a positive integer.")
  }
  if (!is(maxiter, "numeric")) {
    stop("Input (maxiter) is of wrong class.")
  }
  if (length(maxiter) != 1) {
    stop("Input (maxiter) is of wrong length.")
  }
  if (is.na(maxiter)) {
    stop("Input (maxiter) is not a positive integer.")
  }
  if (maxiter < 0) {
    stop("Input (num_gene) is not a positive integer.")
  }
  if (!is(ncpus, "numeric")) {
    stop("Input (ncpus) is of wrong class.")
  }
  if (length(ncpus) != 1) {
    stop("Input (ncpus) is of wrong length.")
  }
  if (is.na(ncpus)) {
    stop("Input (ncpus) is not a positive integer.")
  }
  if (ncpus < 0) {
    stop("Input (ncpus) is not a positive integer.")
  }
  if(method == "shrinkbayes"){
    cond <- dds@colData@listData[[contrast[1]]]
    data <- round(counts(dds,
                         normalized=TRUE),
                  0)[1:num_gene, cond%in%contrast[2:3]]
    groupfac <- factor(as.numeric(factor(cond[cond%in%contrast[2:3]]))
    )
    #Formula as used by INLA
    form=y ~ groupfac
    #Simultaneous shrinkage for 'groupfac' 
    shrink <- ShrinkSeq(form=form, 
                        dat=data,
                        shrinkfixed="groupfac",
                        ncpus=ncpus, 
                        maxiter=maxiter
    )
    #Fit all using the priors resulting from ShrinkSeq 
    fit <- FitAllShrink(forms=form,
                        dat=data,
                        fams=family,
                        shrinksimul=shrink,
                        ncpus=ncpus
    )
    #Find nonparametric prior for differences
    npprior <- NonParaUpdatePrior(fitall=fit,
                                  modus="fixed", 
                                  shrinkpara="groupfac",
                                  ncpus=ncpus, 
                                  maxiter=3
    )
    #Update posteriors using the nonparametric prior
    nppostshr <- NonParaUpdatePosterior(fit,
                                        npprior,
                                        ncpus=ncpus
                                        
    )
    #Compute local fdrs 
    lfdr <- SummaryWrap(nppostshr,
                        thr=0.1
    )
    #Compute Bayesian FDRs
    BFDRs <- BFDR(lfdr)
    out <- data.frame(results(dds, contrast)[1:num_gene,],
                      BFDRs)
    colnames(out) <- c("baseMean",
                       "log2FoldChange",
                       "lfcSE", 
                       "stat", 
                       "pvalue", 
                       "padj",
                       "BFDRs"
    )
    return(out)
    
  } else if(method == "deseq2"){
    res <- results(dds, contrast)
  }else{
    stop("Choose correct method to be used.")
  }
}

