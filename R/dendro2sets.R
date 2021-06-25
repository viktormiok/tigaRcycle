###############################################################################
#
#         format converting function - create set names from the dendrogram
#
###############################################################################

dendro2sets <- function(hc){
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