#########################################################################################
#
#       function for KEGG enrichment analysis of selected pathways
#
#########################################################################################

enrichSelectPath <- function(genes, 
                             selectPathway, 
                             cut_off, 
                             name, 
                             orgdb,
                             organism){
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
