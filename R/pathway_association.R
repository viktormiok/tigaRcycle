#########################################################################################
#
#       Association of selected pathways with peripheral data
#
#########################################################################################

pathway_association <- function(x_data, 
                                y_data, 
                                pathways, 
                                pathways_genes){
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
