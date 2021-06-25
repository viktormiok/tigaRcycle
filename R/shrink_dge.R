#########################################################################################
#
#       function for differential expression analysis at each time point using schrinkage
#
#########################################################################################

shrink_dge <- function(dds, 
                       method="shrinkbayes",
                       contrast,
                       family="nb",
                       num_gene=100,
                       maxiter=3,
                       ncpus=2){
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
