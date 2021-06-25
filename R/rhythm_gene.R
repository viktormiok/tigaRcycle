#########################################################################################
#
#       detect rhythmic gene expression over time
#
#########################################################################################

rhythm_gene <- function(infile,
                        filestyle,
                        method,
                        family,
                        timepoints,
                        knot,
                        deg,
                        ncpus){
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
