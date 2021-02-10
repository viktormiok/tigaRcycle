#########################################################################################
#
#       data summary per time point for each diet and tissue - function for plotting 
#
#########################################################################################

data_summary <- function(data,
                         varname,
                         groupnames){
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
                  c(mean = mean(x[[col]], 
                                na.rm = TRUE),
                    sd = sqrt(var(x[[col]],
                                  na.rm = TRUE)/length(x[[col]])))/2
                }
                data_sum<-ddply(data,
                                groupnames,
                                .fun = summary_func,
                                varname
                )
                return(data_sum)
}

