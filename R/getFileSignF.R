#########################################################################################
#
#         extract the field separator character(FILE_SEP), the set of quoting
#         characters(FILE_QUOTE), the character used for decimal points(FILE_DEC)
#
#########################################################################################
getFileSignF <- function(filestyle){
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
