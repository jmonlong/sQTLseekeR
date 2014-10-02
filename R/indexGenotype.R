##' Index genotype file. The file must be ordered by chr, start, end position. Moreover,
##' the first three columns must be 'chr', 'start' and 'end' information.
##' @title Index genotype file
##' @param file the name of the genotype file
##' @return the name of the indexed genotype file.
##' @author Jean Monlong
indexGenotype <- function(file){
    Rsamtools::bgzip(file, overwrite=TRUE)
    Rsamtools::indexTabix(paste0(file,".bgz"), format="bed")
    return(paste0(file,".bgz"))
}
