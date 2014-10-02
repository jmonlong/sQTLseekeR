##' Index genotype file. The file must be ordered by chr, start, end position. Moreover,
##' the first three columns must be 'chr', 'start' and 'end' information.
##' @title Index genotype file
##' @param file the name of the genotype file
##' @return the name of the indexed genotype file.
##' @author Jean Monlong
##' @export
indexGenotype <- function(file){
    file.final = Rsamtools::bgzip(file, dest=paste0(sub(".gz","",file),".bgz"),overwrite=TRUE)
    Rsamtools::indexTabix(file.final, format="bed")
    return(file.final)
}
