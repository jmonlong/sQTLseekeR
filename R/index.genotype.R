##' Index genotype file. The file must be tab-separated and ordered by chr, start, end
##' position. Moreover, the first three columns must be 'chr', 'start' and 'end'
##' information. Then each additional column is named using a sample ID and represents
##' its genotype. The genotype is coded as follow: 0 for ref/ref, 1 for ref/mut, 2 for
##' mut/mut and -1 for missing values. The input file should be a text file.
##' @title Index the genotype file
##' @param file the name of the genotype file
##' @return the name of the indexed genotype file.
##' @author Jean Monlong
##' @export
index.genotype <- function(file){
    file.final = Rsamtools::bgzip(file, dest=paste0(sub(".gz","",file),".bgz"),overwrite=TRUE)
    Rsamtools::indexTabix(file.final, format="bed")
    return(file.final)
}
