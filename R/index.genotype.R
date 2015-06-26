##' Index genotype file. The file must be tab-separated and ordered by chr, start, end
##' position. Moreover, the first four columns must be 'chr', 'start', 'end' and 'snpId'
##' information. Then each additional column is named using a sample ID and represents
##' its genotype. The genotype is coded as follow: 0 for ref/ref, 1 for ref/mut, 2 for
##' mut/mut and -1 for missing values. The input file should be a text file.
##' @title Index the genotype file
##' @param file the name of the genotype file
##' @return the name of the indexed genotype file.
##' @author Jean Monlong
##' @export
index.genotype <- function(file){
  if(!file.exists(file)){
    stop("File not found: ", file)
  }

  ## Check column names and order
  snp50 = read.table(file, sep="\t", header=TRUE, quote="", comment.char = "", as.is=TRUE, nrows=50)
  if(!all(colnames(snp50)[1:4] == c("chr","start","end","snpId"))){
    stop("Missing column or in incorrect order. The first 4 columns must be 'chr', 'start', 'end' and 'snpId'.")
  }
  ## Check that genotypes are discrete
  if(length(unique(unlist(snp50[,-(1:4)])))>10){
    stop("Discrete genotypes are required but the inputed genotypes look continuous.")
  }
  
  file.final = Rsamtools::bgzip(file, dest=paste0(sub(".gz","",file),".bgz"),overwrite=TRUE)
  Rsamtools::indexTabix(file.final, format="bed")
  return(file.final)
}
