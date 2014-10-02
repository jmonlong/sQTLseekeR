##' Compute the relative expression of transcripts.
##' @title Relative expression computation
##' @param x a vector with of transcript expression
##' @param min.gene.exp the minimum gene expression in the sample. If the
##' gene expression if too low it is usually safer to remove from the
##' analysis.
##' @return a vector with the relative expression.
##' @author Jean Monlong
##' @keywords internal
relativize <- function(x, min.gene.exp=.01){
    x = as.numeric(x)
    if (!any(is.na(x)) && sum(x) > min.gene.exp) {
        x/sum(x)
    } else {
        return(rep(NA,length(x)))
    }
}
