##' Compute the number of different splicing ratios. Too few splicing ratios is
##' a concern for the permutation process. Indeed if many samples exactly the same
##' splicing ratios, the permutation might create many identical permuted F scores.
##' @title Number of different splicing ratio points
##' @param sr a data.frame with the splicing ratios (transcript x sample).
##' @return the number of different splicing ratios
##' @author Jean Monlong
##' @keywords internal
nbDiffPt <- function(sr){
    length(unique(apply(sr,2,paste,collapse="-")))
}
