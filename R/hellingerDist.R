##' Computes Hellinger distance for each pair of columns from the input
##' data.
##' @title Hellinger distance computation
##' @param sr a matrix or data.frame with the transcript relative
##' expression for each sample (column).
##' @return a matrix with the pairwise distance between samples.
##' @author Jean Monlong
##' @keywords internal
hellingerDist <- function(sr){
    hellingerDist.p <- function(x1, x2) {
        a <- (sqrt(x1) - sqrt(x2))
        b <- sqrt(sum(a*a))  
        return(b)
    }
    k <- dim(sr)[2]
    interdist <- matrix(data=0,nrow=k,ncol=k)
    for(i in 1:k) {
        for(j in 1:i) {
                interdist[i,j] <- hellingerDist.p(sr[,i],sr[,j])          
            }
    }
    interdist = interdist + t(interdist)
    colnames(interdist) = rownames(interdist) = colnames(sr)
    non.na = !is.na(diag(interdist))
    interdist = interdist[non.na, non.na]
    return(as.dist(interdist))
}
