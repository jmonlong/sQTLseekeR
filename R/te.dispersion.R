##' Compute the dispersion from a distance matrix. The dispersion is computed as
##' the average distance to the centroid.
##' @title Dispersion computation
##' @param d a distance matrix with pairwise distance between the samples.
##' @return a value for the dispersion
##' @author Jean Monlong
##' @keywords internal
te.dispersion <-  function(d){
    if(class(d) != "dist"){
        if(any(na.ind <- is.na(diag(d)))){
            d = d[!na.ind,!na.ind]
        }
        d = as.dist(d)
    }
    bd = tryCatch(vegan::betadisper(d,rep(1,ncol(as.matrix(d)))),
        error = function(e)return(list(distances=0)))
    return(mean(bd$distances))
}
