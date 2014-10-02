##' Compute the null distribution for the F scores. 
##' @title F score null distribution
##' @param F.lead the highest F score to drive the permutation optimization.
##' @param dis the distance object of the transcript relative expression.
##' @param nb.gp the number of group to simulate.
##' @param min.nb.ext.scores the minimum number of permuted score higher than
##' 'F.lead' to allow the computation to stop. Default is 1000.
##' @param nb.perm.max the maximum number of permutations. Default is 1e6.
##' @param svQTL should svQTL test be performed instead of sQTL. Default is FALSE.
##' @param approx should the asymptotic distribution be used instead of permutations.
##' Default is TRUE.
##' @return a list with
##' \item{FP}{a vector of the permuted F scores.}
##' \item{nbP.tot}{the number of permuted F scores.}
##' @author Jean Monlong
##' @keywords internal
compute.null <- function(F.lead,dis,nb.gp,min.nb.ext.scores=1e3,nb.perm.max=1e6,svQTL=FALSE,approx=TRUE){
        
    estNbPerm <- function(pv,min.nb.ext.scores=1000,nb.perm.max=1e6){
        return(min(ceiling(min.nb.ext.scores / pv  + min.nb.ext.scores / 10),nb.perm.max + min.nb.ext.scores / 10))
    }

    ado.null <- function(dist.o,nb.null,nb.gp,svQTL=FALSE,approx=TRUE){
        nb.tot = attr(dist.o,"Size")
        groups.f = factor(sample(1:nb.gp, nb.tot, TRUE))
        adonis.comp(dist.o,groups.f,permutations=nb.null,f.perms=TRUE,svQTL=svQTL,approx=approx)
    }
    
    pv.lead = 1
    nbP.tot = 0
    FP = NULL
    while( ( pv.lead * nbP.tot < min.nb.ext.scores ) && ( nbP.tot < nb.perm.max ) ){
        nbP.new = estNbPerm(pv.lead,min.nb.ext.scores,nb.perm.max) - nbP.tot
        if(nbP.new > 0){
            FP = c(FP, ado.null(dis,nbP.new,nb.gp,svQTL=svQTL,approx=approx))
            nbP.tot = nbP.tot + nbP.new
            pv.lead = (sum(FP >= F.lead) + 1) / (nbP.tot + 1)
        }
    }
    return(list(FP=FP,nbP.tot=nbP.tot))
}

