##' Compute F score using the \code{vegan} package.
##' @title Compute single F score
##' @param dis the distance object of the transcript relative expression.
##' @param groups a factor with the group information.
##' @param permutations the number of permutations.
##' @param f.perms should the permuted F scores be returned instead of the
##' real F score. Default is FALSE.
##' @param svQTL should svQTL test be performed instead of sQTL. Default is FALSE.
##' @param approx should the asymptotic distribution be used instead of permutations.
##' Default is TRUE.
##' @return a vector with the F score or the permuted F scores.
##' @author Jean Monlong
##' @keywords internal
adonis.comp <- function(dis,groups,permutations=99,f.perms=FALSE,svQTL=FALSE,approx=TRUE){
    if(svQTL){
        bd <- vegan::betadisper(dis, groups,type="centroid")
        bd.perm <- permutest.betadisper(bd,control = permute::how(nperm = permutations)) 
        if(f.perms){
            return(bd.perm$f.perms)
        } else {
            return(bd.perm$F)
        }
    } else {
        
        if(f.perms){
            if(approx){
                eigenG <- function (interdist,tol=10^-12) {
                    A <- (- 0.5) * interdist^2
                    n <- ncol(A)
                    I <- diag (1,nrow=n) 
                    J <- matrix(1/n,ncol=n,nrow=n)
                    Aux <- I - J
                    G <- Aux %*% A %*% Aux
                    e <- eigen(G,symmetric=T,only.values=T)$values
                    index <- abs(e) > tol
                    return(e[index])
                }

                approx.dist <- function(dist,nb.mont,nb.gp){
                    dist = as.matrix(dist)
                    e <- eigenG(dist)
                    n = ncol(dist)
                    eigenStats <- c(length (e), sum(e>0), sum(e<0))
                    if (eigenStats[3]>0)  e <- abs(e) 
                    randomChisqN <- matrix(rchisq(nb.mont*eigenStats[1],df=nb.gp-1),
                                           nrow=eigenStats[1],ncol=nb.mont)
                    randomChisqD <- matrix(rchisq(nb.mont*eigenStats[1],df=n-nb.gp),
                                           nrow=eigenStats[1],ncol=nb.mont)
                    asymptNume   <- e %*% randomChisqN 
                    asymptDeno   <- e %*% randomChisqD 
                    asymptF      <- asymptNume / asymptDeno * (n-nb.gp) / (nb.gp-1)
                    return(asymptF)
                }

                return(approx.dist(dis,permutations,nlevels(groups)))
            } else {
                res = vegan::adonis(dis ~ groups,permutations=permutations)
                return(as.numeric(res$f.perms[,1]))
            }
        } else {
            res = vegan::adonis(dis ~ groups,permutations=1)
            return(res$aov.tab[1,4])
        }
    }
}
