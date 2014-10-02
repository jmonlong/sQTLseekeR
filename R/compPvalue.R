##' Compute P-values from F scores.
##' @title P-values computation
##' @param res.df a data.frame with the F scores and number of groups.
##' @param tre.dist the distance object of the transcript relative expression.
##' @param min.nb.ext.scores the minimum number of permuted score higher than
##' 'F.lead' to allow the computation to stop. Default is 1000.
##' @param nb.perm.max the maximum number of permutations. Default is 1e6.
##' @param svQTL should svQTL test be performed instead of sQTL. Default is FALSE.
##' @param approx should the asymptotic distribution be used instead of permutations.
##' Default is TRUE.
##' @return an updated data.frame with new columns 'pv' and 'nb.perms'.
##' @author Jean Monlong
##' @keywords internal
compPvalue <- function(res.df, tre.dist, min.nb.ext.scores=1e3, nb.perm.max=1e6, svQTL, approx){

    perm.F = compute.null(max(res.df$F),tre.dist,res.df$nb.groups,min.nb.ext.scores,nb.perm.max,svQTL=svQTL,approx=approx)
    res.df$nb.perms = perm.F$nbP.tot

    compute.pv <- function(F,F.perms){
        if(length(F)>1){
            F.f = factor(F)
            FP.c = cut(F.perms,c(-Inf,levels(F.f),Inf),right=FALSE)
            FP.sum = summary(FP.c,maxsum=nlevels(FP.c))
            FP.cs = cumsum(FP.sum)[-length(FP.sum)]
            names(FP.cs) = levels(F.f)
            FP.cs.f = FP.cs[F.f] ## Pb ? as.character, luckily not
            names(FP.cs.f) = names(F)
            pv = 1 - ( FP.cs.f / (length(F.perms)+1) )
        } else {
            pv = (sum(F <= F.perms) + 1) / (length(F.perms) + 1)
            names(pv) = names(F)
        }
        return(pv)
    }
    
    res.df$pv = compute.pv(res.df$F,perm.F$FP)
    return(res.df)
}
