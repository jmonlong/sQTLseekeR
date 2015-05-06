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
compPvalue <- function(res.df, tre.dist, min.nb.ext.scores=1e3, nb.perm.max=1e6, svQTL=FALSE, approx=TRUE){
    nb.gp = res.df$nb.groups[1]
    
    estNbPerm <- function(pv,min.nb.ext.scores=1000,nb.perm.max=1e6){
      return(min(ceiling(min.nb.ext.scores / pv  + min.nb.ext.scores / 10),nb.perm.max + min.nb.ext.scores / 10))
    }
    
    ado.null <- function(dist.o,nb.null,nb.gp,svQTL=FALSE,approx=TRUE){
      nb.tot = attr(dist.o,"Size")
      groups.f = factor(sample(1:nb.gp, nb.tot, TRUE))
      adonis.comp(dist.o,groups.f,permutations=nb.null,f.perms=TRUE,svQTL=svQTL,approx=approx)
    }
    
    if(svQTL){
      F = res.df$F.svQTL
    } else {
      F= res.df$F
    }
    
    if(length(F)>1){
      F.f = factor(F)
      F.nb.ext.FP = rep(0,length(F.f))
      pv.lead = 1
      nbP.tot = 0
      while( ( pv.lead * nbP.tot < min.nb.ext.scores ) && ( nbP.tot < nb.perm.max ) ){
        nbP.new = min(c(1e6, estNbPerm(pv.lead,min.nb.ext.scores,nb.perm.max)-nbP.tot))
        if(nbP.new > 0){
          FP = ado.null(tre.dist,nbP.new,nb.gp,svQTL=svQTL,approx=approx)
          FP.c = cut(FP,c(-Inf,levels(F.f),Inf),right=FALSE)
          FP.sum = table(FP.c)
          FP.cs = cumsum(FP.sum)[-length(FP.sum)]
          names(FP.cs) = levels(F.f)
          F.nb.ext.FP = F.nb.ext.FP + FP.cs[F.f] ## Pb ? as.character, luckily not
          nbP.tot = nbP.tot + nbP.new
          pv = 1 - ( F.nb.ext.FP / (nbP.tot+1) )
          pv.lead = min(pv)
        }
      }
    } else {
      F.nb.ext.FP = 0
      pv = 1
      nbP.tot = 0
      while( ( pv * nbP.tot < min.nb.ext.scores ) && ( nbP.tot < nb.perm.max ) ){
        nbP.new = estNbPerm(pv,min.nb.ext.scores,1e6)
        if(nbP.new > 0){
          FP = ado.null(tre.dist,nbP.new,nb.gp,svQTL=svQTL,approx=approx)
          F.nb.ext.FP = F.nb.ext.FP + sum(F <= FP)
          nbP.tot = nbP.tot + nbP.new
          pv = (F.nb.ext.FP + 1) / (nbP.tot + 1)
        }
      }
    }

    if(svQTL){
        res.df$nb.perms.svQTL = nbP.tot
        res.df$pv.svQTL = as.numeric(pv)
    } else {
        res.df$nb.perms = nbP.tot
        res.df$pv = as.numeric(pv)
    }

    return(res.df)
}
