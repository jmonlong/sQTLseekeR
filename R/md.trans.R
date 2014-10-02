##' Compute the maximum difference in transcript relative expression
##' between genotype groups.
##' @title MaxDiff splicing ratios computation
##' @param sr.o a matrix or data.frame with the transcript relative expression.
##' @param groups.o a factor with the genotype groups for each sample.
##' @param samples a vector with the names of samples used in the genotype
##' groups information.
##' @return a list with:
##' \item{md}{the maximum difference in splicing ratios between genotype groups.}
##' \item{tr.first, tr.second}{the two transcripts that change the most.}
##' @author Jean Monlong
##' @keywords internal
md.trans <- function(sr.o, groups.o, samples){
    mTrans = apply(sr.o[,samples], 1, function(sr.r)tapply(sr.r, groups.o, mean, na.rm=TRUE))
    lr = nrow(mTrans)
    ind1 = rep(1:(lr-1),(lr-1):1)
    ind2 = NULL
    for(ii in 2:lr)
        ind2 = c(ind2,ii:lr)
    MDtrans = apply(mTrans,2,function(r)diff(rbind(r[ind1],r[ind2])))
    if(!is.matrix(MDtrans))
        MDtrans = matrix(MDtrans,1)
    gpMD = apply(MDtrans,1,function(e)max(abs(e)))
    gpMD.max = which.max(gpMD)
    tr.first = which.max(abs(MDtrans[gpMD.max,]))
    tr.second = which.max(-sign(MDtrans[gpMD.max,tr.first])*MDtrans[gpMD.max,])
    return(list(md=max(gpMD, na.rm=TRUE), tr.first= sr.o$trId[tr.first], tr.second=sr.o$trId[tr.second]))
}
