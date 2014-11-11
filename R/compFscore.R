##' Compute the F score, max diff ratio difference.
##' @title F score computation
##' @param geno.df a data.frame of one row with the genotype information for each sample.
##' @param tre.dist a distance object from the transcript relative expression.
##' @param tre.df a data.frame with the transcript relative expression. 
##' @param svQTL should svQTL test be performed in addition to sQTL. Default is FALSE.
##' @return a data.frame with columns:
##' \item{F}{the F score.}
##' \item{nb.groups}{the number of groups created by the genotypes.}
##' \item{md}{the maximum difference in splicing ratios between genotype groups.}
##' \item{tr.first, tr.second}{the two transcripts that change the most.}
##' @author Jean Monlong
##' @keywords internal
compFscore <- function(geno.df, tre.dist, tre.df,svQTL=FALSE){
  if(nrow(geno.df)>1){
    stop(geno.df$snpId[1], " SNP is duplicated in the genotype file.")
  }
  geno.snp = geno.df[,labels(tre.dist)]
  if(any(geno.snp==-1)){
    non.na = geno.snp > -1
    geno.snp = geno.snp[non.na]
    tre.dist = as.dist(as.matrix(tre.dist)[non.na, non.na])
  }
  groups.snp.f = factor(as.numeric(geno.snp))
  F.snp = adonis.comp(tre.dist,groups.snp.f,permutations=2,svQTL=FALSE)
  mdt = md.trans(tre.df, groups.snp.f, labels(tre.dist))
  res.df = data.frame(F=F.snp,
      nb.groups=nlevels(groups.snp.f) ,
      md=mdt$md,
      tr.first=mdt$tr.first,
      tr.second=mdt$tr.second,
      stringsAsFactors=FALSE)
  if(svQTL){
      res.df$F.svQTL = adonis.comp(tre.dist,groups.snp.f,permutations=2,svQTL=TRUE)
  }
  return(res.df)
}
