##' \code{sqtl.seeker} is the main function of \code{sQTLseekeR} package. From transcript relative expression,
##' prepared using \code{prepare.trans.exp}, information about the gene location and the path to an ordered genotype
##' file, indexed by function \code{index.genotype}, association between each SNP and the transcript relative
##' expression is tested. Eventually, svQTL, i.e. SNPs affecting splicing variability can also be tested to pinpoint
##' potential false sQTL (see Details).
##'
##' A set of filters is automatically used to remove SNPs which are unpractical or not informative. Precisely, these
##' filters remove SNP with :
##' \itemize{
##' \item{more than 3 unknown genotypes}
##' \item{less than 5 samples in any genotype group}
##' \item{less than 5 different splicing pattern (needed for permutation efficiency) in any genotype group}}
##'
##' Testing difference in transcript relative expression between genotype groups assumes homogeneity of the variances
##' in these groups. Testing this assumption is more complex and computationnally intensive but if needed the user
##' can choose to test for svQTL (splicing variability QTL), i.e. gene/SNPs where this assumption is violated,
##' by using the \code{svQTL=TRUE}. This test is run in parallel to the sQTL tests, but the computation time
##' will be considerably higher. For this reason, another parameter can be tweaked, \code{nb.perm.max.svQTL},
##' to reduce the number of permutation for the svQTL tests if needed for feasibility reasons.
##'
##' The permutation process is optimized by computing one permuted distribution per gene and using a number
##' of permutation depending on how extreme the true scores are compared to the permuted ones. To decrease
##' even more the computation time, an approximation of the null F distribution was given by Anderson &
##' Robinson (2003), as a misture of Chi-square distributions whoose parameters are derived from the eigen
##' values of the distance matrix.
##'
##' In addition to the F score and P-value, the maximum difference(MD) in relative expression between genotype
##' groups is reported. This is to be used as a measure of the size of the effect. For example, if 'md' is 0.2
##' there is one transcript whose relative expression shifted by 20% between two genotype groups.
##' @title sQTL seeker
##' @param tre.df a data.frame with transcript relative expression
##' produced by 'prepare.trans.exp'.
##' @param genotype.f the name of the genotype file. This file need to
##' be ordered by position, compressed and indexed using 'index.genotype' or externally using tabix (samtools).
##' Must have column 'snpId'.
##' @param gene.loc a data.frame with the genes location. Columns 'chr', 'start',
##' 'end' and 'geneId' are required.
##' @param genic.window the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
##' @param min.nb.ext.scores the minimum number of permuted score higher than
##' the highest true score to allow the computation to stop. Default is 1000.
##' @param nb.perm.max the maximum number of permutations. Default is 1e6.
##' @param nb.perm.max.svQTL the maximum number of permutations for the svQTL computation. Default is 1e4.
##' @param svQTL should svQTLs test be performed in addition to sQTLs. Default is FALSE. Warning:
##' computation of svQTLs cannot rely on asymptotic approximation, hence the heavy permutations will
##' considerably increase the running time.
##' @param approx should the asymptotic distribution be used instead of permutations.
##' Default is TRUE.
##' @param verbose Should the gene IDs be outputed when analyzed. Default is TRUE. Mainly for debugging.
##' @return a data.frame with columns
##' \item{geneId}{the gene name}
##' \item{snpId}{the SNP name}
##' \item{F}{the F score}
##' \item{nb.groups}{the number of genotype groups (2 or 3)}
##' \item{md}{the maximum difference in relative expression between genotype groups (see Details)}
##' \item{tr.first/tr.second}{the transcript IDs of the two transcripts that change the most (and symetrically).}
##' \item{pv}{the P-value}
##' \item{nb.perms}{the number of permutation used for the P-value computation}
##' \item{F.svQTL/pv.svQTL/nb.perms.svQTL}{idem for svQTLs, if 'svQTL=TRUE'.}
##' @author Jean Monlong
##' @export
sqtl.seeker <- function(tre.df,genotype.f, gene.loc, genic.window=5e3, min.nb.ext.scores=1000,nb.perm.max=1000000,nb.perm.max.svQTL=1e4,svQTL=FALSE,approx=TRUE, verbose=TRUE){

  . = nb.groups = snpId = NULL ## Uglily appease R checks (dplyr)

  ## Check if:
  ## - less than 3 missing genotype values
  ## - more than 5 samples per genotype group
  ## - more than 5 different splicing pts per genotype group
  ## Return: TRUE if snp passed, FALSE if not.
  check.genotype <- function(geno.df, tre.df){
    apply(geno.df, 1, function(geno.snp){
      if(sum(as.numeric(geno.snp)==-1)>2){
        return("Missing genotype")
      }
      geno.snp.t = table(geno.snp[geno.snp>-1])
      if(sum(geno.snp.t >= 5) < 2){
        return("One group of >5 samples")
      }
      nb.diff.pts = sapply(names(geno.snp.t)[geno.snp.t>1], function(geno.i){
        nbDiffPt(tre.df[,which(geno.snp==geno.i)])
      })
      if(sum(nb.diff.pts >= 5) < 2){
        return("One group of >5 different splicing")
      }
      return("PASS")
    })
  }

  analyze.gene.f <- function(tre.gene){
    if(verbose) message(tre.gene$geneId[1])
    ## Load genotype
    gr.gene = with(gene.loc[which(gene.loc$geneId==tre.gene$geneId[1]),], GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    if(genic.window>0){
      gr.gene = GenomicRanges::resize(gr.gene, GenomicRanges::width(gr.gene)+2*genic.window, fix="center")
    }
    if(length(gr.gene)>0){
      ## Remove samples with non expressed genes
      tre.gene = tre.gene[,!is.na(tre.gene[1,])]
      ## Focus on common samples
      genotype.headers = as.character(read.table(genotype.f, as.is=TRUE, nrows=1))
      com.samples = intersect(colnames(tre.gene),genotype.headers)
      tre.dist = hellingerDist(tre.gene[,com.samples])

      res.df = data.frame()
      gr.gene.spl = gr.gene
      if(any(GenomicRanges::width(gr.gene)>2e4)){
        gr.gene.spl = gr.gene[which(GenomicRanges::width(gr.gene)<=2e4)]
        for(ii in unique(which(GenomicRanges::width(gr.gene)>2e4))){
          pos.breaks = unique(round(seq(GenomicRanges::start(gr.gene[ii]), GenomicRanges::end(gr.gene[ii]),length.out=floor(GenomicRanges::width(gr.gene[ii])/1e4)+1)))
          gr.gene.spl.ii = rep(gr.gene[ii], length(pos.breaks)-1)
          GenomicRanges::start(gr.gene.spl.ii) = pos.breaks[-length(pos.breaks)]
          pos.breaks[length(pos.breaks)] = pos.breaks[length(pos.breaks)]+1
          GenomicRanges::end(gr.gene.spl.ii) = pos.breaks[-1]-1
          gr.gene.spl = c(gr.gene.spl, gr.gene.spl.ii)
        }
      }
      
      res.df = lapply(1:length(gr.gene.spl), function(ii){
        if(verbose){message("  Sub-range ",ii)}
        genotype.gene = read.bedix(genotype.f, gr.gene.spl[ii])
        if(verbose & is.null(genotype.gene)){message("    No SNPs in the genomic range.")}
        if(!is.null(genotype.gene)){
          ## Filter SNP with not enough power
          snps.to.keep = check.genotype(genotype.gene[,com.samples], tre.gene[,com.samples])
          if(verbose){
            snps.to.keep.t = table(snps.to.keep)
            message("    ",paste(names(snps.to.keep.t), snps.to.keep.t,sep=":",collapse=", "))
          }
          if(any(snps.to.keep=="PASS")){
            genotype.gene = genotype.gene[snps.to.keep=="PASS", ]
            res.df = dplyr::do(dplyr::group_by(genotype.gene, snpId), compFscore(., tre.dist, tre.gene, svQTL=svQTL))
            ## res.df = lapply(unique(genotype.gene$snpId), function(snpId){
            ##   data.frame(snpId=snpId, compFscore(genotype.gene[which(genotype.gene$snpId==snpId),], tre.dist, tre.gene, svQTL=svQTL))
            ## })
            ## res.df = plyr::ldply(res.df)
          }
        }
        return(res.df)
      })
      res.df = do.call(rbind, res.df)
      
      if(nrow(res.df)>0){
        res.df = dplyr::do(dplyr::group_by(res.df, nb.groups), compPvalue(., tre.dist, approx=approx, min.nb.ext.scores=min.nb.ext.scores, nb.perm.max=nb.perm.max))
        if(svQTL){
          res.df = dplyr::do(dplyr::group_by(res.df, nb.groups), compPvalue(., tre.dist, svQTL=TRUE, min.nb.ext.scores=min.nb.ext.scores, nb.perm.max=nb.perm.max.svQTL))
        }
        ## res.df = lapply(unique(res.df$nb.groups), function(nbgp.i){
        ##   res.f = res.df[which(res.df$nb.groups==nbgp.i),]
        ##   res.f = compPvalue(res.f, tre.dist, approx=approx, min.nb.ext.scores=min.nb.ext.scores, nb.perm.max=nb.perm.max)
        ##   if(svQTL){
        ##     res.f = compPvalue(res.f, tre.dist, svQTL=TRUE, min.nb.ext.scores=min.nb.ext.scores, nb.perm.max=nb.perm.max.svQTL)
        ##   }
        ##   res.f
        ## })
        ## res.df = plyr::ldply(res.df)
        return(data.frame(done=TRUE,res.df))
      }
    } else {
      if(verbose){
        warning("Issue with the gene location.")
      }
    }
    return(data.frame(done=FALSE))
  }

  ret.df = plyr::ldply(lapply(unique(tre.df$geneId), function(gene.i){
    df = tre.df[which(tre.df$geneId==gene.i),]
    data.frame(geneId=gene.i, analyze.gene.f(df))
  }))

  if(any(ret.df$done)){
    ret.df = ret.df[which(ret.df$done),]
    ret.df$done=NULL
    return(ret.df)
  } else {
    return(NULL)
  }
}
