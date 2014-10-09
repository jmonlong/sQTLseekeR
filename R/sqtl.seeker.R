##' \code{sqtl.seeker} is the main function of \code{sQTLseekeR} package. From transcript relative expression,
##' prepared using \code{prepare.trans.exp}, information about the gene location and the path to an ordered genotype
##' file, indexed by function \code{index.genotype}, potential association between each SNP and the transcript
##' relative expression is tested. Potentially, svQTL, i.e. SNPs affecting splicing variability can also be tested
##' to pinpoint potential false sQTL (see Details).
##'
##' A set of filters is automatically used to remove SNPs which are unpractical or not informative. Precisely, these
##' filters remove SNP :
##' \itemize{
##' \item{with more than 3 unknown genotypes}
##' \item{with less than 5 samples in any genotype group}
##' \item{with less than 5 different splicing pattern (needed for permutation efficiency) in any genotype group}}
##'
##' Testing difference in transcript relative expression between genotype groups assumes homogeneity of the variances
##' in these groups. Testing this assumption is more complex and computationnally intensive but if needed the user
##' can choose to test for svQTL, i.e. gene/SNPs where this assumption is violated, by using the \code{svQTL=TRUE}. This
##' test is run in parallel to the sQTL tests, but the copmutation time will be considerably higher. For this reason, another
##' parameter can be tweaked, \code{nb.perm.max.svQTL}, to reduce the number of permutation for the svQTL tests if needed for
##' feasibility reasons.
##'
##' DETAILS ON APPROXIMATION
##'
##' DETAILS ON MD OUTPUT
##' @title sQTL seeker
##' @param tre.df a data.frame with transcript relative expression
##' produced by 'prepare.trans.exp'. 
##' @param genotype.f the name of the genotype file. This file need to
##' be ordered by position, compressed and indexed using tabix (samtools).
##' Must have column 'snpId'.
##' @param gene.loc a data.frame with the genes location. Columns 'chr', 'start',
##' 'end' and 'geneId' are required.
##' @param genic.window the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
##' @param min.nb.ext.scores the minimum number of permuted score higher than
##' 'F.lead' to allow the computation to stop. Default is 1000.
##' @param nb.perm.max the maximum number of permutations. Default is 1e6.
##' @param nb.perm.max.svQTL the maximum number of permutations for the svQTL computation. Default is 1e5.
##' @param svQTL should svQTLs test be performed in addition to sQTLs. Default is FALSE. Warning:
##' computation of svQTLs cannot rely on asymptotic approximation, hence the heavy permutations will
##' considerably increase the running time. 
##' @param approx should the asymptotic distribution be used instead of permutations.
##' Default is TRUE.
##' @return a data.frame with columns
##' \item{snpId}{the SNP name}
##' \item{geneId}{the gene name}
##' \item{F}{the F score}
##' \item{nb.groups}{the number of genotype groups}
##' \item{pv}{the P-value}
##' \item{nb.perms}{the number of permutation used for the P-value computation}
##' @author Jean Monlong
##' @export
sqtl.seeker <- function(tre.df,genotype.f, gene.loc, genic.window=5e3, min.nb.ext.scores=1000,nb.perm.max=1000000,nb.perm.max.svQTL=1e4,svQTL=FALSE,approx=TRUE){

    ## Check if:
    ## - less than 3 missing genotype values
    ## - more than 5 samples per genotype group
    ## - more than 5 different splicing pts per genotype group
    ## Return: TRUE if snp passed, FALSE if not.
    check.genotype <- function(geno.df, tre.df){
        apply(geno.df, 1, function(geno.snp){
            if(sum(as.numeric(geno.snp)==-1)>2){
                return(FALSE)
            } 
            geno.snp.t = table(geno.snp[geno.snp>-1])
            if(sum(geno.snp.t >= 5) < 2){
                return(FALSE)
            }
            nb.diff.pts = sapply(names(geno.snp.t)[geno.snp.t>1], function(geno.i){
                nbDiffPt(tre.df[,which(geno.snp==geno.i)])
            })
            if(sum(nb.diff.pts >= 5) < 2){
                return(FALSE)
            }
            return(TRUE)
        })
    }
    
    analyze.gene.f <- function(tre.gene){
        cat(tre.gene$geneId[1],"\n")
        ## Load genotype
        gr.gene = with(subset(gene.loc, geneId==tre.gene$geneId[1]),
            GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
        gr.gene = GenomicRanges::resize(gr.gene, GenomicRanges::width(gr.gene)+2*genic.window, fix="center")
        if(length(gr.gene)>0){
            genotype.gene = read.bedix(genotype.f, gr.gene)
            
            if(!is.null(genotype.gene)){
                ## Remove samples with non expressed genes
                tre.gene = tre.gene[,!is.na(tre.gene[1,])]
                ## Focus on common samples
                com.samples = intersect(colnames(tre.gene),colnames(genotype.gene))
                ## Filter SNP with not enough power
                snps.to.keep = check.genotype(genotype.gene[,com.samples], tre.gene[,com.samples])
                if(any(snps.to.keep)){
                    genotype.gene = genotype.gene[snps.to.keep, ]
                    
                    tre.dist = hellingerDist(tre.gene[,com.samples])
                    res.df = dplyr::do(dplyr::group_by(genotype.gene, snpId), compFscore(., tre.dist, tre.gene, svQTL=svQTL))
                    res.df = dplyr::do(dplyr::group_by(res.df, nb.groups), compPvalue(., tre.dist, approx=approx, min.nb.ext.scores=min.nb.ext.scores, nb.perm.max=nb.perm.max))
                    if(svQTL){
                        res.df = dplyr::do(dplyr::group_by(res.df, nb.groups), compPvalue(., tre.dist, svQTL=TRUE, min.nb.ext.scores=min.nb.ext.scores, nb.perm.max=nb.perm.max.svQTL))
                    }
                    return(data.frame(done=TRUE,res.df))
                }
            }
        }
        return(data.frame(done=FALSE))
    }
    
    ret.df = dplyr::do(dplyr::group_by(tre.df, geneId), analyze.gene.f(.))
    if(any(ret.df$done)){
        ret.df = subset(ret.df, done)
        ret.df$done=NULL
        return(ret.df)
    } else {
        return(NULL)
    }
}
