##' Run sQTLseekeR on a set of genes.
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
##' @param svQTL should svQTL test be performed instead of sQTL. Default is FALSE.
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
sqtl.seeker <- function(tre.df,genotype.f, gene.loc, genic.window=5e3, min.nb.ext.scores=1000,nb.perm.max=1000000,svQTL=FALSE,approx=TRUE){

    ## Check if:
    ## - less than 3 missing genotype values
    ## - more than 5 samples per genotype group
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
        
        if(nrow(genotype.gene)>0){
          ## Remove samples with non expressed genes
          tre.gene = tre.gene[,!is.na(tre.gene[1,])]
          ## Focus on common samples
          com.samples = intersect(colnames(tre.gene),colnames(genotype.gene))
          ## Filter SNP with not enough power
          snps.to.keep = check.genotype(genotype.gene[,com.samples], tre.gene[,com.samples])
          if(any(snps.to.keep)){
            genotype.gene = genotype.gene[snps.to.keep, ]
            
            tre.dist = hellingerDist(tre.gene[,com.samples])
            res.df = dplyr::do(dplyr::group_by(genotype.gene, snpId), compFscore(., tre.dist, tre.gene))
            res.df = dplyr::do(dplyr::group_by(res.df, nb.groups), compPvalue(., tre.dist))
            return(res.df)
          }
        }
      }
    }
    return(data.frame())
    
    dplyr::do(dplyr::group_by(tre.df, geneId), analyze.gene.f(.))
}
