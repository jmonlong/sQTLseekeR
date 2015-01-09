##' Prepare the transcript expression matrix. Transcript with low
##' expression are removed. Then gene expressing only one transcript
##' are removed because not informative for splicing analysis. Finally
##' the relative expression of the transcripts is retrieved and genes with
##' low dispersion are removed. Also genes with just a few (<25)
##' different splicing patterns are removed as they don't fit the
##' permutation process.
##'
##' After removal of inappropriate transcripts/genes, transcript
##' relative expression is computed by dividing its expression
##' by the total expression of all transcripts for the gene.
##' @title Transcript expression preparation
##' @param te.df a data.frame with the transcript expression. The first
##' two columns are the transcript and gene ids, named 'trId' and
##' 'geneId', then the values for each sample the other columns.
##' @param min.transcript.exp the minimum transcript expression. Transcript
##' with lower expression in all samples are removed.
##' @param min.gene.exp the minimum gene expression. Samples with lower
##' gene expression are removed from the analysis of this gene.
##' @param min.dispersion the minimum dispersion of the transcript
##' relative expression. Genes with lower dispersion are removed.
##' @return a data.frame with the relative transcript expression for the
##' genes to study.
##' @author Jean Monlong
##' @export
prepare.trans.exp <- function(te.df, min.transcript.exp=.01,min.gene.exp=.01, min.dispersion=.1){
    samples = setdiff(colnames(te.df), c("chr","start","end","geneId","trId"))
    samples.sub = sample(samples, min(40,length(samples)))
    trans.to.keep = apply(te.df[,samples],1,function(r)any(r>min.transcript.exp))
    te.df = te.df[trans.to.keep,]
    nb.trans = table(te.df$geneId)
    trans2 = names(which(nb.trans>1))
    te.df = subset(te.df, geneId %in% trans2)

    relativize.filter.dispersion <- function(df){
        df[,samples] = apply(df[,samples], 2,relativize, min.gene.exp=min.gene.exp)
        disp = te.dispersion(hellingerDist(df[,samples.sub]))
        if(disp > min.dispersion & nbDiffPt(df[,samples])>25){
            return(df)
        } else {
            return(data.frame())
        }
    }
    
    te.df = plyr::ldply(lapply(unique(te.df$geneId), function(gene.i){
        df = subset(te.df, geneId==gene.i)
        relativize.filter.dispersion(df)
    }), identity)
    
    return(te.df)
}
