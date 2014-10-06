##' Retrieves sQTLs after multiple test correction and, potentially removing svQTLs.
##' @title sQTLs retrieval
##' @param res.df a data.frame with the results from 'sqtl.seeker' with the P-values
##' for each gene/SNP test.
##' @param FDR the False Discovery Rate to call an association significant.
##' @param md.min the minimum MD (Maximum Difference) in relative expression desired. Maximum difference in relative expression (MD) gives an idea of the effect size of the association.
##' @param out.pdf the name of the pdf file to create. If NULL (default), no pdf output
##' will be created. If non-NULL, the distribution of the P-values and a semi-volcano plot
##' (P-value vs MD) will be shown.
##' @param svQTL.removal if TRUE (and column 'pv.svQTL' is present in 'res.df') significant
##' sQTL which are also significant svQTLs are not reported.
##' @return a subset of the input data.frame with only significant sQTLs.
##' @author Jean Monlong
##' @export
sqtls <- function(res.df, FDR=.01, md.min=.01, out.pdf=NULL, svQTL.removal=TRUE){
    res.df$qv = qvalue::qvalue(res.df$pv)$qvalues

    if(svQTL.removal & any(colnames(res.df)=="pv.svQTL")){
        res.df$qv.svQTL = qvalue::qvalue(res.df$pv.svQTL)$qvalues
        res.df = subset(res.df, qv.svQTL <= FDR)
    }

    if(!is.null(out.pdf)){
        pdf(out.pdf, 8,6)
        print(ggplot::ggplot(res.df, ggplot::aes(x=pv)) +
              ggplot::histogram() + ggplot::theme_bw() +
              ggplot::xlab("P-value") + ggplot::ylab("number of gene/SNP")
              )
    }

    res.df = subset(res.df, qv<=FDR & md>=md.min)

    if(!is.null(out.pdf)){
        if(nrow(res.df)>0){
            print(ggplot::ggplot(res.df, ggplot::aes(y=-log10(pv),x=md)) +
                  ggplot::point() + ggplot::theme_bw() +
                  ggplot::ylab("-log10(P-value)") + ggplot::ylab("Relative expression maximum difference")
                  )
        }
        dev.off()
    }
    return(res.df)
}
