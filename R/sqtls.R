##' Retrieves sQTLs after multiple testing correction and potentially removing svQTLs.
##' The distribution of the P-values and a semi-volcano plot can also be shown.
##'
##' Multiple testing correction is performed using \code{qvalue} package with default parameters.
##'
##' If \code{svQTL.removal=TRUE} and svQTLs were tested in 'sqtl.seeker', gene/SNPs with
##' significant svQTL association (after multiple testing correction and similar FDR threshold)
##' are removed from the final set of significant sQTLs.
##' @title sQTLs retrieval
##' @param res.df a data.frame, output of 'sqtl.seeker' with the P-values
##' for each gene/SNP test.
##' @param FDR the False Discovery Rate to call an association significant. Default is 0.01.
##' @param md.min the minimum MD (Maximum Difference) in relative expression desired. Maximum
##' difference in relative expression (MD) gives an idea of the effect size of the association.
##' Default is 0.01.
##' @param out.pdf the name of the pdf file to create. If NULL (default), no pdf output
##' will be created. If non-NULL, the distribution of the P-values and a semi-volcano plot
##' (P-value vs MD) will be shown.
##' @param svQTL.removal if TRUE (and column 'pv.svQTL' is present in 'res.df') significant
##' sQTL which are also significant svQTLs are not reported.
##' @param FDR.svQTL the False Discovery Rate to call a svQTL, that possibly removed from the final set of sQTLs.
##' @return a subset of the input data.frame with only significant sQTLs and FDR estimates.
##' @author Jean Monlong
##' @export
sqtls <- function(res.df, FDR=.01, md.min=.01, out.pdf=NULL, svQTL.removal=TRUE, FDR.svQTL=.01){
    res.df$qv = qvalue::qvalue(res.df$pv)$qvalues

    if(!is.null(out.pdf)){
        pdf(out.pdf, 8,6)
        print(ggplot2::ggplot(res.df, ggplot2::aes(x=pv)) +
              ggplot2::geom_histogram() + ggplot2::theme_bw() +
              ggplot2::xlab("P-value") + ggplot2::ylab("number of gene/SNP")
              )
    }
    
    if(any(colnames(res.df)=="pv.svQTL")){
        res.df$qv.svQTL = qvalue::qvalue(res.df$pv.svQTL)$qvalues
        if(svQTL.removal){
            res.df = subset(res.df, qv.svQTL >= FDR.svQTL)
            if(!is.null(out.pdf)){
                print(ggplot2::ggplot(res.df, ggplot2::aes(x=pv)) +
                      ggplot2::geom_histogram() + ggplot2::theme_bw() +
                      ggplot2::xlab("P-value") + ggplot2::ylab("number of gene/SNP") +
                      ggplot2::ggtitle("After svQTL removal")
                      )
            }
        }
    }
    
    res.df = subset(res.df, qv<=FDR & md>=md.min)

    if(!is.null(out.pdf)){
        if(nrow(res.df)>0){
            print(ggplot2::ggplot(res.df, ggplot2::aes(y=pv,x=md)) +
                  ggplot2::geom_bin2d(bins=100) + ggplot2::theme_bw() +
                  ggplot2::ylab("P-value") + ggplot2::xlab("Relative expression maximum difference") + 
                  ggplot2::scale_y_log10()
                  )
        }
        dev.off()
    }
    return(res.df)
}
