##
#### Example of an analysis ran locally on a computer
##
## Due to the amount of genes/snps tested for a typical QTL analysis, it is recommended
## to use other computing resources. A parctical example, using BatchJobs package, is
## presented in 'run-example-BatchJobs.R' file.
##
## This example is here as a clearer summary of the workflow and for debugging purpose
## to check is the package is correctly installed and the inputs correctly formatted.
##

library(sQTLseekeR)

## Input files: transcript expression, gene location and genotype information
trans.exp.f = "transExpression.tsv.gz"
gene.bed.f = "genes.bed"
genotype.f = "snps-012coded.tsv"

## 1) Index the genotype file (if not done externally before)
genotype.indexed.f = index.genotype(genotype.f)

## 2) Prepare transcript expression
te.df = read.table(trans.exp.f, as.is=TRUE, header=TRUE, sep="\t")
tre.df = prepare.trans.exp(te.df)

## 3) Test gene/SNP associations
gene.bed = read.table(gene.bed.f, as.is=TRUE, sep="\t")
colnames(gene.bed) = c("chr","start","end","geneId")
res.df = sqtl.seeker(tre.df, genotype.indexed.f, gene.bed, svQTL=TRUE)

## Optional: write a file with all Pvalues
res.f = "sQTLs-all.tsv"
write.table(res.df, file=res.f, quote=FALSE, row.names=FALSE, sep="\t")
## Optional

## 4) Get significant sQTLs
sqtls.df = sqtls(res.df, FDR=.01, out.pdf="sQTLs-FDR01.pdf")
