##
## Example of an analysis ran locally on a computer
##

library(sQTLseekeR)

## Input files: transcript expression, gene location and genotype information
trans.exp.f = "../Data/trExp.tsv.gz"
gene.bed.f = "../Data/genes.bed"
genotype.f = "../Data/snps.tsv.gz"

## Getting the IDs of samples in CEU population
groups = read.table("../Data/sample-groups.tsv", header=TRUE, as.is=TRUE)
ceu.samples = subset(groups,group=="CEU")$sample

## 1) Index the genotype file (if not done externally before)
genotype.indexed.f = index.genotype(genotype.f)

## 2) Prepare transcript expression
te.df = read.table(trans.exp.f, as.is=TRUE, header=TRUE, sep="\t")
colnames(te.df)[1:2] = c("trId", "geneId")
te.df = te.df[,c("trId", "geneId", ceu.samples)]
tre.df = prepare.trans.exp(te.df)

## 3) Test gene/SNP associations
gene.bed = read.table(gene.bed.f, as.is=TRUE, sep="\t")
colnames(gene.bed) = c("chr","start","end","geneId")
res.df = sqtl.seeker(tre.df, genotype.indexed.f, gene.bed)

## Optional: write a file with all Pvalues
res.f = "sQTLs-CEU-all.tsv"
write.table(res.df, file=res.f, quote=FALSE, row.names=FALSE, col.names=!file.exists(res.f), append=file.exists(res.f), sep="\t")
## Optional

## 4) Get significant sQTLs
sqtls.df = sqtls(res.df, FDR=.01, out.pdf="sQTLs-FDR01.pdf")
