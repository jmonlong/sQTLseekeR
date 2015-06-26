context("sQTLseekeR main workflow")

samples = paste0("samp",1:40)

## Transcript expression creation
genes = rep(paste0("g",1:10), c(1:5, 1:5))
transcripts = paste(genes,1:30 ,sep="-")
transcript.lambda = runif(length(transcripts), 1,20)
te.mat = matrix(rpois(length(samples)*length(transcripts), transcript.lambda)/10, length(transcripts))
colnames(te.mat) = samples
te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
te.df = prepare.trans.exp(te.df)

## Create genes coordinates
gene.bed = data.frame(chr=sample(1:21,12, replace=TRUE), start=sample.int(1e6,12, replace=TRUE))
gene.bed$end = gene.bed$start + round(runif(nrow(gene.bed),10,1e3))
gene.bed$geneId = paste0("g",1:12)

## Genotypes creation
snp.df = lapply(unique(te.df$geneId)[-1], function(gene.i){
  data.frame(chr=gene.bed$chr[gene.bed$geneId==gene.i], start=runif(runif(1,10,100),gene.bed$start[gene.bed$geneId==gene.i]-1e4,gene.bed$end[gene.bed$geneId==gene.i]+1e4))
})
snp.df = rbind(do.call(rbind,snp.df),data.frame(chr=sample(1:21,100, replace=TRUE), start=sample.int(1e6,100, replace=TRUE)))
snp.df$end = snp.df$start + 1
snp.df$snpId = paste0("snp",1:nrow(snp.df))
snp.df = snp.df[order(snp.df$chr, snp.df$start),]
geno.mat = matrix(sample(0:2, nrow(snp.df)*length(samples), replace=TRUE), nrow(snp.df))
geno.mat[sample.int(length(geno.mat),10)] = -1
colnames(geno.mat) = samples
snp.df = data.frame(snp.df, as.data.frame(geno.mat))
write.table(snp.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
genotype.f = index.genotype("temp.tsv")

test_that("Non-NULL output",{
  expect_true(!is.null(sqtl.seeker(te.df, genotype.f, gene.bed)))
})

file.remove(c("temp.tsv","temp.tsv.bgz","temp.tsv.bgz.tbi"))
