##
#### Example of an analysis ran on a computing cluster
##
## Here BatchJobs parameters are just examples, for the CRG clusters. For a different
## cluster platform, 'queue' will change. 'node' or 'cores' may also change. These
## parameters are exactly the ones defined in the template, that has to be written
## beforehands (see 'cluster.tmpl' file).
##

library(BatchJobs)
library(sQTLseekeR)

## Input files: transcript expression, gene location and genotype information
trans.exp.f = "transExpression.tsv.gz"
gene.bed.f = "genes.bed"
genotype.f = "snps-012coded.tsv"

## 1) Index the genotype file (if not done externally before)
## system("rm -rf indexGeno-files") ## Run to clean previous computations
indexGeno.reg <- makeRegistry(id="indexGeno", seed=123)
batchMap(indexGeno.reg, index.genotype,genotype.f)
submitJobs(indexGeno.reg, 1, resources=list(walltime="3:0:0", nodes="1", cores="1",queue="rg-el6"), wait=function(retries) 100, max.retries=10)
showStatus(indexGeno.reg)
genotype.indexed.f = loadResult(indexGeno.reg,1)
## Note: If the genotype file is already indexed just put the filename in 'genotype.indexed.f'.

## 2) Prepare transcript expression
## system("rm -rf prepTE-files") ## Run to clean previous computations
prepTE.reg <- makeRegistry(id="prepTE", seed=123)
prepTE.f <- function(te.file, samples){
    library(sQTLseekeR)
    te.df = read.table(te.file, as.is=TRUE, header=TRUE, sep="\t")
    prepare.trans.exp(te.df)
}
batchMap(prepTE.reg, prepTE.f,trans.exp.f)
submitJobs(prepTE.reg, 1, resources=list(walltime="6:0:0", cores="1",queue="rg-el6"), wait=function(retries) 100, max.retries=10)
showStatus(prepTE.reg)

## 1) and 2) can be run in parallel

## 3) Test gene/SNP associations
#### Potentially, this part could be run on an interactive node instead of a login node (to be a nice user)
gene.bed = read.table(gene.bed.f, as.is=TRUE, sep="\t")
colnames(gene.bed) = c("chr","start","end","geneId")
tre.df = loadResult(prepTE.reg, 1)
gene.bed = subset(gene.bed, geneId %in% tre.df$geneId)
nb.gene.per.chunk = 30
gene.chunks = tapply(gene.bed$geneId, rep(1:ceiling(nrow(gene.bed)/nb.gene.per.chunk),each=nb.gene.per.chunk)[1:nrow(gene.bed)], identity)
imF = "sQTL-BJ-temp.RData"
save(gene.chunks,tre.df,genotype.indexed.f, gene.bed, file=imF)
#### End of the potential interactive node part
##  system("rm -rf sQTL-files") ## Run to clean previous computations
sQTL.reg <- makeRegistry(id="sQTL", seed=123)
sQTL.f <- function(chunk.id, imF){
    load(imF)
    genes = gene.chunks[[chunk.id]]
    tre.df = subset(tre.df, geneId %in% genes)
    gene.bed = subset(gene.bed, geneId %in% genes)
    library(sQTLseekeR)
    sqtl.seeker(tre.df, genotype.indexed.f, gene.bed, svQTL=TRUE)
}
batchMap(sQTL.reg, sQTL.f, 1:length(gene.chunks), more.args=list(imF=imF))
submitJobs(sQTL.reg, findNotDone(sQTL.reg), resources=list(walltime="20:0:0", cores="1",queue="rg-el6"), wait=function(retries) 100, max.retries=10)
showStatus(sQTL.reg)

## Optional: write a file with all Pvalues
res.f = "sQTLs-all.tsv"
if(file.exists(res.f)) file.remove(res.f)
tmp = reduceResultsList(sQTL.reg, fun=function(job, res){
  write.table(res, file=res.f, quote=FALSE, row.names=FALSE, col.names=!file.exists(res.f), append=file.exists(res.f), sep="\t")
})
## Optional

## 4) Get significant sQTLs
## system("rm -rf getSig-files") ## Run to clean previous computations
getSig.reg <- makeRegistry(id="getSig", seed=123)
getSig.f <- function(FDR, out.pdf, sQTL.reg){
  res.df = do.call(rbind,reduceResultsList(sQTL.reg, fun=function(job, res)res))
  library(sQTLseekeR)
  sqtls(res.df, FDR=FDR, out.pdf=out.pdf)
}
batchMap(getSig.reg, getSig.f,.01, more.args=list(out.pdf="sQTLs-FDR01.pdf", sQTL.reg=sQTL.reg))
submitJobs(getSig.reg, 1, resources=list(walltime="1:0:0", cores="1",queue="rg-el6"), wait=function(retries) 100, max.retries=10)
showStatus(getSig.reg)
sqtls.df = loadResult(getSig.reg, 1)
