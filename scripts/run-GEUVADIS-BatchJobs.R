if(FALSE){ ## To upgrade to the latest version of sQTLseekeR
  devtools::install_github("jmonlong/sQTLseekeR")
}
library(BatchJobs)
library(sQTLseekeR)

trans.exp.f = "GD667.TrQuantFlux.GeneENSG.rpkm.noRepl.ProtCod.ourf.sampNames.txt.gz"
genotype.f = "snps-Geuvadis.txt"
gene.bed.f = "genes.bed.gz"

groups = read.table("sample-groups.tsv", header=TRUE, as.is=TRUE)
ceu.samples = subset(groups,group=="CEU")$sample

## 1) Index the genotype file (if not done externally before)
system("rm -rf indexGeno")
indexGeno.reg <- makeRegistry(id="indexGeno", seed=123, file.dir="indexGeno")
batchMap(indexGeno.reg, index.genotype,genotype.f)
submitJobs(indexGeno.reg, 1, resources=list(walltime="30:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(indexGeno.reg)
genotype.indexed.f = loadResult(indexGeno.reg,1)

## 2) Prepare transcript expression
## system("rm -rf prepTE")
prepTE.reg <- makeRegistry(id="prepTE", seed=123, file.dir="prepTE")
prepTE.f <- function(te.file, samples){
    te.df = read.table(te.file, as.is=TRUE, header=TRUE, sep="\t")
    colnames(te.df)[1:2] = c("trId", "geneId")
    te.df = te.df[,c("trId", "geneId", samples)]
    library(sQTLseekeR)
    prepare.trans.exp(te.df)
}
batchMap(prepTE.reg, prepTE.f,trans.exp.f, more.args=list(samples=ceu.samples))
submitJobs(prepTE.reg, 1, resources=list(walltime="6:0:0", cores="1",queue="short"), wait=function(retries) 100, max.retries=10)
showStatus(prepTE.reg)

## 1) and 2) can be run in parallel

## 3) Run sQTLseekeR call
##  system("rm -rf sQTL")
sQTL.reg <- makeRegistry(id="sQTL", seed=123, file.dir="sQTL")
#### Eventually this part could be run on an interactive node instead of a login node
gene.bed = read.table(gene.bed.f, as.is=TRUE, sep="\t")
colnames(gene.bed) = c("chr","start","end","geneId")
tre.df = loadResult(prepTE.reg, 1)
gene.bed = subset(gene.bed, geneId %in% tre.df$geneId)
nb.gene.per.chunk = 200
gene.chunks = tapply(gene.bed$geneId, rep(1:ceiling(nrow(gene.bed)/nb.gene.per.chunk),each=nb.gene.per.chunk)[1:nrow(gene.bed)], identity)
sQTL.f <- function(chunk.id, imF){
    load(imF)
    genes = gene.chunks[[chunk.id]]
    tre.df = subset(tre.df, geneId %in% genes)
    gene.bed = subset(gene.bed, geneId %in% genes)
    library(sQTLseekeR)
    sqtl.seeker(tre.df, genotype.indexed.f, gene.bed, svQTL=TRUE)
}
imF = "sQTL-BJ-temp.RData"
save(gene.chunks,tre.df,genotype.indexed.f, gene.bed, file=imF)
batchMap(sQTL.reg, sQTL.f, 1:length(gene.chunks), more.args=list(imF=imF))
#### End of the eventual interactive node part
submitJobs(sQTL.reg, findNotDone(sQTL.reg), resources=list(walltime="20:0:0", nodes="1", cores="1",queue="sw"), wait=function(retries) 100, max.retries=10)
showStatus(sQTL.reg)

res.f = "Results/sQTLs-CEU-all.tsv"
if(file.exists(res.f)) file.remove(res.f)
tmp = reduceResultsList(sQTL.reg, fun=function(job, res){
  write.table(res, file=res.f, quote=FALSE, row.names=FALSE, col.names=!file.exists(res.f), append=file.exists(res.f), sep="\t")
})

## 4) Get significant sQTLs

## 5) svQTLs

## 6) Final sQTLs
