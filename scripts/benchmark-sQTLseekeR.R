## library(sQTLseekeR)
devtools::load_all("..")

library(plyr)
library(ggplot2)

## Input files: transcript expression, gene location and genotype information
trans.exp.f = "../Data/trExp.tsv.gz"
gene.bed.f = "../Data/genes.bed"
genotype.f = "../Data/snps.tsv.gz"

## 1) Index the genotype file (if not done externally before)
genotype.indexed.f = index.genotype(genotype.f)

## 2) Prepare transcript expression
te.df = read.table(trans.exp.f, as.is=TRUE, header=TRUE, sep="\t")
samples = colnames(te.df)[-(1:2)]
colnames(te.df)[1:2] = c("trId", "geneId")

ctime.prep.te = ldply(seq(10,300,10),function(nb.samp){
    data.frame(nb.samp=nb.samp,elapsed=system.time({tre.df = prepare.trans.exp(te.df[,c("trId", "geneId", samples[1:nb.samp])])})["elapsed"])
})

glm.o = glm(elapsed~poly(nb.samp,1), data=ctime.prep.te)
ctime.prep.te$model.pred = predict(glm.o)
##save(ctime.prep.te, file="../Data/ctimePrepTe-noOpt.RData")

ggplot(ctime.prep.te, aes(x=nb.samp, y=elapsed)) + geom_point() + xlab("sample size") + ylab("elapsed time for 5 genes (s)") + theme_bw() + ggtitle("prepare.trans.exp benchmark") + geom_line(aes(y=model.pred), linetype=2)


## Dispersion benchmark
df = subset(te.df[,c("trId", "geneId", samples[1:300])], geneId == "ENSG00000153162.7")
ctime.disp = ldply(seq(10,300,40),function(nb.samp){
    st = system.time({disp = te.dispersion(hellingerDist(df[,samples[1:nb.samp]]))})
    data.frame(nb.samp=nb.samp,elapsed=st["elapsed"], disp=disp)
})

ggplot(ctime.disp, aes(x=nb.samp, y=elapsed)) + geom_point() + xlab("sample size") + ylab("elapsed time for one gene (s)") + theme_bw() + ggtitle("dispersion computation benchmark")
ggplot(ctime.disp, aes(x=nb.samp, y=disp)) + geom_point() + xlab("sample size") + ylab("elapsed time for one gene (s)") + theme_bw() + ggtitle("dispersion computation benchmark")
