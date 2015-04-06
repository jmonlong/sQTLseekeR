context("Prepare transcript expression")

genes = rep(paste0("g",1:10), c(1:5, 1:5))
transcripts = paste(genes,1:30 ,sep="-")
samples = paste0("samp",1:40)
transcript.lambda = runif(length(transcripts), 1,20)
te.mat = matrix(rpois(length(samples)*length(transcripts), transcript.lambda)/10, length(transcripts))
colnames(te.mat) = samples

test_that("Checks for missing columns", {
  te.df = data.frame(gene = genes, trans=transcripts, as.data.frame(te.mat))
  expect_error(prepare.trans.exp(te.df), "Missing column")
  te.df = data.frame(geneId = genes, tr=transcripts, as.data.frame(te.mat))
  expect_error(prepare.trans.exp(te.df), "Missing column")
  te.df = data.frame(geneid = genes, trId=transcripts, as.data.frame(te.mat))
  expect_error(prepare.trans.exp(te.df), "Missing column")
})

test_that("A subset of the input data is returned", {
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df)
  expect_true(nrow(te.p)>0)
  expect_true(nrow(te.p)<=nrow(te.df))
})

test_that("Transcripts with low expression are filtered", {
  te.mat[4,] = 0.001
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df)
  expect_true(all(te.p$trId != te.df$trId[4]))  
  te.mat[4,1] = 0.1
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df, min.dispersion = 0)
  expect_true(any(te.p$trId == te.df$trId[4]))  
  te.p = prepare.trans.exp(te.df, min.transcript.exp = .2)
  expect_true(all(te.p$trId != te.df$trId[4]))  
})

test_that("Genes with low dispersion are filtered", {
  te.mat[4:6,] = runif(3*length(samples), 2, 2.1)
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df)
  expect_true(all(te.p$geneId != te.df$geneId[4]))  
  te.p = prepare.trans.exp(te.df, min.dispersion = .0001)
  expect_true(any(te.p$geneId == te.df$geneId[4]))  
  te.mat[4:6,] = runif(3*length(samples), 2, 6)
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df)
  expect_true(any(te.p$geneId == te.df$geneId[4]))  
})

test_that("Genes with low expression are filtered", {
  te.mat[4:6,1] = .001
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df)
  expect_true(is.na(subset(te.p, geneId == te.df$geneId[4])[1,samples[1]]))
})

test_that("Single-transcript genes are removed", {
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df)
  expect_true(all(te.p$trId != te.df$trId[1]))  
})

test_that("Output are ratios", {
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  te.p = prepare.trans.exp(te.df)
  expect_true(all(as.numeric(unlist(te.p[,samples])) <=1))
  expect_true(all(as.numeric(unlist(te.p[,samples])) >=0))
  expect_true(all(round(apply(te.p[which(te.p$geneId==te.p$geneId[1]),samples], 2, sum, na.rm=TRUE), 3)==1))
})

test_that("Stops if not enough samples", {
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat[,1:3]))
  expect_error(prepare.trans.exp(te.df), "Not enough samples")
})

test_that("Error if no suitable gene", {
  te.df = data.frame(geneId = genes, trId=transcripts, as.data.frame(te.mat))
  expect_error(prepare.trans.exp(te.df, min.dispersion=Inf), "No genes")
})


