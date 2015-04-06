context("Index genotype file")

bin.df =  data.frame(chr=sample(1:22,1e3, replace=TRUE), start=sample.int(1e6,1e3, replace=TRUE))
bin.df$end = bin.df$start + 1
bin.df = bin.df[order(bin.df$chr, bin.df$start),]

samples = paste0("samp", 1:10)
geno.mat = matrix(sample(-1:2, 1e4, replace=TRUE), 1e3)
colnames(geno.mat) = samples
bin.df = data.frame(bin.df, as.data.frame(geno.mat))

test_that("Checks that file exists", {
  expect_error(index.genotype("fakeFile.tsv"), "File not found")
  write.table(bin.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_equal(length(index.genotype("temp.tsv")),1)
  expect_true(all(file.remove(c("temp.tsv","temp.tsv.bgz","temp.tsv.bgz.tbi"))))
})

test_that("Read entire file", {
  write.table(bin.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_equal(length(index.genotype("temp.tsv")),1)
  in.df = read.bedix("temp.tsv.bgz")
  expect_equal(bin.df$chr, in.df$chr)
  expect_equal(bin.df[,8], in.df[,8])
  expect_true(all(file.remove(c("temp.tsv","temp.tsv.bgz","temp.tsv.bgz.tbi"))))
})

test_that("Reads subset of file", {
  write.table(bin.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_equal(length(index.genotype("temp.tsv")),1)
  bin.df = bin.df[sort(sample.int(nrow(bin.df), 100)),]
  in.df = read.bedix("temp.tsv.bgz", bin.df)
  expect_equal(bin.df$chr, in.df$chr)
  expect_equal(bin.df[,8], in.df[,8])
  expect_true(all(file.remove(c("temp.tsv","temp.tsv.bgz","temp.tsv.bgz.tbi"))))
})
