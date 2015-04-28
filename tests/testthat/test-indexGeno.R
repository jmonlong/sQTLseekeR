context("Index genotype file")

snp.df =  data.frame(chr=sample(1:22,1e3, replace=TRUE), start=sample.int(1e6,1e3, replace=TRUE))
snp.df$end = snp.df$start + 1
snp.df$snpId = paste0("snp",1:nrow(snp.df))
snp.df = snp.df[order(snp.df$chr, snp.df$start),]

samples = paste0("samp", 1:10)
geno.mat = matrix(sample(-1:2, 1e4, replace=TRUE), 1e3)
colnames(geno.mat) = samples
snp.df = data.frame(snp.df, as.data.frame(geno.mat))

test_that("Checks that file exists", {
  expect_error(index.genotype("fakeFile.tsv"), "File not found")
  write.table(snp.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_equal(length(index.genotype("temp.tsv")),1)
  expect_true(all(file.remove(c("temp.tsv","temp.tsv.bgz","temp.tsv.bgz.tbi"))))
})

test_that("Checks that required columns are present in the correct order", {
  bad.df = snp.df[,-1]
  write.table(bad.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_error(index.genotype("temp.tsv"), "Missing column")
  bad.df = snp.df[,-3]
  write.table(bad.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_error(index.genotype("temp.tsv"), "Missing column")
  bad.df = snp.df[,-2]
  write.table(bad.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_error(index.genotype("temp.tsv"), "Missing column")
  bad.df = snp.df[,-(2:3)]
  write.table(bad.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_error(index.genotype("temp.tsv"), "Missing column")
  bad.df = snp.df
  colnames(bad.df)[4] = "SNP"
  write.table(bad.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_error(index.genotype("temp.tsv"), "Missing column")
  bad.df = snp.df
  colnames(bad.df)[2:3] = c("end","start")
  write.table(bad.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_error(index.genotype("temp.tsv"), "order")
  file.remove("temp.tsv")
})

test_that("Read entire file", {
  write.table(snp.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_equal(length(index.genotype("temp.tsv")),1)
  in.df = read.bedix("temp.tsv.bgz")
  expect_equal(snp.df$chr, in.df$chr)
  expect_equal(snp.df[,8], in.df[,8])
  expect_true(all(file.remove(c("temp.tsv","temp.tsv.bgz","temp.tsv.bgz.tbi"))))
})

test_that("Reads subset of file", {
  write.table(snp.df, file="temp.tsv", row.names=FALSE, quote=FALSE, sep="\t")
  expect_equal(length(index.genotype("temp.tsv")),1)
  snp.df = snp.df[sort(sample.int(nrow(snp.df), 100)),]
  in.df = read.bedix("temp.tsv.bgz", snp.df)
  expect_equal(snp.df$chr, in.df$chr)
  expect_equal(snp.df[,8], in.df[,8])
  expect_true(all(file.remove(c("temp.tsv","temp.tsv.bgz","temp.tsv.bgz.tbi"))))
})
