context("F score computation for one test")

tre.mat = rmultinom(30,100, c(.1,.3,.6))/100
samples = paste0("samp",1:30)
colnames(tre.mat) = samples
tre.df = data.frame(geneId="gX",trId=paste0("t",1:3), tre.mat)
tre.dist = dist(t(tre.mat))
geno.mat = matrix(sample(0:2, 30*3, replace=TRUE), 3)
colnames(geno.mat) = samples
geno.df = data.frame(snpId=paste0("snp",1:3), geno.mat)

test_that("Only one genotype is present", {
  geno.df = rbind(geno.df, geno.df)
  expect_error(compFscore(subset(geno.df, snpId=="snp1"), tre.dist, tre.df), "SNP is duplicated")
})

test_that("Unknown genotypes are removed", {
  f1 = compFscore(geno.df[1,], tre.dist, tre.df)
  expect_true(nrow(f1)>0)
  geno.df[1,5:9] = -1
  f2 = compFscore(geno.df[1,], tre.dist, tre.df)
  expect_true(f1$F != f2$F)
})

test_that("Results are reproducible", {
  f1 = compFscore(geno.df[1,], tre.dist, tre.df)
  f2 = compFscore(geno.df[1,], tre.dist, tre.df)
  f3 = compFscore(geno.df[1,], tre.dist, tre.df)
  expect_true(all(c(f1$F, f2$F, f3$F) == f1$F))
})

test_that("Inconsistent sample names throw an error", {
  colnames(geno.df)[-1] = paste0("s",1:30)
  expect_error(compFscore(geno.df[1,], tre.dist, tre.df), "No common samples")
})

test_that("Distance object as input", {
  expect_error(compFscore(geno.df[1,], as.matrix(tre.dist), tre.df), "must be a distance")  
  expect_true(nrow(compFscore(geno.df[1,], tre.dist, tre.df))>0)
})

test_that("svQTL can be computed", {
  f1 = compFscore(geno.df[1,], tre.dist, tre.df)
  expect_true(nrow(f1)>0)
  expect_true(!is.null(f1$F))
  expect_true(!is.na(f1$F))
  expect_true(is.null(f1$F.svQTL))
  f1 = compFscore(geno.df[1,], tre.dist, tre.df, svQTL=TRUE)
  expect_true(nrow(f1)>0)
  expect_true(!is.null(f1$F.svQTL))
  expect_true(!is.na(f1$F.svQTL))
  expect_true(!is.null(f1$F))
  expect_true(!is.na(f1$F))
})

test_that("Other metrics are computed", {
  f1 = compFscore(geno.df[1,], tre.dist, tre.df)
  expect_true(nrow(f1)>0)
  expect_true(!is.null(f1$F))
  expect_true(!is.na(f1$F))
  expect_true(!is.null(f1$md))
  expect_true(!is.na(f1$md))
  expect_true(!is.null(f1$nb.groups))
  expect_true(!is.na(f1$nb.groups))
  expect_true(!is.null(f1$tr.first))
  expect_true(!is.na(f1$tr.first))
  expect_true(!is.null(f1$tr.second))
  expect_true(!is.na(f1$tr.second))
})

