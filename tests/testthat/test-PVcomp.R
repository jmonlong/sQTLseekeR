context("Pvalue computation")

tre.mat = rmultinom(30,100, c(.1,.3,.6))/100
samples = paste0("samp",1:30)
colnames(tre.mat) = samples
tre.df = data.frame(geneId="gX",trId=paste0("t",1:3), tre.mat)
tre.dist = dist(t(tre.mat))
geno.mat = matrix(sample(0:2, 30*3, replace=TRUE), 3)
colnames(geno.mat) = samples
geno.df = data.frame(snpId=paste0("snp",1:3), geno.mat)

test_that("Null pv are null",{})

test_that("Results are reproducible",{})

test_that("Permutation and approximation works",{})

test_that("Pvalues for svQTLs work also",{})

test_that("Minimum pvalue is never below expected",{})

