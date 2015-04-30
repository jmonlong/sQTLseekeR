context("Pvalue computation")

tre.mat = rmultinom(30,100, c(.1,.3,.6))/100
samples = paste0("samp",1:30)
colnames(tre.mat) = samples
tre.df = data.frame(geneId="gX",trId=paste0("t",1:3), tre.mat)
tre.dist = dist(t(tre.mat))
geno.mat = matrix(sample(0:2, 30*3, replace=TRUE), 3)
colnames(geno.mat) = samples
geno.df = data.frame(snpId=paste0("snp",1:3), geno.mat)
res.df = plyr::ldply(1:nrow(geno.df), function(ii){
  compFscore(geno.df[ii,], tre.dist, tre.df)
})

test_that("Results are reproducible",{
  res1 = compPvalue(res.df, tre.dist)
  res2 = compPvalue(res.df, tre.dist)
  expect_true(all(abs(res1$pv-res2$pv)<.1))
})

test_that("Permutation and approximation give similar pvalues",{
  res.perm = compPvalue(res.df, tre.dist, approx=FALSE)
  res.approx = compPvalue(res.df, tre.dist, approx=TRUE)
  expect_true(all(abs(res.perm$pv-res.approx$pv)<.1))
})

test_that("Pvalues for svQTLs work also",{
  res.df = plyr::ldply(1:nrow(geno.df), function(ii){
    compFscore(geno.df[ii,], tre.dist, tre.df, svQTL=TRUE)
  })
  res.df = compPvalue(res.df, tre.dist, svQTL = TRUE, nb.perm.max=200)
  expect_true(!is.null(res.df$pv.svQTL))
})

test_that("Minimum pvalue is never below expected",{
  res.df$F[3] = res.df$F[3] * 1000
  res.df = compPvalue(res.df, tre.dist, min.nb.ext.scores = 1, nb.perm.max = 10)
  expect_true(all(res.df$pv>1/100))
  res.df = compPvalue(res.df, tre.dist, min.nb.ext.scores = 100, nb.perm.max = 1000)
  expect_true(all(res.df$pv>1/10000))
  expect_true(all(res.df$pv[3]<1/100))
  res.df = compPvalue(res.df, tre.dist, min.nb.ext.scores = 1000, nb.perm.max = 10000)
  expect_true(all(res.df$pv>1/100000))
  expect_true(all(res.df$pv[3]<1/1000))
})

