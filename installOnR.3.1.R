source("http://bioconductor.org/biocLite.R")
biocLite(c("Rsamtools", "qvalue", "vegan"))
install.packages("devtools")

## Older devtools
devtoolsurl = "http://cran.r-project.org/src/contrib/Archive/devtools/devtools_1.11.0.tar.gz"
install.packages(devtoolsurl, repos=NULL, type="source")

## Older Rcpp
rcppurl = "http://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_0.12.0.tar.gz"
install.packages(rcppurl, repos=NULL, type="source")

## Older dplyr
plyrurl = "http://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.4.2.tar.gz"
install.packages(plyrurl, repos=NULL, type="source")

library(devtools)
install_git("git://github.com/jmonlong/sQTLseekeR")
## Some ERROR message but it should install anyway
