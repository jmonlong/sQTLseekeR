sQTLseekeR
==========

sQTLseekeR is a R package to detect splicing QTLs (sQTLs), which are variants associated with change in the splicing pattern of a gene. Here, splicing patterns are modeled by the relative expression of the transcripts of a gene.

For more information about the method and performance see article :
Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun. 
5:4698 doi: [10.1038/ncomms5698](http://www.nature.com/ncomms/2014/140820/ncomms5698/full/ncomms5698.html) (2014).

### Installation

To install the latest development version: `devtools::install_github("jmonlong/sQTLseekeR")`. 

This requires `devtools` package (more information [here](https://github.com/hadley/devtools)) 
which can be installed with `install.packages("devtools")`. 

It also requires R 3.1 or higher. 

### Analysis steps

The first step is to prepare the input data. `sQTLseekeR` requires three inputs:
* transcript expression. Column *trId* and *geneId*, corresponding to the transcript and gene ID are required. Then each column represents a sample and is filled with the expression values. Relative expression will be used hence both **read counts or RPKMs** works as the expression measure. However, it is not recommended to use transformed (log, quantile, or any non-linear tranformation) counts or RPKMs because Hellinger distance is suited for relative expression.
* gene location information. In a BED-like format, the range of each gene is explicitly defined in this file.
* genotype information. The genotype of each sample is coded as follow: 0 for ref/ref; 1 for ref/mutated; 2 for mutated/mutated; -1 for missing value. Furthermore the first four columns should gather information about the SNP: *chr*, *start*, *end* and *snpId*. Finally **this file needs to be ordered** per *chr* and *start* position.

When all input files are correctly formatted `sQTLseekeR` prepares the data through functions `prepare.trans.exp` and `index.genotype`.
* `prepare.trans.exp` will :
  * remove transcripts with low expression.
  * remove genes with less than two expressed transcript.
  * remove genes with low splicing dispersion.
  * remove genes with not enough different splicing patterns.
  * flag samples with low gene expression.
* `index.genotype` compresses and indexes the genotype file to optimize further accession of particular regions.

Once the input files are ready, `sqtl.seeker` function will compute the P-values for each pair of gene/SNP testing the association between the genotype and transcript relative expression. Here is a quick description of the parameters that would most likely be tweaked:
* `genic.window` the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
* `svQTL` should svQTLs test be performed in addition to sQTLs (default is FALSE). svQTLs are used to identify potential false positive among the significant sQTLs. svQTLs represents situation where the variance in transcript relative expression is different between genotype groups. In this particular situation identification of sQTLs is less robust as we assume homogeneity of the variance between groups, hence it might be safer to remove svQTLs from the list of reported sQTLs. However computation of svQTLs cannot rely on an asymptotic approximation, hence the heavy permutations will considerably increase the running time. 
* `nb.perm.max` the maximum number of permutation/simulation to compute the P-value. The higher this number, the lower the P-values can potentially get but the longer the computation (especially relevant when `svQTL=TRUE`).

Finally, function `sqtls` is used to retrieve significant associations. The user can manually define a false discovery rate(FDR) or perform further filtering afterwards. Of note, there is a separate FDR threshold for svQTL removal (if svQTLs were computed), which is usually preferred to be low (e.g. around 0.01).

An example of an analysis can be found in folder `scripts`.

### Running on computing clusters

`sQTLseekeR` can be used on a cluster using package `BatchJobs`. An example of an analysis using `BatchJobs` can
be found in folder `scripts`.

`BatchJobs` is a potent package but basic functions are enough in our situation. Here is a quick practical summary of `BatchJobs` commands used in the script:
* `makeRegistry` creates a registry used to manipulate jobs for a particular analysis step.
* `batchMap` adds jobs to a registry. Simply, the user gives a function and a list of parameters. One job per parameter will be created to compute the output of the function using this specific parameter.
* `submitJobs` submits the jobs to the cluster. This is where the queue, ,maximum computation time, number of core can be specified. Moreover, if needed, a subset of the jobs can be sent to the cluster. Functions `findNotDone` and `findErrors` are particularly useful to find which the jobs that didn't finish or were lost in the limbo of the cluster management process.
* `showStatus` outputs the status of the computations.
* `loadResult` retrieves the output of one specific job, while `reduceResultsList` retrieves output for all jobs into a list format.

Another important point about `BatchJobs` is its configuration for the computing cluster. An example of the configuration files can be found in the `scripts` folder:
* If present in the working directory, `.BatchJobs.R` is loaded when the `BatchJobs` package is loaded. It defines which template to use and `BatchJobs` functions. In practice, it loads another R script file (here `makeClusterFunctionsAdaptive.R`) with the functions to use. In `.BatchJobs.R` users would only need to change the email address where to send the log messages to.
* In `makeClusterFunctionsAdaptive.R`, users just need to check/replace `qsub`/`qdel`/`qstat` calls with the correct bash commands (sometimes `msub`/`canceljob`/`showq`). This file should also be in the working directory when `BatchJobs` is loaded.
* Finally `cluster.tmpl` is a template form of a job bash script that would be send to the cluster. There the correct syntax for the resources or parameters of the cluster are defined. This file should also be in the working directory when `BatchJobs` is loaded.
