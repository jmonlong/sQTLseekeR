sQTLseekeR
==========

sQTLseekeR is a R package to detect spicing QTLs (sQTLs), which are variants associated with change in the splicing pattern of a gene. Splicing patterns are modeled by the relative expression of the transcripts of a gene.

For more information about the method and performance see article :
Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun. 
5:4698 doi: [10.1038/ncomms5698](http://www.nature.com/ncomms/2014/140820/ncomms5698/full/ncomms5698.html) (2014).

[TOC]

### Installation

To install the latest development version: `devtools::install_github("jmonlong/sQTLseekeR")`. 

This requires `devtools` package (more information [here](https://github.com/hadley/devtools)) which can be installed with `install.packages("devtools")`. 

### Analysis steps

The first step is to prepare the input data. `sQTLseekeR` requires three inputs:
* transcript expression. Column *trId* and *geneId*, corresponding to the transcript and gene ID are required. Then each column represents a sample and is filled with the expression values. Relative expression will be used hence both read counts or RPKMs works as the expression measure.
* gene location information. In a BED-like format, the range of each gene is explicitly defined in this file.
* genotype information. The genotype of each sample is coded as: 0=normal/normal; 1=normal/mutated; 2=mutated/mutated; -1=missing value. Furthermore the first four columns should gather information about the SNP: *chr*, *start*, *end* and *snpId*. Finally **this file needs to be ordered** per chr and start position.

When all input files are correctly formatted the `sQTLseekeR` prepares the data through functions `prepare.trans.exp` and `index.genotype`.
* `prepare.trans.exp` will :
  * remove transcripts with low expression.
  * remove genes with less than two expressed transcript.
  * remove genes with low splicing dispersion.
  * remove genes with not enough different splicing patterns.
  * flag samples with low gene expression.
* `index.genotype` compresses and indexes the genotype file to optimize further accession of particular regions.

Once the imput files are ready, `sqtl.seeker` function will compute the P-values for each pair of gene/SNP testing the association between the genotype and transcript relative expression. Here is a quick description of the parameters that would most likely be tweaked:
* `genic.window` the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
* `svQTL` should svQTLs test be performed in addition to sQTLs (default is FALSE). svQTLs are used to identify potential false positive among the significant sQTLs. svQTLs represents situtation where the variance in transcript relative expression is different between genotype groups. In this particular situation identification of sQTLs is less robust as we assume homogenetiy of the variance between groups, hence it might be safer to remove svQTLs from the list of reported sQTLs. However computation of svQTLs cannot rely on an asymptotic approximation, hence the heavy permutations will considerably increase the running time. 
* `nb.perm.max` the maximum number of permutation/simulation to compute the P-value. The higher this number, the lower the P-values can potentially get but the longer the computatio (especially relevant when `svQTL=TRUE`).


### Running on computing clusters

`sQTLseekeR` can be used on a cluster using package `BatchJobs`. An example of an analysis using `BatchJobs` can
be found in folder `scripts`.


