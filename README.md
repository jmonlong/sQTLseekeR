sQTLseekeR
==========

sQTLseekeR is a R package to detect spicing QTLs (sQTLs), which are variants associated with change in 
the splicing pattern of a gene. Splicing patterns are modeled by the relative expression of the transcripts
of a gene.

For more information about the method and performance see article :
Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun. 
5:4698 doi: [10.1038/ncomms5698](http://www.nature.com/ncomms/2014/140820/ncomms5698/full/ncomms5698.html) (2014).

### Installation

To install the latest development version: `devtools::install_github("jmonlong/sQTLseekeR")`. 

This requires `devtools` package (more information [here](https://github.com/hadley/devtools)) 
which can be installed with `install.packages("devtools")`. 

### Analysis steps

The first step is to prepare the input data. `sQTLseekeR` requires three inputs:
* transcript expression. Column **trId** and *geneId**, corresponding to the transcript and gene ID are required. Then each column represents a sample and is filled with the expression values. Relative expression will be used hence both read counts or RPKMs works as the expression measure.
* gene location information. In a BED-like format, the range of each gene is explicitly defined in this file.
* genotype information. The genotype of each sample is coded as: 0=normal/normal; 1=normal/mutated; 2=mutated/mutated; -1=missing value. Furthermore the first four columns should gather information about the SNP: **chr**, **start**, **end** and **snpId**. Finally this file needs to be ordered per chr and start position.


### Running on computing clusters

`sQTLseekeR` can be used on a cluster using package `BatchJobs`. An example of an analysis using `BatchJobs` can
be found in folder `scripts`.


