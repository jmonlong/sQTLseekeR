## Get gene location file
ANNOF=gencode_v12.gtf.gz

zcat $ANNOF | perl -ne 'my @l=split(/\t/,$_);if($l[2] eq "gene"){$l[0]=~s/chr//;$_ =~/gene_id "([^"]+)"/;print $l[0]."\t".$l[3]."\t".$l[4]."\t".$1."\n";}' | gzip > genes.bed

