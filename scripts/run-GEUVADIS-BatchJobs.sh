## Get gene location file
ANNOF=/users/rg/projects/Geuvadis/reference/gencode_v12.gtf.gz
INDAT=/users/rg/projects/Geuvadis/quantification/transcript/GD667.TrQuantFlux.GeneENSG.rpkm.noRepl.ProtCod.ourf.sampNames.txt.gz

zcat $ANNOF | perl -ne 'my @l=split(/\t/,$_);if($l[2] eq "gene"){$l[0]=~s/chr//;$_ =~/gene_id "([^"]+)"/;print $l[0]."\t".$l[3]."\t".$l[4]."\t".$1."\n";}' | gzip > Results/genes.bed

