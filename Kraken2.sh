
#Author:Daizee Talukdar
for i in *.R1.fastq.gz
do kraken2 --db /Databases/k2_standard_20221209 --threads 126 --output krakenOutput$(basename $i .R1.fastq.gz)_k2.output --paired $i $(basename $i .R1.fastq.gz).R2.fastq.gz --gzip-compressed --use-mpa-style --report krakenOutput$(basename $i .R1.fastq.gz)_k2.report --unclassified-out $(basename $i .R1.fastq.gz)_unclassified.fastq
done
