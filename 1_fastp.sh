#Author:Daizee Talukdar
mkdir trimmed_reads
for a in *_R1_001.fastq.gz
do mkdir trimmed_reads/$(basename $a _R1_001.fastq.gz)
fastp -i $a -o trimmed_reads/$(basename $a _R1_001.fastq.gz)_trimmed_R1_001.fastq.gz \
  -I $(basename $a _R1_001.fastq.gz)_R2_001.fastq.gz \
  -O trimmed_reads/$(basename $a _R1_001.fastq.gz)_trimmed_R2_001.fastq.gz \
  -w 96 -M 20 -5 -3 -j trimmed_reads/$(basename $a _R1_001.fastq.gz)/fastp.json -h trimmed_reads/$(basename $a _R1_001.fastq.gz)/fastp.html
done
