#Author:Daizee Talukdar
for i in *.1.fastq.gz
do spades.py --meta -1 $i -2 $(basename $i .1.fastq.gz).2.fastq.gz --threads 126 -o $(basename $i .1.fastq.gz) 
done
