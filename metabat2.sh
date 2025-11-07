#Author: Daizee Talukdar
for i in /scaffolds/*scaffold.fasta
do
  base=$(basename "$i" .scaffold.fasta)
  outdir="/scaffolds/metabat2_output/bin"
  mkdir -p "$outdir"
  metabat2 -i "$i" -o "$outdir"
done


