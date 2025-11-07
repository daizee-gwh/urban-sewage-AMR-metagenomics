#Author:Daizee Talukdar
# Remove human dna contamination
#hisat2-build  /Databases/humanRefSeq/Homo_sapiensGRCH38

mkdir -p unaligned Aligned

# hisat2 aligner
for i in *_trimmed_R1_001.fastq.gz
do
    base=$(basename "$i" _trimmed_R1_001.fastq.gz)
    hisat2 -x /Databases/humanRefSeq/Homo_sapiensGRCH38 \
        -1 "$i" \
        -2 "${base}_trimmed_R2_001.fastq.gz" \
        --un-conc-gz "unaligned/${base}_unaligned.fastq.gz" \
        -S "Aligned/${base}.sam" \
        -p 96
done

