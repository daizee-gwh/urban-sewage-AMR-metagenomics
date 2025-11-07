# Author:Daizee Talukdar

# Create the base output directory (if not exists)
mkdir -p /fastq/eggnog_output

# Define input and output locations
INPUT_DIR="/fastq/eggnog_input"
OUTPUT_BASE="/fastq/eggnog_output"
CPU=90

# Loop through all .fasta files in the input directory
for FILE in "$INPUT_DIR"/*.fasta; do
    BASENAME=$(basename "$FILE" .fasta)
    OUTPUT_DIR="$OUTPUT_BASE/$BASENAME"

    echo "Processing $BASENAME..."

    # Create output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIR"

    # Run emapper.py with override
    emapper.py -i "$FILE" \
        --itype metagenome \
        --cpu "$CPU" \
        --output "$BASENAME" \
        --output_dir "$OUTPUT_DIR" \
        --override \
        --query_cover 80 \
        --subject_cover 80 \
        --pident 80 \
        -m diamond \
        --evalue 0.001 \
        --score 60 \
        --genepred prodigal \
        --tax_scope auto \
        --target_orthologs all \
        --go_evidence non-electronic \
        --pfam_realign none \
        --report_orthologs \
        --decorate_gff yes
done

