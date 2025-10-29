#!/bin/bash

output_file="all_sequences_combined.fasta"
> "$output_file"

echo "Starting merge..."

shopt -s nullglob

# Process .fasta files (but NOT the output file)
for fasta in *.fasta; do
    # Skip the output file itself!
    if [ "$fasta" = "$output_file" ]; then
        continue
    fi
    
    isolate_name="${fasta%.*}"
    awk -v isolate="$isolate_name" '
        /^>/ {print $0 " [isolate:" isolate "]"; next}
        {print}
    ' "$fasta" >> "$output_file"
    echo "Processed: $fasta"
done

# Process .fa files
for fasta in *.fa; do
    isolate_name="${fasta%.*}"
    awk -v isolate="$isolate_name" '
        /^>/ {print $0 " [isolate:" isolate "]"; next}
        {print}
    ' "$fasta" >> "$output_file"
    echo "Processed: $fasta"
done

echo ""
echo "Done! Combined file created: $output_file"
ls -lh "$output_file"
grep -c "^>" "$output_file"

# Run this from the parent directory containing all fasta files