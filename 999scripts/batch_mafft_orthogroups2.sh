#!/bin/bash
#SBATCH --job-name=mafft_2orthoII
#SBATCH --account=project_2009813
#SBATCH --error=/scratch/project_2009813/JOHANNA/MAFFT/000logs/mafft_ortho2_%A_%a_err.txt
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=small
#SBATCH --output=mafft_orthogroups2_%j.log

module load biokit

# Create output directory for alignments
mkdir -p alignments

echo "=========================================="
echo "MAFFT Alignment of Orthogroups"
echo "=========================================="
echo "Start time: $(date)"
echo ""

# Counter
total=0
success=0
failed=0

# Loop through all orthogroup FASTA files
for fasta in orthogroups/_OG*.fasta; do
    if [ -f "$fasta" ]; then
        # Extract orthogroup ID from filename and remove leading underscore
        og_id=$(basename "$fasta" .fasta | sed 's/^_//')
        
        # Count sequences in file
        n_seqs=$(grep -c "^>" "$fasta")
        
        echo "Processing $og_id ($n_seqs sequences)..."
        
        # Output file
        output="alignments/${og_id}_aligned.fasta"
        
        # Run MAFFT
        # Use --auto for automatic algorithm selection based on dataset size
        if mafft --auto --thread 4 "$fasta" > "$output" 2>> mafft_errors.log; then
            echo "  ✓ Alignment complete: $output"
            ((success++))
        else
            echo "  ✗ Alignment failed for $og_id"
            ((failed++))
        fi
        
        ((total++))
        echo ""
    fi
done

echo "=========================================="
echo "ALIGNMENT SUMMARY"
echo "=========================================="
echo "Total orthogroups: $total"
echo "Successful: $success"
echo "Failed: $failed"
echo "End time: $(date)"
echo "=========================================="
echo ""
echo "Aligned files are in: alignments/"
echo "Check mafft_errors.log for any error messages"