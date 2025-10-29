#!/usr/bin/env python3
"""
Split a FASTA file into separate files by orthogroup identifier
Now with proper taxon and conservation parsing
"""

import os
import re
from collections import defaultdict

def parse_header_info(header):
    """Extract taxon, conservation, and orthogroup from header"""
    # Extract taxon
    taxon = "unknown"
    taxon_match = re.search(r'TAXON=([^,\s]+)', header)
    if taxon_match:
        taxon = taxon_match.group(1)
    
    # Extract conservation pattern
    conservation = "unknown"
    if 'Conserved_in_gram-' in header:
        conservation = "gram-"
    elif 'Conserved_in_gram+' in header or 'Conserved_in_grampos' in header:
        conservation = "gram+"
    elif 'Conserved_in_both' in header:
        conservation = "both"
    
    # Extract orthogroup
    og_match = re.search(r'_OG(\d+)', header)
    og_id = og_match.group(0) if og_match else None
    
    return taxon, conservation, og_id

def parse_fasta_by_orthogroup(fasta_file, output_dir='orthogroups'):
    """Parse FASTA and split sequences by orthogroup"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    orthogroups = defaultdict(list)
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous sequence
                if current_header:
                    taxon, conservation, og_id = parse_header_info(current_header)
                    if og_id:
                        orthogroups[og_id].append({
                            'header': current_header,
                            'sequence': ''.join(current_seq),
                            'taxon': taxon,
                            'conservation': conservation
                        })
                
                # Start new sequence
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # Don't forget last sequence
        if current_header:
            taxon, conservation, og_id = parse_header_info(current_header)
            if og_id:
                orthogroups[og_id].append({
                    'header': current_header,
                    'sequence': ''.join(current_seq),
                    'taxon': taxon,
                    'conservation': conservation
                })
    
    return orthogroups

def write_orthogroup_files(orthogroups, output_dir='orthogroups'):
    """Write separate FASTA file for each orthogroup"""
    
    stats = []
    
    for og_id, sequences in sorted(orthogroups.items()):
        output_file = os.path.join(output_dir, f'{og_id}.fasta')
        
        # Get unique taxa and conservation pattern
        taxa = set(seq['taxon'] for seq in sequences)
        conservations = set(seq['conservation'] for seq in sequences)
        conservation_pattern = list(conservations)[0] if len(conservations) == 1 else "mixed"
        
        with open(output_file, 'w') as f:
            for seq_data in sequences:
                f.write(seq_data['header'] + '\n')
                # Write sequence in 60-character lines
                sequence = seq_data['sequence']
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + '\n')
        
        stats.append({
            'og_id': og_id,
            'n_sequences': len(sequences),
            'n_taxa': len(taxa),
            'conservation': conservation_pattern,
            'filename': output_file
        })
    
    return stats

def main():
    input_file = 'ConservedInGramposGramnegBoth2.fa'
    output_dir = 'orthogroups'
    
    print("="*70)
    print("SPLITTING FASTA BY ORTHOGROUP")
    print("="*70)
    
    print(f"\nReading: {input_file}")
    orthogroups = parse_fasta_by_orthogroup(input_file, output_dir)
    
    print(f"Found {len(orthogroups)} unique orthogroups")
    
    print(f"\nWriting individual FASTA files to: {output_dir}/")
    stats = write_orthogroup_files(orthogroups, output_dir)
    
    # Summary statistics
    print("\n" + "="*70)
    print("ORTHOGROUP SUMMARY")
    print("="*70)
    print(f"{'Orthogroup':<15} {'Seqs':<6} {'Taxa':<6} {'Conservation':<15} {'Filename'}")
    print("-"*70)
    
    for stat in stats:
        print(f"{stat['og_id']:<15} {stat['n_sequences']:<6} {stat['n_taxa']:<6} "
              f"{stat['conservation']:<15} {stat['filename']}")
    
    total_sequences = sum(s['n_sequences'] for s in stats)
    print("-"*70)
    print(f"{'TOTAL':<15} {total_sequences:<6} {'':<6} {'':<15} {len(stats)} files")
    
    # Conservation pattern summary
    print(f"\n{'='*70}")
    print("CONSERVATION PATTERN SUMMARY")
    print("="*70)
    
    conservation_counts = defaultdict(int)
    for stat in stats:
        conservation_counts[stat['conservation']] += 1
    
    for pattern, count in sorted(conservation_counts.items()):
        print(f"  {pattern}: {count} orthogroups")
    
    # Save summary
    summary_file = 'orthogroup_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("Orthogroup\tSequences\tTaxa\tConservation\tFilename\n")
        for stat in stats:
            f.write(f"{stat['og_id']}\t{stat['n_sequences']}\t{stat['n_taxa']}\t"
                   f"{stat['conservation']}\t{stat['filename']}\n")
    
    print(f"\nSummary saved to: {summary_file}")
    print("\n" + "="*70)
    print("DONE!")
    print("="*70)

if __name__ == "__main__":
    main()