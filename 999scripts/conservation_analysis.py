#!/usr/bin/env python3
"""
Analyze conservation in aligned FASTA sequences and identify conserved regions
"""

from collections import Counter
import re

def parse_aligned_fasta(fasta_file):
    """Parse aligned FASTA file"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                # Extract ID and isolate
                header = line[1:].strip()
                seq_id = header.split()[0]
                isolate = "unknown"
                if '[isolate:' in header:
                    isolate_match = re.search(r'\[isolate:([^\]]+)\]', header)
                    if isolate_match:
                        isolate = isolate_match.group(1)
                
                current_id = f"{seq_id} [{isolate}]"
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def calculate_conservation(sequences):
    """Calculate conservation score for each position"""
    
    if not sequences:
        return [], []
    
    # Get alignment length
    seq_list = list(sequences.values())
    alignment_length = len(seq_list[0])
    num_sequences = len(seq_list)
    
    conservation_scores = []
    consensus = []
    
    for pos in range(alignment_length):
        # Get all amino acids at this position
        column = [seq[pos] for seq in seq_list if pos < len(seq)]
        
        # Count occurrences
        counts = Counter(column)
        
        # Remove gaps from consideration
        if '-' in counts:
            del counts['-']
        
        if not counts:
            conservation_scores.append(0.0)
            consensus.append('-')
            continue
        
        # Most common amino acid
        most_common_aa, most_common_count = counts.most_common(1)[0]
        consensus.append(most_common_aa)
        
        # Conservation score = frequency of most common (excluding gaps)
        total_non_gap = sum(counts.values())
        conservation = most_common_count / total_non_gap if total_non_gap > 0 else 0
        conservation_scores.append(conservation)
    
    return conservation_scores, consensus

def find_conserved_regions(conservation_scores, consensus, threshold=0.8, min_length=5):
    """Find stretches of high conservation"""
    
    regions = []
    in_region = False
    start = 0
    
    for i, score in enumerate(conservation_scores):
        if score >= threshold:
            if not in_region:
                start = i
                in_region = True
        else:
            if in_region:
                length = i - start
                if length >= min_length:
                    regions.append({
                        'start': start,
                        'end': i - 1,
                        'length': length,
                        'avg_conservation': sum(conservation_scores[start:i]) / length,
                        'sequence': ''.join(consensus[start:i])
                    })
                in_region = False
    
    # Don't forget last region
    if in_region:
        length = len(conservation_scores) - start
        if length >= min_length:
            regions.append({
                'start': start,
                'end': len(conservation_scores) - 1,
                'length': length,
                'avg_conservation': sum(conservation_scores[start:]) / length,
                'sequence': ''.join(consensus[start:])
            })
    
    return regions

def analyze_alignment(fasta_file, output_file='conservation_analysis.txt', 
                     conservation_threshold=0.8, min_region_length=5):
    """Complete conservation analysis"""
    
    print("="*70)
    print("CONSERVATION ANALYSIS")
    print("="*70)
    
    # Parse alignment
    print(f"\nParsing alignment from: {fasta_file}")
    sequences = parse_aligned_fasta(fasta_file)
    print(f"Found {len(sequences)} sequences")
    
    if not sequences:
        print("ERROR: No sequences found!")
        return
    
    # Show sequence names
    print("\nSequences in alignment:")
    for i, seq_id in enumerate(list(sequences.keys())[:10], 1):
        print(f"  {i}. {seq_id}")
    if len(sequences) > 10:
        print(f"  ... and {len(sequences) - 10} more")
    
    # Get alignment length
    alignment_length = len(list(sequences.values())[0])
    print(f"\nAlignment length: {alignment_length} positions")
    
    # Calculate conservation
    print("\nCalculating conservation scores...")
    conservation_scores, consensus = calculate_conservation(sequences)
    
    # Overall statistics
    avg_conservation = sum(conservation_scores) / len(conservation_scores)
    high_conservation_count = sum(1 for s in conservation_scores if s >= conservation_threshold)
    
    print(f"\nOVERALL STATISTICS:")
    print(f"  Average conservation: {avg_conservation:.2%}")
    print(f"  Positions with ≥{conservation_threshold:.0%} conservation: {high_conservation_count} ({high_conservation_count/len(conservation_scores):.1%})")
    
    # Find conserved regions
    print(f"\nSearching for conserved regions (≥{conservation_threshold:.0%} conservation, ≥{min_region_length} aa)...")
    regions = find_conserved_regions(conservation_scores, consensus, 
                                    conservation_threshold, min_region_length)
    
    print(f"Found {len(regions)} conserved regions")
    
    # Write detailed output
    with open(output_file, 'w') as f:
        f.write("CONSERVATION ANALYSIS RESULTS\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Sequences analyzed: {len(sequences)}\n")
        f.write(f"Alignment length: {alignment_length}\n")
        f.write(f"Average conservation: {avg_conservation:.2%}\n")
        f.write(f"High conservation threshold: {conservation_threshold:.0%}\n\n")
        
        if regions:
            f.write(f"CONSERVED REGIONS (≥{min_region_length} amino acids):\n")
            f.write("="*70 + "\n\n")
            
            for i, region in enumerate(regions, 1):
                f.write(f"Region {i}:\n")
                f.write(f"  Position: {region['start']+1}-{region['end']+1} (1-indexed)\n")
                f.write(f"  Length: {region['length']} amino acids\n")
                f.write(f"  Average conservation: {region['avg_conservation']:.1%}\n")
                f.write(f"  Consensus sequence: {region['sequence']}\n")
                f.write("\n")
        else:
            f.write("No conserved regions found with current thresholds.\n")
            f.write("Try lowering the conservation threshold (e.g., 0.6 or 0.7)\n")
        
        # Position-by-position details
        f.write("\n" + "="*70 + "\n")
        f.write("POSITION-BY-POSITION CONSERVATION\n")
        f.write("="*70 + "\n\n")
        f.write("Position\tConsensus\tConservation\n")
        
        for i, (cons_aa, score) in enumerate(zip(consensus, conservation_scores), 1):
            f.write(f"{i}\t{cons_aa}\t{score:.2%}\n")
    
    print(f"\nDetailed results written to: {output_file}")
    
    # Display conserved regions
    if regions:
        print(f"\n{'='*70}")
        print("CONSERVED REGIONS:")
        print("="*70)
        
        for i, region in enumerate(regions, 1):
            print(f"\nRegion {i}:")
            print(f"  Position: {region['start']+1}-{region['end']+1}")
            print(f"  Length: {region['length']} aa")
            print(f"  Conservation: {region['avg_conservation']:.1%}")
            print(f"  Sequence: {region['sequence']}")
    else:
        print(f"\n{'='*70}")
        print("NO CONSERVED REGIONS FOUND")
        print("="*70)
        print(f"\nNo regions with ≥{conservation_threshold:.0%} conservation and ≥{min_region_length} amino acids.")
        print("\nTry running with lower thresholds:")
        print("  python3 conservation_analysis.py --threshold 0.6 --min-length 5")
    
    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print("="*70)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze conservation in aligned sequences')
    parser.add_argument('--input', default='orthogroups_aligned.fasta', 
                       help='Input aligned FASTA file')
    parser.add_argument('--output', default='conservation_analysis.txt',
                       help='Output text file')
    parser.add_argument('--threshold', type=float, default=0.8,
                       help='Conservation threshold (0-1, default: 0.8)')
    parser.add_argument('--min-length', type=int, default=5,
                       help='Minimum conserved region length (default: 5)')
    
    args = parser.parse_args()
    
    analyze_alignment(args.input, args.output, args.threshold, args.min_length)
