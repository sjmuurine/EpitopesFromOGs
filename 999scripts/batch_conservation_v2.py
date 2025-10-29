#!/usr/bin/env python3
"""
Enhanced conservation analysis with:
- Taxon counting (not sequence counting)
- Gram classification breakdown
- No max-length restriction
- Proper gap handling
"""

import os
import glob
import re
from collections import Counter, defaultdict

def parse_aligned_fasta_with_metadata(fasta_file):
    """Parse aligned FASTA and extract taxon and conservation info"""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header:
                    # Save previous sequence
                    seq_data = parse_sequence_metadata(current_header, ''.join(current_seq))
                    sequences[seq_data['id']] = seq_data
                
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_header:
            seq_data = parse_sequence_metadata(current_header, ''.join(current_seq))
            sequences[seq_data['id']] = seq_data
    
    return sequences

def parse_sequence_metadata(header, sequence):
    """Extract metadata from header"""
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
    
    # Extract sequence ID
    seq_id = header.split()[0]
    
    return {
        'id': seq_id,
        'taxon': taxon,
        'conservation': conservation,
        'sequence': sequence,
        'full_header': header
    }

def calculate_conservation(sequences, gap_threshold=0.5):
    """Calculate conservation with proper gap handling"""
    if not sequences:
        return [], []
    
    seq_list = [s['sequence'] for s in sequences.values()]
    num_sequences = len(seq_list)
    alignment_length = len(seq_list[0])
    
    conservation_scores = []
    consensus = []
    
    for pos in range(alignment_length):
        column = [seq[pos] for seq in seq_list if pos < len(seq)]
        
        counts = Counter(column)
        
        # Count gaps
        gap_count = counts.get('-', 0)
        gap_fraction = gap_count / num_sequences
        
        # If too many gaps, position is not conserved
        if gap_fraction > gap_threshold:
            conservation_scores.append(0.0)
            consensus.append('-')
            continue
        
        # Remove gaps from conservation calculation
        if '-' in counts:
            del counts['-']
        
        if not counts:
            conservation_scores.append(0.0)
            consensus.append('-')
            continue
        
        most_common_aa, most_common_count = counts.most_common(1)[0]
        consensus.append(most_common_aa)
        
        # Conservation relative to all sequences
        conservation = most_common_count / num_sequences
        conservation_scores.append(conservation)
    
    return conservation_scores, consensus

def merge_nearby_regions(regions, sequences_dict, conservation_scores, consensus, 
                        merge_distance=10):
    """
    Merge conserved regions that are close together
    
    merge_distance: maximum gap (in amino acids) between regions to merge them
    """
    if len(regions) <= 1:
        return regions
    
    merged = []
    current = regions[0]
    
    for next_region in regions[1:]:
        gap = next_region['start'] - current['end']
        
        # If gap is small enough, merge the regions
        if gap <= merge_distance:
            # Extend current region to include the gap and next region
            new_end = next_region['end']
            new_length = new_end - current['start']
            
            # Recalculate average conservation for merged region
            start_idx = current['start'] - 1  # Convert to 0-indexed
            end_idx = new_end
            merged_conservation = sum(conservation_scores[start_idx:end_idx]) / new_length
            merged_sequence = ''.join(consensus[start_idx:end_idx])
            
            # Recalculate taxa for merged region
            region_taxa = set()
            for seq_data in sequences_dict.values():
                seq = seq_data['sequence']
                region_seq = seq[start_idx:end_idx] if end_idx <= len(seq) else seq[start_idx:]
                
                non_gap_count = sum(1 for aa in region_seq if aa != '-')
                if non_gap_count >= new_length * 0.5:
                    region_taxa.add(seq_data['taxon'])
            
            # Update current region with merged data
            current = {
                'start': current['start'],
                'end': new_end,
                'length': new_length,
                'avg_conservation': merged_conservation,
                'sequence': merged_sequence,
                'taxa': region_taxa,
                'n_taxa': len(region_taxa)
            }
        else:
            # Gap too large, keep current and move to next
            merged.append(current)
            current = next_region
    
    # Don't forget the last region
    merged.append(current)
    
    return merged

def find_conserved_regions(conservation_scores, consensus, sequences_dict, 
                          threshold=0.6, min_length=10):
    """
    Find conserved regions and track which taxa are present
    Fixed: Taxa are counted based on having NON-GAP sequence in the region,
    not based on having the exact consensus amino acid
    """
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
                    # Count taxa that have non-gap sequence in this region
                    region_taxa = set()
                    for seq_data in sequences_dict.values():
                        seq = seq_data['sequence']
                        region_seq = seq[start:i] if i <= len(seq) else seq[start:]
                        
                        # Taxa is present if it has mostly non-gap residues in region
                        non_gap_count = sum(1 for aa in region_seq if aa != '-')
                        if non_gap_count >= length * 0.5:  # At least 50% non-gap
                            region_taxa.add(seq_data['taxon'])
                    
                    regions.append({
                        'start': start + 1,
                        'end': i,
                        'length': length,
                        'avg_conservation': sum(conservation_scores[start:i]) / length,
                        'sequence': ''.join(consensus[start:i]),
                        'taxa': region_taxa,
                        'n_taxa': len(region_taxa)
                    })
                in_region = False
    
    if in_region:
        length = len(conservation_scores) - start
        if length >= min_length:
            region_taxa = set()
            for seq_data in sequences_dict.values():
                seq = seq_data['sequence']
                region_seq = seq[start:]
                
                non_gap_count = sum(1 for aa in region_seq if aa != '-')
                if non_gap_count >= length * 0.5:
                    region_taxa.add(seq_data['taxon'])
            
            regions.append({
                'start': start + 1,
                'end': len(conservation_scores),
                'length': length,
                'avg_conservation': sum(conservation_scores[start:]) / length,
                'sequence': ''.join(consensus[start:]),
                'taxa': region_taxa,
                'n_taxa': len(region_taxa)
            })
    
    return regions

def analyze_single_orthogroup(aligned_file, threshold=0.6, min_length=10, 
                             gap_threshold=0.5, merge_distance=10):
    """Analyze conservation for a single orthogroup with taxon tracking and region merging"""
    sequences = parse_aligned_fasta_with_metadata(aligned_file)
    
    if not sequences:
        return None
    
    # Get conservation pattern (should be same for all seqs in orthogroup)
    conservation_patterns = set(s['conservation'] for s in sequences.values())
    og_conservation = list(conservation_patterns)[0] if len(conservation_patterns) == 1 else "mixed"
    
    # Get unique taxa
    unique_taxa = set(s['taxon'] for s in sequences.values())
    
    conservation_scores, consensus = calculate_conservation(sequences, gap_threshold)
    
    avg_conservation = sum(conservation_scores) / len(conservation_scores) if conservation_scores else 0
    high_conservation_count = sum(1 for s in conservation_scores if s >= threshold)
    
    # Find initial regions
    regions = find_conserved_regions(conservation_scores, consensus, sequences,
                                    threshold, min_length)
    
    # Merge nearby regions if merge_distance > 0
    if merge_distance > 0 and len(regions) > 1:
        regions = merge_nearby_regions(regions, sequences, conservation_scores, 
                                      consensus, merge_distance)
    
    return {
        'n_sequences': len(sequences),
        'n_taxa': len(unique_taxa),
        'taxa': unique_taxa,
        'og_conservation': og_conservation,
        'alignment_length': len(conservation_scores),
        'avg_conservation': avg_conservation,
        'high_conservation_positions': high_conservation_count,
        'conserved_regions': regions
    }

def batch_analyze(alignments_dir='alignments', output_dir='conservation_results',
                 threshold=0.6, min_length=10, gap_threshold=0.5, merge_distance=10):
    """Batch analyze all orthogroups"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    aligned_files = glob.glob(os.path.join(alignments_dir, '*_aligned.fasta'))
    aligned_files.sort()
    
    print("="*70)
    print("ENHANCED CONSERVATION ANALYSIS")
    print("="*70)
    print(f"Aligned files: {len(aligned_files)}")
    print(f"Conservation threshold: {threshold:.0%}")
    print(f"Min region length: {min_length} aa")
    print(f"Gap threshold: {gap_threshold:.0%}")
    print(f"Merge distance: {merge_distance} aa")
    print("")
    
    results_summary = []
    
    for aligned_file in aligned_files:
        og_id = os.path.basename(aligned_file).replace('_aligned.fasta', '')
        print(f"Analyzing {og_id}...", end=' ')
        
        result = analyze_single_orthogroup(aligned_file, threshold, min_length, 
                                          gap_threshold, merge_distance)
        
        if result:
            print(f"✓ ({result['n_sequences']} seqs, {result['n_taxa']} taxa, "
                  f"{len(result['conserved_regions'])} regions)")
            
            # Write detailed output
            output_file = os.path.join(output_dir, f'{og_id}_conservation.txt')
            with open(output_file, 'w') as f:
                f.write(f"CONSERVATION ANALYSIS: {og_id}\n")
                f.write("="*70 + "\n\n")
                f.write(f"Orthogroup conservation: {result['og_conservation']}\n")
                f.write(f"Sequences: {result['n_sequences']}\n")
                f.write(f"Taxa represented: {result['n_taxa']}\n")
                f.write(f"Taxa: {', '.join(sorted(result['taxa']))}\n")
                f.write(f"Alignment length: {result['alignment_length']}\n")
                f.write(f"Average conservation: {result['avg_conservation']:.2%}\n")
                f.write(f"Positions ≥{threshold:.0%}: {result['high_conservation_positions']}\n\n")
                
                if result['conserved_regions']:
                    f.write(f"CONSERVED REGIONS (≥{min_length} aa, ≥{threshold:.0%}):\n")
                    f.write("="*70 + "\n\n")
                    
                    for i, region in enumerate(result['conserved_regions'], 1):
                        f.write(f"Region {i}:\n")
                        f.write(f"  Position: {region['start']}-{region['end']}\n")
                        f.write(f"  Length: {region['length']} aa\n")
                        f.write(f"  Conservation: {region['avg_conservation']:.1%}\n")
                        f.write(f"  Taxa with this region: {region['n_taxa']} ({', '.join(sorted(region['taxa']))})\n")
                        f.write(f"  Sequence: {region['sequence']}\n\n")
                else:
                    f.write("No conserved regions found.\n")
            
            results_summary.append({
                'og_id': og_id,
                'n_sequences': result['n_sequences'],
                'n_taxa': result['n_taxa'],
                'og_conservation': result['og_conservation'],
                'avg_conservation': result['avg_conservation'],
                'n_regions': len(result['conserved_regions']),
                'regions': result['conserved_regions']
            })
        else:
            print("✗ (failed)")
    
    # Write master summary
    write_master_summary(results_summary, output_dir, threshold, min_length)
    
    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE!")
    print("="*70)

def write_master_summary(results_summary, output_dir, threshold, min_length):
    """Write enhanced summary with gram breakdown"""
    
    summary_file = os.path.join(output_dir, 'SUMMARY_all_orthogroups.txt')
    
    with open(summary_file, 'w') as f:
        f.write("CONSERVATION ANALYSIS SUMMARY - ALL ORTHOGROUPS\n")
        f.write("="*70 + "\n\n")
        f.write(f"Parameters:\n")
        f.write(f"  Conservation threshold: {threshold:.0%}\n")
        f.write(f"  Minimum region length: {min_length} aa\n\n")
        
        # Sort by number of regions, then by taxa count
        results_summary.sort(key=lambda x: (x['n_regions'], x['n_taxa']), reverse=True)
        
        f.write("="*70 + "\n")
        f.write(f"{'Orthogroup':<15} {'Gram':<8} {'Taxa':<6} {'Seqs':<6} {'AvgCons':<10} {'Regions'}\n")
        f.write("="*70 + "\n")
        
        for result in results_summary:
            f.write(f"{result['og_id']:<15} {result['og_conservation']:<8} "
                   f"{result['n_taxa']:<6} {result['n_sequences']:<6} "
                   f"{result['avg_conservation']:<10.1%} {result['n_regions']}\n")
        
        f.write("="*70 + "\n\n")
        
        # Top candidates with conserved regions
        top_candidates = [r for r in results_summary if r['n_regions'] > 0]
        
        f.write(f"\nTOP CANDIDATES ({len(top_candidates)} orthogroups with conserved regions):\n")
        f.write("="*70 + "\n\n")
        
        for result in top_candidates:
            f.write(f"\n{result['og_id']} (Conserved in {result['og_conservation']}):\n")
            f.write(f"  Taxa: {result['n_taxa']}, Sequences: {result['n_sequences']}\n")
            f.write(f"  Average conservation: {result['avg_conservation']:.1%}\n")
            f.write(f"  Conserved regions: {result['n_regions']}\n")
            
            for i, region in enumerate(result['regions'], 1):
                f.write(f"    Region {i}: pos {region['start']}-{region['end']} "
                       f"({region['length']} aa, {region['avg_conservation']:.1%})\n")
                f.write(f"      Taxa ({region['n_taxa']}): {', '.join(sorted(region['taxa']))}\n")
                f.write(f"      {region['sequence']}\n")
        
        # Summary by gram classification
        f.write(f"\n{'='*70}\n")
        f.write("SUMMARY BY GRAM CONSERVATION PATTERN:\n")
        f.write("="*70 + "\n\n")
        
        gram_summary = defaultdict(lambda: {'count': 0, 'with_regions': 0, 'total_regions': 0})
        
        for result in results_summary:
            gram = result['og_conservation']
            gram_summary[gram]['count'] += 1
            if result['n_regions'] > 0:
                gram_summary[gram]['with_regions'] += 1
                gram_summary[gram]['total_regions'] += result['n_regions']
        
        for gram in ['gram+', 'gram-', 'both', 'mixed', 'unknown']:
            if gram in gram_summary:
                stats = gram_summary[gram]
                f.write(f"{gram}:\n")
                f.write(f"  Total orthogroups: {stats['count']}\n")
                f.write(f"  With conserved regions: {stats['with_regions']}\n")
                f.write(f"  Total conserved regions: {stats['total_regions']}\n\n")
    
    print(f"\nMaster summary: {summary_file}")
    
    # Print top candidates to screen
    if top_candidates:
        print(f"\n{'='*70}")
        print(f"TOP CANDIDATES ({len(top_candidates)} orthogroups)")
        print("="*70)
        
        for result in top_candidates[:10]:
            print(f"\n{result['og_id']} ({result['og_conservation']}): "
                  f"{result['n_regions']} region(s), {result['n_taxa']} taxa")
            for region in result['regions'][:2]:  # Show first 2 regions
                print(f"  {region['sequence'][:50]}{'...' if len(region['sequence']) > 50 else ''} "
                      f"({region['length']} aa, {region['n_taxa']} taxa)")
        
        if len(top_candidates) > 10:
            print(f"\n... and {len(top_candidates) - 10} more (see summary file)")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Enhanced conservation analysis')
    parser.add_argument('--alignments-dir', default='alignments')
    parser.add_argument('--output-dir', default='conservation_results')
    parser.add_argument('--threshold', type=float, default=0.6)
    parser.add_argument('--min-length', type=int, default=10)
    parser.add_argument('--gap-threshold', type=float, default=0.5)
    parser.add_argument('--merge-distance', type=int, default=10,
                       help='Merge regions within X amino acids (default: 10, use 0 to disable)')
    
    args = parser.parse_args()
    
    batch_analyze(args.alignments_dir, args.output_dir, 
                 args.threshold, args.min_length, args.gap_threshold, args.merge_distance)