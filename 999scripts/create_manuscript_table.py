#!/usr/bin/env python3
"""
Create manuscript table by expanding orthogroups with their conserved regions
"""

import re
import csv

# Your original data
original_data = [
    {'Orthogroup_ID': 'OG0000089', 'Conserved_in': 'Gram+', 'Region_analysis_rank_for_OG': 1, 
     'OG_analysis': 'Mentioned in top candidates', 'Selection_rationale': 'Conserved OG with conserved regions',
     'OG_median_in_all': '0.43%', 'OG_median_in_gram+': '0.41%', 'OG_median_in_gram-': '0.57%',
     'Location_notes': 'Cytoplasmic'},
    {'Orthogroup_ID': 'OG0000067', 'Conserved_in': 'Both', 'Region_analysis_rank_for_OG': 2,
     'OG_analysis': 'Mentioned in top candidates', 'Selection_rationale': 'Conserved OG with conserved regions',
     'OG_median_in_all': '0.41%', 'OG_median_in_gram+': '0.41%', 'OG_median_in_gram-': '0.38%',
     'Location_notes': 'Cytoplasmic'},
    {'Orthogroup_ID': 'OG0000090', 'Conserved_in': 'Gram+', 'Region_analysis_rank_for_OG': 3,
     'OG_analysis': 'Not mentioned', 'Selection_rationale': 'Conserved regions',
     'OG_median_in_all': '0.21%', 'OG_median_in_gram+': '1.20%', 'OG_median_in_gram-': '0.03%',
     'Location_notes': 'Cytoplasmic'},
    {'Orthogroup_ID': 'OG0000034', 'Conserved_in': 'Both', 'Region_analysis_rank_for_OG': 4,
     'OG_analysis': 'Mentioned in top candidates', 'Selection_rationale': 'Conserved OG with conserved regions',
     'OG_median_in_all': '0.51%', 'OG_median_in_gram+': '1.76%', 'OG_median_in_gram-': '0.15%',
     'Location_notes': 'Cytoplasmic'},
    {'Orthogroup_ID': 'OG0000019', 'Conserved_in': 'Both', 'Region_analysis_rank_for_OG': 5,
     'OG_analysis': 'Mentioned in top candidates', 'Selection_rationale': 'Conserved OG with conserved regions',
     'OG_median_in_all': '1.96%', 'OG_median_in_gram+': '1.96%', 'OG_median_in_gram-': '0.41%',
     'Location_notes': 'Cytoplasmic'},
    {'Orthogroup_ID': 'OG0000083', 'Conserved_in': 'Gram-', 'Region_analysis_rank_for_OG': 14,
     'OG_analysis': 'Mentioned in top candidates', 'Selection_rationale': 'Highest intensity of conserved in gram- and detected regions',
     'OG_median_in_all': '0.05%', 'OG_median_in_gram+': '0.05%', 'OG_median_in_gram-': '0.23%',
     'Location_notes': 'Cytoplasmic'},
    {'Orthogroup_ID': 'OG0000072', 'Conserved_in': 'Gram+', 'Region_analysis_rank_for_OG': 9,
     'OG_analysis': 'Not mentioned', 'Selection_rationale': 'Highest rank in reg. Analysis of cytoplasmic membrane proteins',
     'OG_median_in_all': '0.03%', 'OG_median_in_gram+': '0.08%', 'OG_median_in_gram-': '0.01%',
     'Location_notes': 'Cytoplasmic_membrane'},
    {'Orthogroup_ID': 'OG0000080', 'Conserved_in': 'Gram+', 'Region_analysis_rank_for_OG': 16,
     'OG_analysis': 'Not mentioned', 'Selection_rationale': '2nd highest rank in reg. Analysis of cytoplasmic membrane proteins',
     'OG_median_in_all': '0.08%', 'OG_median_in_gram+': '0.28%', 'OG_median_in_gram-': '0.01%',
     'Location_notes': 'Cytoplasmic_membrane'}
]

def parse_regions_file(filepath):
    """Parse the orthogroup_regions.txt file"""
    
    regions_by_og = {}
    current_og = None
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Split by orthogroup
    og_sections = re.split(r'=== (OG\d+) ===', content)
    
    for i in range(1, len(og_sections), 2):
        og_id = og_sections[i]
        og_content = og_sections[i+1]
        
        regions = []
        
        # Parse each region
        region_pattern = r'Region \d+:\s+Position: (\d+-\d+)\s+Length: (\d+) aa\s+Conservation: ([\d.]+)%\s+Taxa with this region: (\d+) \(([^)]+)\)\s+Sequence: ([^\n]+)'
        
        for match in re.finditer(region_pattern, og_content):
            regions.append({
                'position': match.group(1),
                'length': int(match.group(2)),
                'conservation': float(match.group(3)),
                'taxa_count': int(match.group(4)),
                'taxa_names': match.group(5),
                'sequence': match.group(6).strip()
            })
        
        regions_by_og[og_id] = regions
    
    return regions_by_og

def create_expanded_table(original_data, regions_by_og, output_file='manuscript_table.csv'):
    """Create expanded table with one row per region"""
    
    rows = []
    
    for og_data in original_data:
        og_id = og_data['Orthogroup_ID']
        
        if og_id in regions_by_og and regions_by_og[og_id]:
            # Create one row per region
            for region in regions_by_og[og_id]:
                row = {
                    'Orthogroup_ID': og_id,
                    'Position_range': region['position'],
                    'Consensus_Sequence': region['sequence'],
                    'Length_aa': region['length'],
                    'Conservation_%': f"{region['conservation']:.1f}",
                    'Taxa_Count': f"{region['taxa_count']}/8",
                    'Taxa_Names': region['taxa_names'],
                    'Conserved_in': og_data['Conserved_in'],
                    'Region_analysis_rank_for_OG': og_data['Region_analysis_rank_for_OG'],
                    'OG_analysis': og_data['OG_analysis'],
                    'Selection_rationale': og_data['Selection_rationale'],
                    'OG_median_in_all': og_data['OG_median_in_all'],
                    'OG_median_in_gram+': og_data['OG_median_in_gram+'],
                    'OG_median_in_gram-': og_data['OG_median_in_gram-'],
                    'Location_notes': og_data['Location_notes'],
                    'Jalview_notes': ''  # To be filled in manually
                }
                rows.append(row)
        else:
            # No regions found, keep original row structure
            row = {
                'Orthogroup_ID': og_id,
                'Position_range': 'No regions found',
                'Consensus_Sequence': '',
                'Length_aa': '',
                'Conservation_%': '',
                'Taxa_Count': '',
                'Taxa_Names': '',
                'Conserved_in': og_data['Conserved_in'],
                'Region_analysis_rank_for_OG': og_data['Region_analysis_rank_for_OG'],
                'OG_analysis': og_data['OG_analysis'],
                'Selection_rationale': og_data['Selection_rationale'],
                'OG_median_in_all': og_data['OG_median_in_all'],
                'OG_median_in_gram+': og_data['OG_median_in_gram+'],
                'OG_median_in_gram-': og_data['OG_median_in_gram-'],
                'Location_notes': og_data['Location_notes'],
                'Jalview_notes': ''
            }
            rows.append(row)
    
    # Write to CSV
    fieldnames = ['Orthogroup_ID', 'Position_range', 'Consensus_Sequence', 'Length_aa', 
                  'Conservation_%', 'Taxa_Count', 'Taxa_Names', 'Conserved_in',
                  'Region_analysis_rank_for_OG', 'OG_analysis', 'Selection_rationale',
                  'OG_median_in_all', 'OG_median_in_gram+', 'OG_median_in_gram-',
                  'Location_notes', 'Jalview_notes']
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"Expanded table written to: {output_file}")
    print(f"Total rows: {len(rows)}")
    
    # Print summary
    print("\nSummary by orthogroup:")
    for og_data in original_data:
        og_id = og_data['Orthogroup_ID']
        n_regions = len(regions_by_og.get(og_id, []))
        print(f"  {og_id}: {n_regions} regions")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        regions_file = sys.argv[1]
    else:
        regions_file = 'orthogroup_regions.txt'
    
    print(f"Reading regions from: {regions_file}")
    regions_by_og = parse_regions_file(regions_file)
    
    print("\nCreating expanded table...")
    create_expanded_table(original_data, regions_by_og)
    
    print("\nDone! Open manuscript_table.csv in Excel for Jalview review.")