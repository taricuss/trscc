
import pandas as pd
from Bio.Seq import Seq

# Load data
golden_list = pd.read_csv("results/golden_list_final_full_scores_with_rlt_real_trcss_real_chromatin.csv")
flanking = pd.read_csv("results/flanking_sequences_200bp.csv")
merged = pd.merge(golden_list, flanking, on=['gene', 'chromosome', 'position'], suffixes=('', '_flanking'))

# Process each variant to create editseq
batch_data = []
seen = set()  # Track unique sequence_name
for idx, row in merged.iterrows():
    # Get info
    gene = row['gene']
    hgvs = row['hgvs_c']
    ref = row['ReferenceAlleleVCF'] if pd.notna(row['ReferenceAlleleVCF']) else ''
    alt = row['AlternateAlleleVCF'] if pd.notna(row['AlternateAlleleVCF']) else ''
    flanking_seq = row['flanking_sequence']
    
    # Determine edit type
    if 'del' in hgvs and 'ins' not in hgvs:  # Deletion
        # For deletions: (X/-) where X is the deleted sequence
        # Find the deleted part: ref is longer than alt
        deleted = ref.replace(alt, '') if len(ref) > len(alt) else ref
        # Our flanking sequence is 100bp upstream, then ref, then 100bp downstream
        upstream = flanking_seq[:100]
        downstream = flanking_seq[100+len(ref):]
        editseq = f"{upstream}({deleted}/-){downstream}"
    elif 'ins' in hgvs:  # Insertion
        # For insertions: (-+X) where X is inserted
        inserted = alt.replace(ref, '') if len(alt) > len(ref) else alt
        upstream = flanking_seq[:100]
        downstream = flanking_seq[100:]
        editseq = f"{upstream}(-+{inserted}){downstream}"
    else:  # Assume other types (not present)
        continue
    
    # Create batch entry
    sequence_name = f"{gene}_{hgvs.split(':')[-1].split(' ')[0]}"
    if sequence_name not in seen:
        seen.add(sequence_name)
        batch_data.append({
            'sequence_name': sequence_name,
            'editseq': editseq
        })

# Save batch file
batch_df = pd.DataFrame(batch_data)
batch_df.to_csv("data/pridict2/input/batch_planti.csv", index=False)
print("Saved PRIDICT2 batch file to data/pridict2/input/batch_planti.csv")
