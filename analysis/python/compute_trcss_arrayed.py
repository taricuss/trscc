
import gzip
import pandas as pd
from pathlib import Path
import numpy as np

def load_okseq_data(file_path):
    print(f"Loading OK-seq data from {file_path}...")
    df = pd.read_csv(
        file_path,
        sep="\t",
        compression="gzip"
    )
    df['chr'] = df['chr'].astype(str).str.replace('chr', '')
    return df

def compute_okseq_score(df, chrom, pos, window=50000):
    chrom = str(chrom).replace('chr', '')
    chrom_df = df[df['chr'] == chrom].copy()
    if len(chrom_df) == 0:
        return 0.0
    window_df = chrom_df[
        (chrom_df['pos'] >= pos - window) & 
        (chrom_df['pos'] <= pos + window)
    ]
    if len(window_df) == 0:
        return 0.0
    avg_w = window_df['w'].mean()
    avg_c = window_df['c'].mean()
    total = avg_w + avg_c
    if total == 0:
        return 0.0
    fd = (avg_w - avg_c) / total
    return fd

if __name__ == "__main__":
    print("=== Computing TRCSS for Arrayed Editing Data ===")
    
    base_dir = Path(__file__).parent
    
    # Load all necessary files
    arrayed = pd.read_csv(base_dir / "results/arrayed_editing_pe.csv")
    tx_dir = pd.read_csv(base_dir / "results/arrayed_editing_tx_dir.csv")
    
    # Drop duplicates from tx_dir to avoid multiplying rows
    tx_dir_unique = tx_dir.drop_duplicates(subset=['chromosome', 'position'])
    
    # Merge to get tx_dir, keeping original arrayed order
    merged = pd.merge(arrayed, tx_dir_unique, on=['chromosome', 'position'], how='left')
    print("Merged data shape:", merged.shape)
    
    # Load OK-seq data
    okseq_file = base_dir / "data/geo_data/GSE114017/GSM3130725_REP1_rpe_edu.txt.gz"
    okseq_df = load_okseq_data(okseq_file)
    
    # Compute RFD (okseq_fd) and TRCSS for each row
    rfd_list = []
    trcss_list = []
    for _, row in merged.iterrows():
        chrom = row['chromosome']
        pos = row['position']
        tx_dir = row['transcription_direction']
        
        # Compute RFD
        rfd = compute_okseq_score(okseq_df, chrom, pos)
        rfd_list.append(rfd)
        
        # Compute TRCSS
        trcss = (1 - rfd * tx_dir) / 2
        trcss_list.append(trcss)
    
    # Add to merged dataframe
    merged['RFD'] = rfd_list
    merged['TRCSS'] = trcss_list
    
    print("\n=== Summary ===")
    print("Number of variants with non-zero RFD:", (merged['RFD'] != 0).sum())
    print("TRCSS summary stats:")
    print(merged['TRCSS'].describe())
    
    # Save the final data
    output_path = base_dir / "results/arrayed_editing_full_trcss.csv"
    merged.to_csv(output_path, index=False)
    print("\nSaved full data to", output_path)
