import pandas as pd
import numpy as np
import gtfparse
import os
import pyranges as pr

# =====================
# Configuration
# =====================
BIN_SIZE = 32  # DO NOT CHANGE (model resolution)
GTF_DIR = "/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out"
OUTPUT_DIR = f"{GTF_DIR}/multi_task_targets"
DOMINANT_GTF_DIR = f"{GTF_DIR}/dominant_isoforms"
REFERENCE_GTF = "/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"

CELL_TYPES = ['B_Cell', 'T_Cell', 'Macrophage']
GTF_PATHS = {
    ct: os.path.join(GTF_DIR, f"{ct}_aggregate.gtf")
    for ct in CELL_TYPES
}

CLASS_MAP = {
    'desert': 0,
    'intron': 1,
    'exon': 2,
    'UTR': 3
}

PRIORITY_MAP = {
    'exon': 3,
    'UTR': 2,
    'intron': 1,
    'desert': 0
}

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(DOMINANT_GTF_DIR, exist_ok=True)

print(f"[INFO] Target bin size: {BIN_SIZE} bp")

# =====================
# Genome boundaries
# =====================
print("[INFO] Loading reference GTF to determine chromosome lengths...")
ref_df = gtfparse.read_gtf(REFERENCE_GTF).to_pandas()[['seqname', 'end']]

chrom_lengths = (
    ref_df
    .groupby('seqname', observed=True)['end']
    .max()
    .astype(int)
    .to_dict()
)

print(f"[INFO] Found {len(chrom_lengths)} chromosomes.")

# =====================
# Helpers
# =====================
def generate_bins_for_chrom(chrom, chrom_len, bin_size):
    starts = np.arange(1, chrom_len + 1, bin_size)
    ends = np.minimum(starts + bin_size - 1, chrom_len)

    return pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": ends
    })


def get_dominant_features(gtf_path, cell_type):
    print(f"[INFO] Processing {cell_type} GTF")

    full_gtf = gtfparse.read_gtf(gtf_path).to_pandas()
    full_gtf['TPM'] = pd.to_numeric(full_gtf['TPM'], errors='coerce')

    tx_df = full_gtf[
        (full_gtf['feature'] == 'transcript') &
        (full_gtf['TPM'] > 0)
    ].copy()

    dominant_tx = tx_df.loc[
        tx_df.groupby('gene_id', observed=True)['TPM'].idxmax()
    ]

    dominant_tx_ids = dominant_tx['transcript_id'].tolist()
    print(f"\t- Identified {len(dominant_tx_ids)} dominant transcripts")

    feats = full_gtf[
        full_gtf['transcript_id'].isin(dominant_tx_ids) &
        full_gtf['feature'].isin(['exon', 'five_prime_utr', 'three_prime_utr'])
    ].copy()

    feats['class_feature'] = feats['feature'].map({
        'exon': 'exon',
        'five_prime_utr': 'UTR',
        'three_prime_utr': 'UTR'
    })

    # ---- derive introns ----
    exons = feats[feats['feature'] == 'exon'].sort_values(
        ['transcript_id', 'start']
    )

    introns = []
    for tx_id, grp in exons.groupby('transcript_id', observed=True):
        next_start = grp['start'].shift(-1)
        intron_start = grp['end'] + 1
        intron_end = next_start - 1

        valid = intron_start < intron_end
        if valid.any():
            intr = grp.loc[valid].copy()
            intr['start'] = intron_start[valid]
            intr['end'] = intron_end[valid]
            intr['class_feature'] = 'intron'
            introns.append(intr)

    if introns:
        feats = pd.concat([feats, pd.concat(introns)], ignore_index=True)

    feats = feats[[
        'seqname', 'start', 'end', 'class_feature'
    ]].rename(columns={'seqname': 'chrom'}).drop_duplicates()

    return feats, dominant_tx_ids


def write_dominant_gtf(src_path, keep_tx_ids, dest_path):
    keep = set(keep_tx_ids)
    with open(src_path) as inp, open(dest_path, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            if 'transcript_id "' not in line:
                continue
            tx = line.split('transcript_id "')[1].split('"')[0]
            if tx in keep:
                out.write(line)


# =====================
# Load dominant features (once per cell type)
# =====================
features_by_cell = {}

for ct in CELL_TYPES:
    gtf = GTF_PATHS[ct]
    if not os.path.exists(gtf):
        print(f"[WARN] Missing GTF for {ct}, skipping.")
        continue

    feats, tx_ids = get_dominant_features(gtf, ct)

    dom_gtf = os.path.join(DOMINANT_GTF_DIR, f"{ct}_dominant.gtf")
    write_dominant_gtf(gtf, tx_ids, dom_gtf)

    feats['priority'] = feats['class_feature'].map(PRIORITY_MAP)
    feats['class'] = feats['class_feature'].map(CLASS_MAP)

    features_by_cell[ct] = {
        chrom: df.copy()
        for chrom, df in feats.groupby('chrom', observed=True)
    }

# =====================
# Chromosome-streamed binning + labeling
# =====================
output_path = os.path.join(
    OUTPUT_DIR, "multi_task_classification_targets.tsv"
)

write_header = True

for chrom, chrom_len in chrom_lengths.items():
    print(f"[INFO] Labeling chromosome {chrom}")

    bins = generate_bins_for_chrom(chrom, chrom_len, BIN_SIZE)
    bins['bin_id'] = np.arange(len(bins))

    chrom_labels = bins[['chrom', 'start', 'end']].copy()

    bin_gr = pr.PyRanges(
        bins.rename(columns={
            "chrom": "Chromosome",
            "start": "Start",
            "end": "End"
        })
    )

    for ct in CELL_TYPES:
        feats = features_by_cell.get(ct, {}).get(chrom)

        if feats is None:
            chrom_labels[f"{ct}_Class"] = CLASS_MAP['desert']
            continue

        feat_gr = pr.PyRanges(
            feats.rename(columns={
                "chrom": "Chromosome",
                "start": "Start",
                "end": "End"
            })
        )

        overlaps = (
            bin_gr.join(feat_gr, how="left")
            .as_df()
            .dropna(subset=["priority"])
        )

        if overlaps.empty:
            chrom_labels[f"{ct}_Class"] = CLASS_MAP['desert']
            continue

        maxp = overlaps.loc[
            overlaps.groupby('bin_id')['priority'].idxmax()
        ]

        label_map = maxp.set_index('bin_id')['class']

        chrom_labels[f"{ct}_Class"] = (
            chrom_labels.index.map(label_map)
            .fillna(CLASS_MAP['desert'])
            .astype(int)
        )

    chrom_labels.to_csv(
        output_path,
        sep="\t",
        index=False,
        header=write_header,
        mode="w" if write_header else "a"
    )

    write_header = False

print(f"[DONE] Multi-task target matrix written to {output_path}")
