# Imm_cell_rMATS

RNA-seq splicing data preprocessing pipeline for training splice-aware sequence models.
Produces genome-wide PSI-labeled NPZ shards from aligned BAM files across immune cell types (B cell, T cell, Macrophage).

## Pipeline Overview

```
RNA-seq BAM files
      │
      ├─────────────────────────────────────┐
      │  Branch A: Annotation               │  Branch B: Event Quantification
      ▼                                     ▼
 StringTie (per replicate)              rMATS (per cell type)
      │                                     │
      ▼                                     ▼
 Dominant isoform selection           Event tables (SE, A5SS, A3SS, RI, MXE)
      │                                     │
      ▼                                     ▼
 Transcripts TSV                      Per-sample PSI BED files (soft labels)
      │                                     │
      └──────────────────┬──────────────────┘
                         │
                         ▼
              Master BED aggregation (per cell type)
                         │
                         ▼
              Multitask NPZ shard conversion
              (sequence + PSI + mask, binned at 32 bp, 3 cell-type tasks)
                         │
                         ▼
              train / val / test NPZ shards
                         │
                         ▼
              Single-cell NPZ extraction
              (one task column per cell type)
```

## Repository Structure

```
├── config/                     # Input configuration
│   ├── pipeline_config.sh      # Central path/parameter config — edit this first
│   ├── B_Cell_bams.txt         # BAM file lists per cell type
│   ├── T_Cell_bams.txt
│   ├── Macrophage_bams.txt
│   ├── manifest.tsv            # Hand-curated per-sample BED manifest
│   └── example_manifest.tsv   # Template for manifest format
│
├── src/                        # Python source (all scripts importable from same dir)
│   ├── io_utils.py                   # Shared BED/TSV loaders
│   ├── gtf_to_tsv.py                 # Step 1 – GTF → transcript TSV
│   ├── generate_splicing_labels.py   # Core PSI labeling logic
│   ├── label_splicing_events.py      # Step 2 – rMATS PSI → per-sample BEDs + master BEDs
│   ├── add_confidence_to_bed.py      # Confidence score assignment
│   ├── emit_soft_labels_per_sample.py
│   ├── create_master_beds.py         # Aggregate replicates → master BED per cell type
│   ├── make_manifest.py              # Build/update BED manifest TSV
│   ├── dominant_isoform_labler.py    # Select dominant isoform from StringTie GTF
│   ├── data_conversion_fast.py       # Step 3 – BED + FASTA → NPZ shards
│   ├── filter_npz_by_celltype.py     # Step 4 – Extract single cell type from multitask NPZ
│   └── validation_tool.py
│
├── scripts/                    # Shell pipeline launchers
│   ├── run_full_pipeline.sh    # Top-level orchestrator (Steps 1–3)
│   ├── 01_rmats/               # Step 0 – rMATS splice event detection
│   │   ├── launch_rmats.sh         # Pairwise comparisons (B vs T, etc.)
│   │   └── launch_rmats_solo.sh    # Single cell-type mode (for PSI labeling)
│   ├── 02_annotation/          # Step A – Data-driven annotation
│   │   ├── launch_merge_stringtie.sh
│   │   ├── merge_cell_type_bams.sh
│   │   ├── stringtie_annotation.sh
│   │   ├── stringtie_annotation_parallel.sh
│   │   └── stringtie_per_replicate.sh
│   ├── 03_labeling/            # Step 2 – PSI labeling
│   │   ├── launch_label_splicing_events.sh  # Main labeling launcher
│   │   ├── launch_gtf_to_tsv.sh
│   │   ├── launch_emit_soft_per_sample.sh
│   │   ├── launch_add_confidence_to_bed.sh
│   │   ├── launch_manifest_maker.sh
│   │   ├── launch_data_conversion.sh
│   │   ├── launch_data_conversion_test.sh
│   │   └── launch_debug_data_conversion.sh
│   └── 04_dataset/             # Steps 3–4 – NPZ shard generation
│       ├── launch_convert_parallel.sh          # Step 3 – multitask NPZ
│       ├── launch_convert_parallel_aug4_eonly.sh
│       ├── launch_convert_parallel_aug8.sh
│       └── make_single_cell_type.sh            # Step 4 – single-cell NPZ extraction
│
├── analysis/                   # Visualization and QC
│   ├── check_bams.py
│   ├── extract_robust_hits.py
│   ├── plot_rmats_volcano.py
│   ├── plot_rmats_volcano_capped_y.py
│   ├── launch_plotter.sh
│   ├── launch_robust_hits.sh
│   ├── rmats_global_vis.ipynb
│   └── sashimi/
│
└── notebooks/                  # Exploratory analysis
    ├── inspect_labels.ipynb
    └── inspect_labels_plot.ipynb
```

## Quick Start

**1. Configure paths** in [config/pipeline_config.sh](config/pipeline_config.sh) — set `PROJECT_ROOT`, genome FASTA, and output directories.

**2. Run rMATS** (per cell type, for PSI labeling):
```bash
bash scripts/01_rmats/launch_rmats_solo.sh
```

**3. Build annotation** (StringTie per replicate → dominant isoform GTF):
```bash
bash scripts/02_annotation/launch_merge_stringtie.sh
```

**4. Run the labeling + dataset pipeline**:
```bash
bash scripts/run_full_pipeline.sh
# or step-by-step:
bash scripts/run_full_pipeline.sh --only 1   # GTF → TSV
bash scripts/run_full_pipeline.sh --only 2   # rMATS PSI → BEDs (soft labels)
bash scripts/run_full_pipeline.sh --only 3   # BEDs + FASTA → multitask NPZ shards
```

**5. Extract single-cell-type NPZ datasets** (run after Step 3):
```bash
bash scripts/04_dataset/make_single_cell_type.sh
```

## Output Format

Each NPZ shard contains:
- `sequence` — one-hot encoded DNA, shape `(N, window_bp, 4)`
- `targets` — PSI values per bin per cell type, shape `(N, n_bins, n_tasks)`
- `confidence` — confidence weights, shape `(N, n_bins, n_tasks)`
- `mask` — coverage mask (1 = has signal), shape `(N, n_bins, n_tasks)`
- `chrom`, `window_start`, `window_end` — genomic coordinates
- `has_signal` — float flag, 1.0 if window overlaps any event

Default: 524,288 bp windows binned at 32 bp (16,384 bins), 3 cell-type tasks.
