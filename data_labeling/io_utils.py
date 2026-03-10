"""
io_utils.py - Shared IO utilities for the splicing label pipeline.

Provides canonical loaders for BED label files and transcript TSVs,
ensuring a single source of truth for column names and parsing logic.

Used by:
    - test_soft_labeling.py   (soft label generation)
    - data_conversion_fast.py (NPZ shard creation)
"""

import logging
import pandas as pd

logger = logging.getLogger(__name__)

# Canonical BED column names used throughout the pipeline.
BED_COLUMNS = ["Chromosome", "Start", "End", "PSI", "confidence"]


def load_label_intervals(bed_path: str) -> pd.DataFrame:
    """Load a headerless BED file with PSI and confidence values.

    Expected format (tab-separated, no header):
        Chromosome  Start  End  PSI  confidence

    Returns:
        DataFrame with columns [Chromosome, Start, End, PSI, confidence].
        Empty DataFrame on error.
    """
    try:
        df = pd.read_csv(
            bed_path,
            sep="\t",
            names=BED_COLUMNS,
            dtype={"Chromosome": str},
        )
    except Exception as e:
        logger.error("Error reading %s: %s", bed_path, e)
        return pd.DataFrame(columns=BED_COLUMNS)

    # Handle files with fewer than 5 columns (confidence will be NaN)
    if df["confidence"].isna().all():
        df["confidence"] = 0.0
    df["confidence"] = df["confidence"].fillna(0.0)

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["PSI"] = df["PSI"].astype(float)
    df["confidence"] = df["confidence"].astype(float)
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    return df


def load_transcripts_tsv(path: str) -> pd.DataFrame:
    """Load a transcript coordinates TSV (with header row).

    Returns:
        DataFrame with at minimum [Chromosome, Start, End] columns.
    """
    df = pd.read_csv(path, sep="\t", dtype={"Chromosome": str})

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")

    original_len = len(df)
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)
    if len(df) < original_len:
        logger.warning(
            "Dropped %d transcript rows with invalid coordinates.",
            original_len - len(df),
        )

    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    return df
