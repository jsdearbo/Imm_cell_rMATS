"""
Filter NPZ dataset shards to extract single cell type.

This script creates new filtered NPZ shards containing only a specific cell type,
reducing the targets tensor from (N, T, L) to (N, 1, L) where T is the nth cell type.

The cell-type-to-index mapping is auto-detected from the ``celltype_idx_*``
metadata keys that ``data_conversion_fast.py`` writes into every shard.
If the metadata is missing (older datasets), pass ``--cell-type-order`` manually.

Usage:
    python filter_npz_by_celltype.py \
        --input-dir /path/to/dataset \
        --output-dir /path/to/filtered_dataset \
        --cell-type "B_Cell"
"""

import os
import sys
import argparse
import numpy as np
import glob
from pathlib import Path


def detect_celltype_mapping(shard_path):
    """
    Read the celltype_idx_* metadata keys from a single NPZ shard.

    Returns:
        dict  {cell_type_name: int_index}  e.g. {"B_Cell": 0, "Macrophage": 1, "T_Cell": 2}
        None  if no metadata keys are found
    """
    ct_map = {}
    with np.load(shard_path, allow_pickle=False) as data:
        for key in data.files:
            if key.startswith("celltype_idx_"):
                idx = int(key.split("_")[-1])
                ct_map[str(data[key])] = idx
    return ct_map if ct_map else None


def get_celltype_index(cell_type, celltype_order):
    """
    Get the task index for a specific cell type from a comma-separated order string.

    Args:
        cell_type: Name of cell type (e.g., "B_Cell")
        celltype_order: Comma-separated list of cell types in order

    Returns:
        Index of the cell type in the ordering
    """
    celltypes = [ct.strip() for ct in celltype_order.split(",")]
    if cell_type not in celltypes:
        raise ValueError(f"Cell type '{cell_type}' not found in order: {celltypes}")
    return celltypes.index(cell_type)


def filter_npz_shard(input_path, output_path, cell_type_idx):
    """
    Filter a single NPZ shard to extract one cell type.
    
    Args:
        input_path: Path to input NPZ file
        output_path: Path to output NPZ file
        cell_type_idx: Index of the cell type to extract (0-based)
    """
    print(f"Processing: {os.path.basename(input_path)}")
    
    try:
        with np.load(input_path, mmap_mode="r") as data:
            # Load all arrays
            sequence = data["sequence"][:]  # (N, L, 4) or (N, 4, L)
            targets = data["targets"][:]     # (N, T, L) or similar
            mask = data["mask"][:] if "mask" in data else None
            confidence = data["confidence"][:] if "confidence" in data else None
            
            # Extract single cell type
            # targets shape: (N, L, T) where T is number of tasks
            if targets.ndim == 3:
                # (N, L, T) -> extract T [cell_type_idx] -> (N, L, 1)
                targets_filtered = targets[:, :, cell_type_idx:cell_type_idx+1]
            elif targets.ndim == 2:
                # (N, L) - only one task already
                targets_filtered = targets
            else:
                raise ValueError(f"Unexpected targets shape: {targets.shape}")
            
            # Filter mask and confidence similarly
            mask_filtered = None
            if mask is not None:
                if mask.ndim == 3:
                    mask_filtered = mask[:, :, cell_type_idx:cell_type_idx+1]
                elif mask.ndim == 2:
                    mask_filtered = mask
                else:
                    raise ValueError(f"Unexpected mask shape: {mask.shape}")
            
            confidence_filtered = None
            if confidence is not None:
                if confidence.ndim == 3:
                    confidence_filtered = confidence[:, :, cell_type_idx:cell_type_idx+1]
                elif confidence.ndim == 2:
                    confidence_filtered = confidence
                else:
                    raise ValueError(f"Unexpected confidence shape: {confidence.shape}")
            
            # Save filtered NPZ
            save_dict = {
                "sequence": sequence,
                "targets": targets_filtered,
            }
            if mask_filtered is not None:
                save_dict["mask"] = mask_filtered
            if confidence_filtered is not None:
                save_dict["confidence"] = confidence_filtered
            
            # Copy any other arrays (metadata, etc)
            for key in data.files:
                if key not in ["sequence", "targets", "mask", "confidence"]:
                    arr = data[key]
                    # Handle 0-dimensional arrays (scalars) without slicing
                    if arr.ndim == 0:
                        save_dict[key] = arr.item()
                    else:
                        save_dict[key] = arr[:]
            
            np.savez_compressed(output_path, **save_dict)
            
        print(f"  ✓ Saved to: {os.path.basename(output_path)}")
        return True
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


def filter_dataset_split(input_dir, output_dir, split, cell_type_idx):
    """
    Filter all NPZ shards in a split (train, val, test).
    
    Args:
        input_dir: Base input directory
        output_dir: Base output directory
        split: Split name (train, val, test)
        cell_type_idx: Index of cell type to extract
    """
    split_input_dir = os.path.join(input_dir, split)
    split_output_dir = os.path.join(output_dir, split)
    
    if not os.path.exists(split_input_dir):
        print(f"Split '{split}' not found at {split_input_dir}")
        return False
    
    # Create output directory
    os.makedirs(split_output_dir, exist_ok=True)
    
    # Find all NPZ shards
    shard_paths = sorted(glob.glob(os.path.join(split_input_dir, "*.npz")))
    if not shard_paths:
        print(f"No NPZ shards found in {split_input_dir}")
        return False
    
    print(f"\nFiltering '{split}' split ({len(shard_paths)} shards)...")
    success_count = 0
    
    for shard_path in shard_paths:
        shard_name = os.path.basename(shard_path)
        output_path = os.path.join(split_output_dir, shard_name)
        
        if filter_npz_shard(shard_path, output_path, cell_type_idx):
            success_count += 1
    
    print(f"✓ Completed {split}: {success_count}/{len(shard_paths)} shards processed")
    return success_count == len(shard_paths)


def main():
    parser = argparse.ArgumentParser(
        description="Filter NPZ dataset shards to extract a single cell type"
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Base directory containing train/val/test splits"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory where filtered train/val/test will be created"
    )
    parser.add_argument(
        "--cell-type",
        required=True,
        help="Cell type to extract (e.g., 'B_Cell')"
    )
    parser.add_argument(
        "--cell-type-order",
        default=None,
        help="Comma-separated cell types in task-index order. "
             "Auto-detected from NPZ metadata when omitted."
    )
    parser.add_argument(
        "--splits",
        default="train,val,test",
        help="Comma-separated splits to process (default: train,val,test)"
    )
    
    args = parser.parse_args()
    
    # Validate input
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory not found: {args.input_dir}")
        sys.exit(1)

    # ---- resolve cell-type-to-index mapping ----
    if args.cell_type_order is not None:
        # Explicit override
        celltypes = [ct.strip() for ct in args.cell_type_order.split(",")]
        try:
            cell_type_idx = get_celltype_index(args.cell_type, args.cell_type_order)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        # Auto-detect from the first shard we can find
        splits = [s.strip() for s in args.splits.split(",")]
        ct_map = None
        for split in splits:
            split_dir = os.path.join(args.input_dir, split)
            shards = sorted(glob.glob(os.path.join(split_dir, "*.npz")))
            if shards:
                ct_map = detect_celltype_mapping(shards[0])
                if ct_map is not None:
                    break

        if ct_map is None:
            print("Error: Could not auto-detect cell type mapping from NPZ metadata.")
            print("       Pass --cell-type-order explicitly (e.g. 'B_Cell,Macrophage,T_Cell').")
            sys.exit(1)

        if args.cell_type not in ct_map:
            print(f"Error: Cell type '{args.cell_type}' not found in NPZ metadata: {ct_map}")
            sys.exit(1)

        cell_type_idx = ct_map[args.cell_type]
        # Reconstruct ordered list for display
        celltypes = [ct for ct, _ in sorted(ct_map.items(), key=lambda x: x[1])]
    print(f"\n{'='*70}")
    print(f"FILTERING DATASET TO SINGLE CELL TYPE")
    print(f"{'='*70}")
    print(f"Input:      {args.input_dir}")
    print(f"Output:     {args.output_dir}")
    print(f"Cell type:  {args.cell_type} (index {cell_type_idx} of {celltypes})")
    print(f"Splits:     {args.splits}")
    print(f"{'='*70}\n")
    
    # Process each split
    splits = [s.strip() for s in args.splits.split(",")]
    all_success = True
    
    for split in splits:
        success = filter_dataset_split(
            args.input_dir,
            args.output_dir,
            split,
            cell_type_idx
        )
        all_success = all_success and success
    
    if all_success:
        print(f"\n✓ Filtering complete!")
        print(f"Filtered dataset ready at: {args.output_dir}")
        sys.exit(0)
    else:
        print(f"\n✗ Some splits failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
