import numpy as np
import matplotlib.pyplot as plt
import glob
import os

def load_shard(shard_path):
    """Loads a single NPZ shard and returns the data dict."""
    try:
        data = np.load(shard_path, allow_pickle=True)
        return {k: data[k] for k in data.files}
    except Exception as e:
        print(f"Error loading {shard_path}: {e}")
        return None

def decode_seq(one_hot_seq):
    """Decodes a (L, 4) one-hot sequence back to string ACGT."""
    # map 0->A, 1->C, 2->G, 3->T based on script logic
    idx_to_base = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    
    # Argmax over the channel dimension (last dim)
    # Check if row is all zeros (N)
    is_N = one_hot_seq.sum(axis=1) == 0
    indices = one_hot_seq.argmax(axis=1)
    
    seq_chars = []
    for i, idx in enumerate(indices):
        if is_N[i]:
            seq_chars.append('N')
        else:
            seq_chars.append(idx_to_base[idx])
    return "".join(seq_chars)

def inspect_shard(shard_path, example_idx=0):
    """
    Visualizes one example from a shard.
    Checks for the 'Sparse Multi-task' issue.
    """
    data = load_shard(shard_path)
    if data is None: return

    print(f"--- Inspecting {os.path.basename(shard_path)} ---")
    print(f"Total examples in shard: {data['targets'].shape[0]}")
    
    # Extract specific example
    targets = data['targets'][example_idx] # Shape: (bins, tasks)
    seq = data['sequence'][example_idx]    # Shape: (seq_len, 4)
    cell_type = str(data['cell_type'][example_idx])
    chrom = str(data['chrom'][example_idx])
    start = data['window_start'][example_idx]
    end = data['window_end'][example_idx]
    
    print(f"\nExample {example_idx}:")
    print(f"Location: {chrom}:{start}-{end}")
    print(f"Labeled Cell Type: {cell_type}")
    print(f"Target Shape: {targets.shape}")
    
    # 1. Check Sparsity (The Major Flaw Check)
    # Sum across bins for each task
    task_activity = targets.sum(axis=0)
    active_tasks = np.where(task_activity > 0)[0]
    
    print(f"\n[Gut Check] Active Task Indices (Non-zero columns): {active_tasks}")
    
    if len(active_tasks) == 1:
        print("⚠️  WARNING: Only 1 task column is active. This confirms 'Sparse Multi-task' structure.")
        print("    If your model expects dense multi-task labels, this data is problematic.")
    elif len(active_tasks) == 0:
        print("⚠️  WARNING: No active targets found (all zeros). Signal might be missing for this window.")
    else:
        print("✅  Good: Multiple tasks have data (Dense Multi-task).")

    # 2. Visualize Targets
    plt.figure(figsize=(15, 6))
    
    # Plot heatmap of targets (transposed so tasks are rows)
    # Subsetting to center if needed, but showing full here
    plt.subplot(2, 1, 1)
    import seaborn as sns
    sns.heatmap(targets.T, cmap="viridis", cbar_kws={'label': 'PSI'})
    plt.title(f"Target Heatmap (Tasks x Bins) - {cell_type}")
    plt.xlabel("Bins")
    plt.ylabel("Task Index")

    # 3. Plot the specific active track line
    plt.subplot(2, 1, 2)
    # If specific cell type is known, find its index? 
    # The script doesn't save the task map in the NPZ, only the cell_type name.
    # We will plot whatever is non-zero.
    
    found_signal = False
    for t_idx in active_tasks:
        plt.plot(targets[:, t_idx], label=f"Task {t_idx}", alpha=0.7)
        found_signal = True
        
    if not found_signal:
        plt.plot(targets[:, 0], label="Task 0 (All Zeros)")

    plt.title(f"PSI Values across Window ({cell_type} active)")
    plt.xlabel("Bins")
    plt.ylabel("PSI")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # 4. Check Sequence
    decoded = decode_seq(seq)
    print(f"\nSequence Snippet (Center): ...{decoded[len(decoded)//2 - 50 : len(decoded)//2 + 50]}...")
    print(f"Sequence N content: {decoded.count('N') / len(decoded):.2%}")

# Usage Example for your Notebook:
# file_path = "./results/labelled_data/windowed_multitask/train/shard_0000.npz"
# inspect_shard(file_path, example_idx=0)