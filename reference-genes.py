import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

# Configuration
DATA_PATH = "./data/all_sample_reads.csv"
NUM_GENES = 50
MIN_TOTAL_READS = 100


def get_stable_genes(csv_path):
    df = pd.read_csv(csv_path, index_col=0)
    counts = df.select_dtypes(include=[np.number])
    counts = counts[counts.sum(axis=1) >= MIN_TOTAL_READS]
    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    log_cpm = np.log2(cpm + 1)
    results = pd.DataFrame(
        {
            "stdev_log": log_cpm.std(axis=1),
            "mean_cpm": cpm.mean(axis=1),
            "total_reads": counts.sum(axis=1),
        }
    )

    # --- NEW SECTION STARt ---
    plt.figure(figsize=(10, 6))
    plt.scatter(
        results["total_reads"],
        results["stdev_log"],
        alpha=0.5,
        s=10,
    )
    plt.xlim(0, 10000)
    plt.xlabel("Number of Reads")
    plt.ylabel("Log-SD")
    plt.title("Gene Stability vs. Read Abundance")
    plt.grid(True, alpha=0.3)
    plt.show()
    # --- NEW SECTION END ---

    # Sort by stability (lowest SD) and return top N
    return results.sort_values("stdev_log").head(NUM_GENES)

    # Sort by stability (lowest SD) and return top N
    return results.sort_values("stdev_log").head(NUM_GENES)


if __name__ == "__main__":
    print(f"Processing {DATA_PATH}...")
    print("Calculating stability using Log-Variance heuristic (Vectorized)...")

    top_genes = get_stable_genes(DATA_PATH)

    print(f"\nTop {NUM_GENES} Reference Genes (Lowest SD of Log2-CPM)")
    print(
        f"{'Rank':<5} {'Gene ID':<20} {'Log-SD':<10} {'Mean CPM':<12} {'Total Reads'}"
    )
    print("-" * 65)

    # Reset index so we can access Gene ID as a column for printing
    for i, (gene_id, row) in enumerate(top_genes.iterrows(), 1):
        print(
            f"{i:<5} {str(gene_id):<20} "
            f"{row['stdev_log']:<10.3f} "
            f"{row['mean_cpm']:<12.1f} "
            f"{int(row['total_reads']):,}"
        )
