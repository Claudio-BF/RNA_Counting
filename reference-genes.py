import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Configuration
DATA_PATH = "./data/all_sample_reads.csv"
CPM_BUDGET = 50000
MIN_TOTAL_READS = 100
N_RANGE = range(1, 200, 1)

ALPHA = 0.15


def load_and_preprocess(csv_path):
    if not os.path.exists(csv_path):
        print("CSV not found. Generating dummy data...")
        np.random.seed(42)
        genes = [f"Gene_{i}" for i in range(20000)]
        samples = [f"Sample_{i}" for i in range(10)]
        # Generate data with a realistic mean-variance relationship
        # High expression = low noise, Low expression = high noise
        means = np.exp(np.random.uniform(0, 10, 20000))
        data = np.zeros((20000, 10))
        for i, m in enumerate(means):
            # Noise decreases as mean increases
            r = m * 0.5
            p = r / (m + r)
            data[i] = np.random.negative_binomial(n=r, p=p, size=10)
        df = pd.DataFrame(data, index=genes, columns=samples)
    else:
        df = pd.read_csv(csv_path, index_col=0)

    counts = df.select_dtypes(include=[np.number])
    counts = counts[counts.sum(axis=1) >= MIN_TOTAL_READS]
    library_sizes = counts.sum(axis=0)

    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    log_cpm = np.log2(cpm + 1)

    stats = pd.DataFrame(
        {
            "mean_log": log_cpm.mean(axis=1),
            "mean_cpm": cpm.mean(axis=1),
            "stdev_log": log_cpm.std(axis=1),
        }
    )
    return counts, library_sizes, cpm, log_cpm, stats


def select_smart_weighted(stats, n_target, budget, alpha):
    """
    Selects n genes by minimizing a combined cost function.
    Cost = (Stability Z-Score) + alpha * (Distance from Target CPM Z-Score)
    """
    target_cpm = budget / n_target

    # 1. Metric: Distance from Target CPM (in Log space to handle scale)
    log_target = np.log2(target_cpm + 1)

    # Distance: How far is this gene's mean log-cpm from our target?
    dist = stats["mean_log"] - log_target
    dist_sq = np.where(dist > 0, dist**2, 0)

    # 2. Metric: Stability (Variance)
    stability_metric = stats["stdev_log"] ** 2

    # 3. Normalization
    norm_dist = dist_sq / np.median(dist_sq[dist_sq != 0])
    norm_stab = stability_metric / stability_metric.median()

    # 4. Total Cost
    # If a gene is very stable (low norm_stab), it can afford a higher norm_dist.
    total_cost = norm_stab + (alpha * norm_dist)

    # Select top n genes with lowest cost
    candidates = stats.assign(cost=total_cost).nsmallest(n_target, "cost")

    return candidates.index, candidates


def evaluate_quality(selected_genes, log_cpm):
    if len(selected_genes) == 0:
        return float("inf")
    return log_cpm.loc[selected_genes].mean(axis=0).std()


def calculate_extrapolation_error(
    selected_genes, selected_stats, counts, library_sizes
):
    log_raw = np.log(counts.loc[selected_genes] + 1)
    sample_gm_raw = np.exp(log_raw.mean(axis=0))
    ref_gm_cpm = np.exp(np.log(selected_stats["mean_cpm"] + 1).mean())
    extrapolated_totals = (sample_gm_raw / ref_gm_cpm) * 1_000_000
    diff_pct = (extrapolated_totals - library_sizes) / library_sizes * 100
    return diff_pct.abs().mean(), diff_pct, extrapolated_totals


# def calculate_extrapolation_error_linear(
#     selected_genes, selected_stats, counts, library_sizes
# ):
#     sums = counts.sum(axis=0)
#     total = np.mean(sums)
#     extrapolated_totals = (sample_gm_raw / ref_gm_cpm) * 1_000_000
#     return diff_pct.abs().mean(), diff_pct, extrapolated_totals


def main():
    try:
        counts, library_sizes, cpm, log_cpm, stats = load_and_preprocess(DATA_PATH)
    except Exception as e:
        print(e)
        return

    results = []
    print(f"Sweeping n ({N_RANGE.start}-{N_RANGE.stop})...")

    for n in N_RANGE:
        # Use the new Smart Weighted Selector
        selected_genes, selected_stats = select_smart_weighted(
            stats, n, CPM_BUDGET, alpha=ALPHA
        )

        quality = evaluate_quality(selected_genes, log_cpm)
        error, _, _ = calculate_extrapolation_error(
            selected_genes, selected_stats, counts, library_sizes
        )
        results.append((n, quality, error))

    res_df = pd.DataFrame(results, columns=["n", "quality", "error"])
    best_row = res_df.loc[res_df["quality"].idxmin()]
    best_n = int(best_row["n"])

    # Final Selection
    final_genes, final_stats = select_smart_weighted(
        stats, best_n, CPM_BUDGET, alpha=ALPHA
    )

    final_quality = evaluate_quality(final_genes, log_cpm)
    final_error, _, _ = calculate_extrapolation_error(
        final_genes, final_stats, counts, library_sizes
    )

    print(f"\n====== Final Selection Results (n={best_n}) ======")
    print(f"Quality Score (SD of LogCPM Means): {final_quality:.6f}")
    print(f"Average Extrapolation Error: {final_error:.4f}%")
    print(f"Total CPM Sum: {final_stats['mean_cpm'].sum():.2f} (Target: {CPM_BUDGET})")

    print("\nSelected Genes List:")
    print(f"{'Gene ID':<20} | {'Mean CPM':>12} | {'SD Log2 CPM':>12}")
    print("-" * 52)

    # Sort by CPM descending for cleaner reading
    sorted_stats = final_stats.sort_values(by="mean_cpm", ascending=False)
    for gene_id, row in sorted_stats.iterrows():
        print(f"{gene_id:<20} | {row['mean_cpm']:>12.2f} | {row['stdev_log']:>12.4f}")
    print("=" * 52)

    # --- Plotting Results ---

    # Graph 1: Optimization Curve
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(res_df["n"], res_df["quality"], "b-o", markersize=3, label="Quality (SD)")
    ax1.set_xlabel("n (Count)")
    ax1.set_ylabel("Quality Score", color="b")
    ax2 = ax1.twinx()
    ax2.plot(res_df["n"], res_df["error"], "r-x", markersize=3, label="Error %")
    ax2.set_ylabel("Extrapolation Error %", color="r")
    plt.title(f"Optimization (Budget: {CPM_BUDGET}, Alpha: {ALPHA})")
    plt.axvline(best_n, color="k", linestyle="--", alpha=0.5)
    plt.show()

    # Graph 2: The "Cluster" Shape
    plt.figure(figsize=(10, 6))

    # Background: All genes
    plt.scatter(
        stats["mean_log"],
        stats["stdev_log"],
        c="silver",
        s=5,
        alpha=0.2,
        label="All Genes",
    )

    target_log = np.log2((CPM_BUDGET / best_n) + 1)

    plt.scatter(
        final_stats["mean_log"],
        final_stats["stdev_log"],
        c=np.abs(final_stats["mean_log"] - target_log),  # Color by distance
        cmap="viridis_r",
        s=30,
        edgecolor="k",
        linewidth=0.5,
        label=f"Selected (n={best_n})",
    )

    plt.axvline(target_log, color="red", linestyle="--", alpha=0.8, label="Target CPM")
    plt.colorbar(label="Log-Distance from Target")
    plt.xlabel("Mean Log2-CPM")
    plt.ylabel("SD Log2-CPM (Stability)")
    plt.title(
        f"Smart Weighted Selection (n={best_n})\nAllows deviation from target if stability is high"
    )
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()


if __name__ == "__main__":
    main()
