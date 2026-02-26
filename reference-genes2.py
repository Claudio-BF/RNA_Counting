import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Configuration
DATA_PATH = "./data/all_sample_reads.csv"
BUDGET_RANGE = [
    1000,
    1500,
    2000,
    2500,
    3500,
    5000,
    7500,
    10000,
    17000,
    25000,
    30000,
    35000,
    40000,
    50000,
    75000,
    100000,
    120000,
    140000,
    160000,
    180000,
    200000,
]
MIN_TOTAL_READS = 100
N_RANGE = range(5, 400, 1)
ALPHA = 0.15


def load_and_preprocess(csv_path):
    if not os.path.exists(csv_path):
        print("CSV not found. Generating dummy data...")
        np.random.seed(42)
        genes = [f"Gene_{i}" for i in range(20000)]
        samples = [f"Sample_{i}" for i in range(10)]
        means = np.exp(np.random.uniform(0, 10, 20000))
        data = np.zeros((20000, 10))
        for i, m in enumerate(means):
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
            "stdev_log": log_cpm.std(axis=1),
        }
    )
    return counts, library_sizes, stats


def select_smart_weighted(stats, n_target, budget, alpha):
    target_cpm = budget / n_target
    log_target = np.log2(target_cpm + 1)

    dist = stats["mean_log"] - log_target
    dist_sq = np.where(dist > 0, dist**2, 0)
    stability_metric = stats["stdev_log"] ** 2

    norm_dist = dist_sq / np.median(dist_sq[dist_sq != 0])
    norm_stab = stability_metric / stability_metric.median()

    total_cost = norm_stab + (alpha * norm_dist)
    candidates = stats.assign(cost=total_cost).nsmallest(n_target, "cost")
    return candidates.index


def calculate_extrapolation_error_linear(selected_genes, counts, library_sizes):
    aggregate_counts = counts.loc[selected_genes].sum(axis=0)
    avg_aggregate_cpm = (aggregate_counts / library_sizes).mean() * 1_000_000
    extrapolated_totals = (aggregate_counts / avg_aggregate_cpm) * 1_000_000
    diff_pct = (extrapolated_totals - library_sizes) / library_sizes * 100
    return diff_pct.abs().mean()


def main():
    try:
        counts, library_sizes, stats = load_and_preprocess(DATA_PATH)
    except Exception as e:
        print(e)
        return

    best_results = []
    print(f"Sweeping Budgets: {BUDGET_RANGE}")

    for budget in BUDGET_RANGE:
        errors = []
        for n in N_RANGE:
            selected_genes = select_smart_weighted(stats, n, budget, alpha=ALPHA)
            err = calculate_extrapolation_error_linear(
                selected_genes, counts, library_sizes
            )
            errors.append(err)

        # Find the minimum error achieved for this budget
        min_err = min(errors)
        best_n = N_RANGE[errors.index(min_err)]
        best_results.append((budget, min_err, best_n))
        print(f"Budget {budget}: Best Error {min_err:.4f}% (at n={best_n})")

    # Plotting
    res_df = pd.DataFrame(best_results, columns=["Budget", "Error", "Best_N"])

    plt.figure(figsize=(10, 6))
    plt.plot(res_df["Budget"], res_df["Error"], marker="o", linestyle="-", color="b")

    # Annotate the Best N for each point
    for _, row in res_df.iterrows():
        plt.annotate(
            f"n={int(row['Best_N'])}",
            (row["Budget"], row["Error"]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
        )

    plt.xlabel("CPM Budget")
    plt.ylabel("Lowest Extrapolation Error (%)")
    plt.title("Best Normalization Performance vs. CPM Budget")
    plt.grid(True, alpha=0.5)
    plt.xscale("log")  # Log scale often helps visualize wide budget ranges
    plt.show()


if __name__ == "__main__":
    main()
