import csv
from collections import Counter

DATA_PATH = "./data/all_sample_reads.csv"
NUM_GENES = 50


def get_top_genes(csv_path):
    gene_counts = Counter()
    with open(csv_path, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene_id = row.get("geneID")
            total = sum(int(v) if v.isdigit() else 0 for k, v in row.items())
            gene_counts[gene_id] = total
    return gene_counts.most_common(NUM_GENES)


if __name__ == "__main__":
    print(f"Calculating top {NUM_GENES} genes...")
    top = get_top_genes(DATA_PATH)

    print(f"{'Rank':<5} {'Gene ID':<20} {'Total Reads'}")
    print("-" * 40)
    for rank, (gene_id, count) in enumerate(top, 1):
        print(f"{rank:<5} {gene_id:<20} {count:,}")
