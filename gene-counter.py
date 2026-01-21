import re
import csv
from itertools import islice
from itertools import accumulate
import pandas as pd

# Constants
GENE_LIST_PATH = "./data/pca-circadian-filtering-top-genes-ensdargID.txt"
DATA_PATH = "./data/all_sample_reads.csv"
MIN_READS = 100

# reads all the gene names in GENE_LIST_PATH, priority is the order that the first instance appears.
# Tries to add genes in order of priority, conditional on every gene in the list being predicted atleast MIN_READS reads.


def get_favorite_ids(filepath):
    with open(filepath, "r") as f:
        return dict.fromkeys(re.findall(r"ENSDARG\d{11}", f.read()))


def calculate_inclusion(gene_dict, grand_total):
    items = list(gene_dict.items())
    current_total = 0
    min_gene = items[0][0]
    minimum = items[0][1]
    genes, counts = [], []
    for gene_id, total in items:
        current_total += total
        new_min = min(minimum, total)
        if new_min * (grand_total / current_total) >= MIN_READS:
            genes.append(gene_id)
            counts.append(total)
            if total <= minimum:
                minimum = new_min
                min_gene = gene_id
        else:
            current_total -= total

    print(
        f"You should include the following {len(genes)} genes so that every gene is predicted to be read at least {MIN_READS} times:"
    )
    for i in range(len(genes)):
        print(f"{genes[i]:<19} {counts[i]:,}")
    print(
        f"When reading these genes {min_gene} has the minimal number of reads ({minimum}) and is projected to have {minimum * grand_total / current_total:.2f} reads when only reading these genes"
    )


def process_gene_data(gene_list, csv_path):
    results = {gene_id: 0 for gene_id in gene_list}
    found_genes = set()

    with open(csv_path, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        grand_total = 0
        for count, row in enumerate(reader):
            gene_id = row.get("geneID")
            total = sum(int(v) if v.isdigit() else 0 for k, v in row.items())
            grand_total += total
            if gene_id in results:
                results[gene_id] = total
                found_genes.add(gene_id)
    print(
        f"There are a total of {grand_total} reads across all {count + 1} genes in {DATA_PATH.split('/')[-1]}"
    )
    calculate_inclusion(results, grand_total)


if __name__ == "__main__":
    favorites = get_favorite_ids(GENE_LIST_PATH)
    print(
        f"Loaded {len(favorites)} unique favorite genes from {GENE_LIST_PATH.split('/')[-1]} (priority given by first seen order)."
    )
    process_gene_data(favorites, DATA_PATH)
