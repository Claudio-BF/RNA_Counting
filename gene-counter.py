import re
import csv
from itertools import islice
from itertools import accumulate

# Constants
GENE_LIST_PATH = "./data/pca-circadian-filtering-top-genes-ensdargID.txt"
DATA_PATH = "./data/all_sample_reads.csv"
MIN_READS = 10

# reads all the gene names in GENE_LIST_PATH, priority is the order that the first instance appears.
# Will remove genes in reverse order of priority until the gene with the lowest number of reads in the current list is projected to be read at least MIN_READS times.
# After that, adds as many genes as possible which which have more reads than the gene in list with the minimum number of reads.


def get_favorite_ids(filepath):
    with open(filepath, "r") as f:
        return dict.fromkeys(re.findall(r"ENSDARG\d{11}", f.read()))


def try_int(value):
    try:
        return int(value)
    except (ValueError, TypeError):
        return 0


def calulate_inclusion(gene_dict):
    items = list(gene_dict.items())
    keys, values = zip(*items)
    grand_total = sum(values)

    # pre-calculate stats in optimal time complexity (lazy premade library solution)
    totals = list(accumulate(values))
    mins = list(accumulate(values, func=min))

    # iterate backwards through pre-calculated stats
    for i in range(len(items) - 1, -1, -1):
        minimum = mins[i]
        current_total = totals[i]
        if minimum * (grand_total / current_total) >= 10:
            final_idx = i
            min_gene = keys[mins.index(minimum)]

            print(f"You can include up to the first {final_idx + 1} genes.")
            print(
                f"Gene {min_gene} has the least reads among these with {minimum} reads."
            )
            print(
                f"{min_gene} is projected to have {minimum * grand_total / current_total:.2f} reads when reading only these genes."
            )
            # figure out the other genes you can include
            max_total = minimum * grand_total / 10.0
            other_genes = []
            for j in range(i + 1, len(values)):
                if values[j] >= minimum:
                    current_total += values[j]
                    if current_total >= max_total:
                        current_total -= values[j]
                        break
                    else:
                        other_genes.append(keys[j])
            print(
                f"You might also include the following genes, at which point {min_gene} is projected to have {minimum * grand_total / current_total:.2f} reads"
            )
            print(*other_genes, sep="\n")
            break


def process_gene_data(gene_list, csv_path):
    results = {gene_id: 0 for gene_id in gene_list}
    found_genes = set()

    with open(csv_path, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene_id = row.get("geneID")
            if gene_id in results:
                total = sum(try_int(v) for k, v in row.items() if k != "geneID")
                results[gene_id] = total
                found_genes.add(gene_id)

    for gene_id, total in results.items():
        if gene_id in found_genes:
            print(f"{gene_id} {total}")
    calulate_inclusion(results)


if __name__ == "__main__":
    favorites = get_favorite_ids(GENE_LIST_PATH)
    print(f"Loaded {len(favorites)} unique favorite genes.")
    process_gene_data(favorites, DATA_PATH)
