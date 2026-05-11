from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

DATA_DIR = Path("data")
OUTPUT_DIR = Path("outputs")
FIGURE_DIR = Path("figures")


counts = pd.read_csv(DATA_DIR / "airway_rawcounts.csv")
metadata = pd.read_csv(DATA_DIR / "airway_metadata.csv")

print("Counts shape:", counts.shape)
print("Metadata shape:", metadata.shape)

#counts matrix

gene_ids = counts.iloc[:, 0]

count_matrix = counts.iloc[:, 1:]

count_matrix.index = gene_ids

#rows=samples
#columns=genes
count_matrix = count_matrix.T

print("transposed count matrix shape:", count_matrix.shape)

#metadata

metadata = metadata.set_index("id")

metadata = metadata[["dex"]]

print(metadata.head())


count_matrix = count_matrix.loc[metadata.index]

print("\norder verified.")


print("\nfiltering low-count genes...")

keep = count_matrix.sum(axis=0) >= 10

count_matrix = count_matrix.loc[:, keep]

print("remaining genes:", count_matrix.shape[1])

#differential expression

dds = DeseqDataSet(
    counts=count_matrix,
    metadata=metadata,
    design="~ dex",
    refit_cooks=True
)

dds.deseq2()

stats = DeseqStats(
    dds,
    contrast=("dex", "treated", "control")
)

stats.summary()

results = stats.results_df

print("\nresults")
print(results.head())

results.to_csv(
    OUTPUT_DIR / "deseq_results.csv"
)

print("\nsaved differential expression results.")

normalized = np.log1p(count_matrix)

pca = PCA(n_components=2)

pca_result = pca.fit_transform(normalized)

pca_df = pd.DataFrame({
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
    "dex": metadata["dex"].values
})

plt.figure(figsize=(7, 5))

for group in pca_df["dex"].unique():
    subset = pca_df[pca_df["dex"] == group]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of airway RNA-seq samples")
plt.legend()

plt.tight_layout()

plt.savefig(
    FIGURE_DIR / "pca.png",
    dpi=300
)

plt.close()

print("saved PCA plot.")


results = results.dropna()

results["neglog10_padj"] = -np.log10(results["padj"])

significant = (
    (results["padj"] < 0.05)
    &
    (abs(results["log2FoldChange"]) > 1)
)

plt.figure(figsize=(8, 6))

plt.scatter(
    results["log2FoldChange"],
    results["neglog10_padj"],
    s=10,
    alpha=0.7
)

plt.scatter(
    results.loc[significant, "log2FoldChange"],
    results.loc[significant, "neglog10_padj"],
    s=12
)

plt.axhline(
    -np.log10(0.05),
    linestyle="--"
)

plt.axvline(
    1,
    linestyle="--"
)

plt.axvline(
    -1,
    linestyle="--"
)

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 adjusted p-value")
plt.title("Differential Expression Volcano Plot")

plt.tight_layout()

plt.savefig(
    FIGURE_DIR / "volcano.png",
    dpi=300
)

plt.close()

print("saved volcano plot.")


sig_genes = significant.sum()

summary_text = f"""
Airway RNA-seq analysis summary

Samples analyzed: {count_matrix.shape[0]}
Genes analyzed: {count_matrix.shape[1]}

Significant genes:
padj < 0.05 and |log2FC| > 1

Total significant genes: {sig_genes}
"""

with open(
    OUTPUT_DIR / "summary.txt",
    "w"
) as f:
    f.write(summary_text)

print("\nsaved summary.")
