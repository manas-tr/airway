from pathlib import Path
import pandas as pd

COUNTS_URL = "https://raw.githubusercontent.com/bioconnector/workshops/refs/heads/master/data/airway_rawcounts.csv"

META_URL = "https://raw.githubusercontent.com/bioconnector/workshops/refs/heads/master/data/airway_metadata.csv"

data_dir = Path("data")
data_dir.mkdir(exist_ok=True)

counts = pd.read_csv(COUNTS_URL)

meta = pd.read_csv(META_URL)

counts.to_csv(data_dir / "airway_rawcounts.csv", index=False)
meta.to_csv(data_dir / "airway_metadata.csv", index=False)


print("saved!")

print("\ncounts shape:")
print(counts.shape)

print("\nmetadata shape:")
print(meta.shape)

print("\nmetadata preview:")
print(meta.head())