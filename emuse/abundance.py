import pandas as pd


def read_rel_abundance(path):
    abundance = pd.read_csv(path, sep="\t")
    # Filter for wanted columns
    columns_to_keep = ["abundance", "species", "genus", "family", "tax_id", "estimated counts"]
    abundance_filtered = abundance[columns_to_keep]
    # Rename col names
    abundance_filtered = abundance_filtered.rename(
        columns={
            "estimated counts": "estimated read counts",
            "tax_id": "tax id"
        }
    )
    # Sort based on descending abundance
    abundance_ordered = abundance_filtered.sort_values(by="abundance", ascending=False)
    # Re-index the table
    return abundance_ordered.reset_index(drop=True)
