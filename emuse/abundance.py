import pandas as pd


def read_rel_abundance(path):
    abundance = pd.read_csv(path, sep="\t")
    # Filter for wanted columns
    abundance_filtered = abundance.iloc[:, list(range(5)) + [13]]
    # Move the first column (taxid)
    abundance_switched = abundance_filtered[abundance_filtered.columns[1:5]
        .append(abundance_filtered.columns[:1])
        .append(abundance_filtered.columns[5:])]
    # Rename col names
    abundance_switched = abundance_switched.rename(
        columns={
            "estimated counts": "estimated read counts",
            "tax_id": "tax id"
        }
    )
    # Sort based on descending abundance
    abundance_ordered = abundance_switched.sort_values(by="abundance", ascending=False)
    # Re-index the table
    return abundance_ordered.reset_index(drop=True)
