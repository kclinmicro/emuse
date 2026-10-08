from importlib import resources

import pandas as pd


class TaxTranslator(object):
    def __init__(self, taxonomy_path=None):
        if taxonomy_path is None:
            taxonomy_path = resources.files("emuse") / "data" / "taxonomy.tsv"
        self._taxdf = pd.read_csv(taxonomy_path, sep="\t", dtype=str).set_index("tax_id")
        self._taxid_to_label_mapping = {
            tax_id: self.get_best_tax_label(row) for tax_id, row in self._taxdf.iterrows()
        }

    def taxid_to_label(self, taxid):
        return self._taxid_to_label_mapping.get(taxid, taxid)

    def translate_taxids_in_df_columns(self, df):
        df_cols_orig = df.columns.tolist()
        new_headers = [
            self.taxid_to_label(col.strip()) if col.strip().isdigit() else col
            for col in df_cols_orig
        ]
        df.columns = new_headers
        return df

    def get_best_tax_label(self, row):
        """Return the best available taxonomic label from left to right."""
        for level in [
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "clade",
            "superkingdom",
            "subspecies",
            "species subgroup",
            "species group",
        ]:
            val = row.get(level, "")
            if pd.notna(val) and str(val).strip() != "":
                return val.strip()
        return "Unknown"
