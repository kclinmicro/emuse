import argparse
from datetime import date
from importlib import resources
from importlib.metadata import version
from jinja2 import Environment, FileSystemLoader
import json
import yaml
from pathlib import Path
import pandas as pd
import re
import tomllib

from emuse.alignment import get_alignment_metrics

# Bundled package data (templates, CSS, taxonomy mapping, default config)
DATA_DIR = resources.files("emuse") / "data"
DEFAULT_CONFIG = DATA_DIR / "configs" / "config.toml"

def main():
    argp = argparse.ArgumentParser()
    argp.add_argument("-v", "--version", action="version", version=f"%(prog)s {version('emuse')}")
    argp.add_argument("-i", "--input-dir", type=str, required=True, help="Path to the input directory containing results")
    argp.add_argument("-o", "--output-file", type=str, required=True, help="Path to the output report file")
    argp.add_argument("-s", "--sample-name", type=str, required=True, help="Name of the sample")
    argp.add_argument("-n", "--neg-control", type=str, required=True, help="Name of the negative control")
    argp.add_argument("-c", "--config", type=Path, default=DEFAULT_CONFIG, help="Path to config file. Default is the bundled config.toml")
    argp.add_argument("-p", "--prob-score", action="store_true", help="Include probability score in the report")
    argp.add_argument("-m", "--alignment-metrics", action="store_true", help="Include metrics based on the raw alignment of reads to the database (percent identity and percent coverage)")

    args = argp.parse_args()

    # Read CSS file content
    with open(DATA_DIR / "static" / "style.css", "r") as f:
        css_content = f.read()

    env = Environment(loader=FileSystemLoader(str(DATA_DIR / "templates")))
    template = env.get_template("report.html.j2")

    # Set low abundance cutoff value
    LOW_ABUNDANCE_CUTOFF = 0.005

    # Load sample read assignment table
    assignment = pd.read_csv(f"{args.input_dir}/results/{args.sample_name}_downsampled.fastq_read-assignment-distributions.tsv", sep="\t")
    # Select all columns except the first one
    assignment_filtered = assignment.iloc[:, 1:]
    # Compute mean and median for each column
    assignment_summary = assignment_filtered.agg(['median', 'mean']).T.reset_index()
    # Rename columns
    assignment_summary.columns = ['tax id', 'median probability*', 'mean probability*']

    # Load neg control abundance table
    neg_control_abundance = pd.read_csv(f"{args.input_dir}/results/{args.neg_control}_downsampled.fastq_rel-abundance.tsv", sep="\t")
    # Filter for wanted columns
    neg_control_filtered = neg_control_abundance.iloc[:, list(range(5)) + [13]]
    # Move the first column (taxid)
    neg_control_switched = neg_control_filtered[neg_control_filtered.columns[1:5]
        .append(neg_control_filtered.columns[:1])
        .append(neg_control_filtered.columns[5:])]
    # Rename col names
    neg_control_switched = neg_control_switched.rename(
        columns={
            "estimated counts": "estimated read counts",
            "tax_id": "tax id"
        }
    )
    # Sort based on descending abundance
    neg_control_ordered = neg_control_switched.sort_values(by="abundance", ascending=False)
    # Re-index the table
    neg_control_ordered = neg_control_ordered.reset_index(drop=True)
    # Create fake index column for styling purposes (need it to start from 1 instead of 0)
    neg_control_ordered.insert(0, "row", range(1, len(neg_control_ordered) + 1))

    # Load sample abundance table
    abundance = pd.read_csv(f"{args.input_dir}/results/{args.sample_name}_downsampled.fastq_rel-abundance.tsv", sep="\t")
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
    abundance_ordered = abundance_ordered.reset_index(drop=True)

    # Merge abundance and assignment if prob_score is given
    if args.prob_score:
        abundance_assignment = abundance_ordered.merge(assignment_summary, on="tax id", how="left")
    else:
        abundance_assignment = abundance_ordered

    # Create fake index column for styling purposes (need it to start from 1 instead of 0)
    abundance_assignment.insert(0, "row", range(1, len(abundance_ordered) + 1))

    # Merge abundance and alignment based metrics
    if args.alignment_metrics:
        alignment_metrics = get_alignment_metrics(args.sample_name, args.input_dir)
        abundance_assignment = abundance_assignment.merge(alignment_metrics, on="tax id", how="left")
        neg_control_alignment_metrics = get_alignment_metrics(args.neg_control, args.input_dir)
        neg_control_ordered = neg_control_ordered.merge(neg_control_alignment_metrics, on="tax id", how="left")

    # Spike species
    with open(args.config, "rb") as f:
        config = tomllib.load(f)

    highlight = set(config.get("spike_species", []))
    normalising_spike_species = config.get("normalising_spike_species")

    # Define function for spike species
    def highlight_species(row):
        if row["species"] in highlight:
            return ["background-color: #ddd6fe"] * len(row)
        return [""] * len(row)

    # Define function for unique species not found in negative control
    def unique_species(row):
        if row["species"] not in neg_control_ordered["species"].values:
            return ["background-color: #dcfce7"] * len(row)
        return [""] * len(row)

    def low_abundance(row):
        if row["abundance"] < LOW_ABUNDANCE_CUTOFF:
            return ["color: #9ca3af"] * len(row)
        return [""] * len(row)

    # Define function for species normalized against spike
    def normalised_abundance(row):
        # No spike configured -> do nothing
        if not normalising_spike_species:
            return [""] * len(row)

        # Get abundance of spike in sample
        sample_spike= abundance_ordered.loc[
            abundance_ordered["species"] == normalising_spike_species, "abundance"
        ]

        # Get abundance of spike in neg control
        control_spike = neg_control_ordered.loc[
            neg_control_ordered["species"] == normalising_spike_species, "abundance"
        ]

        # Get abundance of species in neg control
        control_match = neg_control_ordered.loc[
            neg_control_ordered["species"] == row["species"], "abundance"
        ]

        # If any of these are empty, we can't do the calculation, so we return no highlight
        if sample_spike.empty or control_spike.empty or control_match.empty:
            return [""] * len(row)

        # Normalise (species / spike) in sample and control, then compare
        sample_ratio = row["abundance"] / sample_spike.iloc[0]
        control_ratio = control_match.iloc[0] / control_spike.iloc[0]

        if control_ratio > 0 and sample_ratio > 25 * control_ratio:
            return ["background-color: #dcfce7"] * len(row)

        return [""] * len(row)

    # Apply functions for spike species and unique species
    styled_abundance = (abundance_assignment.style
        .apply(normalised_abundance, axis=1)
        .apply(unique_species, axis=1)
        .apply(highlight_species, axis=1)
        .apply(low_abundance, axis=1)
        .format({
            "estimated read counts": "{:.0f}",
            "abundance": "{:.2%}",
            "median aligned identity": "{:.2%}",
            "median aligned coverage": "{:.2%}",
            })
        .set_properties(subset=["row"], **{
            "background-color": "#f2f2f2",
            "font-weight": "bold"
        })
        .hide(axis="index")
    )

    # Convert to html table
    html_table = styled_abundance.to_html(index=False, border=0, escape=False)

    # Apply function for spike species
    styled_neg_control = (neg_control_ordered.style
        .apply(highlight_species, axis=1)
        .apply(low_abundance, axis=1)
        .format({
            "estimated read counts": "{:.0f}",
            "abundance": "{:.2%}",
            "median aligned identity": "{:.2%}",
            "median aligned coverage": "{:.2%}",
            })
        .set_properties(subset=["row"], **{
            "background-color": "#f2f2f2",
            "font-weight": "bold"
        })
        .hide(axis="index")
    )
    # Convert to html table
    neg_control_html_table = styled_neg_control.to_html(index=False, border=0, escape=False)

    legend_lines = [
        '<span class="inline-block w-4 h-4 bg-purple-200 mr-2 border"></span> Purple rows indicate spike species<br>',
        '<span class="inline-block w-4 h-4 bg-green-100 mr-2 border"></span> Green rows indicate species not found in negative control or with a normalised abundance ratio &ge; 25.<br>'
    ]
    # Optional line for probability score
    if args.prob_score:
        legend_lines.append(
            '<span class="not-italic text-sm"><span class="font-bold text-black">median/mean probability*</span>: Median/Mean probability of the assigned taxon across reads</span>'
        )

    legend_html = '<p class="text-gray-700 italic mt-2">'+ "".join(legend_lines) + '</p>'

    legend_neg_html = """
    <p class="text-gray-700 italic mt-2">
        <span class="inline-block w-4 h-4 bg-purple-200 mr-2 border"></span>
        Purple rows indicate spike species
    </p>
    """

    # Save date
    today = date.today().strftime("%Y-%m-%d")

    # Load software version
    with open(f"{args.input_dir}/pipeline_info/software_versions.yml") as v:
        software_versions = yaml.safe_load(v)

    # Load MultiQC JSON
    with open(f"{args.input_dir}/multiqc/multiqc_data/multiqc_data.json") as f:
        multiqc_data = json.load(f)

    trana_version = software_versions["Workflow"]["genomic-medicine-sweden/TRANA"]

    html = template.render(
        css = css_content,
        table = html_table,
        neg_control_table = neg_control_html_table,
        legend = legend_html,
        legend_neg = legend_neg_html,
        today = today,
        pipeline_version = trana_version,
        multiqc_data  = multiqc_data,
        input_dir = Path(args.input_dir).name,
        sample_name = args.sample_name,
        neg_control = args.neg_control
    )

    with open(args.output_file, "w") as f:
        f.write(html)


if __name__ == "__main__":
    main()
