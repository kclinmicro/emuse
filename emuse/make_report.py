import argparse
from datetime import date
from importlib import resources
from jinja2 import Environment, FileSystemLoader
import json
import yaml
from pathlib import Path
import pandas as pd
import pysam
import re
import statistics
import tomllib

# Bundled package data (templates, CSS, taxonomy mapping, default config)
DATA_DIR = resources.files("emuse") / "data"

def main():
    argp = argparse.ArgumentParser()
    argp.add_argument("-i", "--input-dir", type=str, required=True, help="Path to the input directory containing results")
    argp.add_argument("-o", "--output-file", type=str, required=True, help="Path to the output report file")
    argp.add_argument("-s", "--sample-name", type=str, required=True, help="Name of the sample")
    argp.add_argument("-n", "--neg-control", type=str, required=True, help="Name of the negative control")
    argp.add_argument("-c", "--config", type=str, default="configs/config.toml", help="Path to config file. Default is configs/config.toml")
    argp.add_argument("-p", "--prob-score", action="store_true", help="Include probability score in the report")
    argp.add_argument("-m", "--alignment-metrics", action="store_true", help="Include metrics based on the raw alignment of reads to the database (percent identity and percent coverage)")

    args = argp.parse_args()
    if args.config is None:
        args.config = str(DATA_DIR / "configs" / "config.toml")

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

    # Load emu log and extract mapped value
    with open(f"{args.input_dir}/results/emu_logs/{args.sample_name}_emu_log.log") as f:
        emu_log = f.read()
    # Use regex to find the value between "mapped" and "sequences"
    match_mapped = re.search(r"mapped (.*?) sequences", emu_log)

    # If a match is found, extract the value; otherwise, set it to "N/A" (avoid errors if the pattern is not found)
    if match_mapped:
        mapped_value = match_mapped.group(1)
    else:
        mapped_value = "N/A"

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
        mapped_value = mapped_value,
        input_dir = Path(args.input_dir).name,
        sample_name = args.sample_name,
        neg_control = args.neg_control
    )

    with open(args.output_file, "w") as f:
        f.write(html)


def get_alignment_metrics(sample_name, input_dir):
    abundance_path = f"{input_dir}/results/{sample_name}_downsampled.fastq_rel-abundance.tsv"
    assignment_path = f"{input_dir}/results/{sample_name}_downsampled.fastq_read-assignment-distributions.tsv"
    alignments_path = f"{input_dir}/results/{sample_name}_downsampled.fastq_emu_alignments.sam"

    df_abundance_unsorted = pd.read_csv(abundance_path, sep="\t")
    df_abundance = df_abundance_unsorted.sort_values("abundance", ascending=False)

    df_reads = load_read_file(assignment_path, df_abundance)
    colnames_taxids = df_reads.columns

    align_file = pysam.AlignmentFile(alignments_path)

    alns_all = {}
    for aln in align_file:
        readid = aln.query_name  # Ex: b9bb144e-eb53-4509-8931-5f4477444a48
        refname = aln.reference_name  # Ex: 562:emu_db:23853
        taxid = str(aln.reference_name).split(":")[0]  # Ex: 562
        refid = aln.reference_id  # Ex: 23853

        if aln.is_secondary or aln.is_supplementary:
            # We don't count these
            continue

        if taxid not in alns_all:
            alns_all[taxid] = []
        alns_all[taxid].append(aln)

    taxtr = TaxTranslator()
    aln_infos = []
    i = 1
    for taxid in colnames_taxids:
        if taxid in alns_all:
            alns = alns_all[taxid]
            identities, coverages = collect_distribution(alns, align_file)
            median_id = statistics.median(identities)
            median_cov = statistics.median(coverages)
            abundance = float(df_abundance[df_abundance["tax_id"] == taxid]["abundance"].values[0])
            taxon = taxtr.taxid_to_label(taxid)
            aln_infos.append(
                {
                    "tax id": taxid,
                    "median aligned identity": median_id,
                    "median aligned coverage": median_cov,
                }
            )

    df_aln_metrics = pd.DataFrame(aln_infos)
    return df_aln_metrics


def load_read_file(readassmt_path, df_abundance):
    df_reads = pd.read_csv(readassmt_path, sep="\t", header=0)
    colnames_sorted = [cn for cn in df_abundance["tax_id"] if cn in df_reads.columns]
    df_reads = df_reads[colnames_sorted]
    return df_reads


def collect_distribution(alns, alignment_file):
    identities = []
    coverages = []
    for aln in alns:
        identity, coverage = get_align_stats(aln, alignment_file)
        identities.append(identity)
        coverages.append(coverage)
    return identities, coverages


def get_align_stats(alignment, alignment_file):
    CIGAROP_MATCH = 0
    CIGAROP_INS = 1
    CIGAROP_DEL = 2
    CIGAROP_REF_SKIP = 3
    CIGAROP_SOFTCLIP = 4
    CIGAROP_HARDCLIP = 5
    CIGAROP_PAD = 6
    CIGAROP_EQUAL = 7
    CIGAROP_DIFF = 8
    CIGAROP_BACK = 9

    cigar_stats = alignment.get_cigar_stats()[0]
    edit_distance = alignment.get_tag("NM") if alignment.has_tag("NM") else None

    cigar_stats = alignment.get_cigar_stats()[0]

    insertions = cigar_stats[CIGAROP_INS]
    deletions = cigar_stats[CIGAROP_DEL]
    soft_clips = cigar_stats[CIGAROP_SOFTCLIP]

    query_len = alignment.query_length
    query_alignment_len = alignment.query_alignment_length

    reference_len = alignment_file.get_reference_length(alignment.reference_name)
    reference_alignment_len = alignment.reference_length

    identity, ref_coverage = calculate_align_stats(query_len,
                                                   query_alignment_len,
                                                   reference_len,
                                                   reference_alignment_len,
                                                   insertions,
                                                   deletions,
                                                   edit_distance,
                                                   soft_clips)

    return identity, ref_coverage



def calculate_align_stats(query_len,
                          query_alignment_len,
                          reference_len,
                          reference_alignment_len,
                          insertions,
                          deletions,
                          edit_distance,
                          soft_clips):
    """
    Return identity and coverage against the reference sequence

    Note that the identity in this calculation differs from the one in BLAST.
    While BLAST's identity measure does not count indels and instead divides
    the number of identical bases over the number of aligned bases, the
    identity metric below includes each aligned column corresponding to a base
    pair position in the reference or read, and counts the number of matches
    (identical bases) across the total aligned length, which is the length
    where gaps in both the reference and query are included.

    Args:
        query_len: The length in base pairs of the query sequence.
        query_alignment_len: The length in base pairs of the part of the query
            sequence covered by the alignment.
        reference_len: The length in base pairs of the reference sequence.
        reference_alignment_len: The length in base pairs of the part of the
            reference sequence covered by the alignment.
        insertions: Number of insertions measured in base pairs.
        deletions: Number of deletions measured in base pairs.
        edit_distance: Edit distance corresponding to the NM tag in SAM files.
        soft_clips: Bases in the ends of the alignment which do not align.

    Returns:
        identity: Identity as matching bases across the alignment length which
            contains both insertions and deletions.
        coverage: Percent of bases of the reference length covered by the
            alignment.
    """

    # ------------------------------------------------
    # Calculate identity
    # ------------------------------------------------
    matches = 0
    mismatches = 0
    if edit_distance is not None:
        mismatches = edit_distance - insertions - deletions
        matches = query_alignment_len - insertions - mismatches

    # We can not normalize over query or reference length, as these
    # won't contain either insertions or deletions
    divisor = matches + mismatches + insertions + deletions

    identity = 0
    if divisor > 0:
        identity = matches / divisor

    # ------------------------------------------------
    # Calculate coverage
    # ------------------------------------------------
    coverage = 0
    if reference_len > 0:
        coverage = reference_alignment_len / reference_len

    return identity, coverage


class TaxTranslator(object):
    def __init__(self, taxonomy_path=None):
        if taxonomy_path is None:
            taxonomy_path = DATA_DIR / "taxonomy.tsv"
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


if __name__ == "__main__":
    main()
