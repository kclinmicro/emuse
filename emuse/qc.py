import json

import yaml


def load_multiqc_data(path):
    with open(path) as f:
        return json.load(f)


def trana_version(software_versions_path):
    with open(software_versions_path) as v:
        software_versions = yaml.safe_load(v)
    return software_versions["Workflow"]["genomic-medicine-sweden/TRANA"]


NANOPLOT_PROCESSED_STATS = {
    "reads_after_filtering": "Number of reads_fastq",
    "mean_read_length": "Mean read length_fastq",
    "mean_read_quality": "Mean read quality_fastq",
    "median_read_length": "Median read length_fastq",
    "median_read_quality": "Median read quality_fastq",
    "read_length_n50": "Read length N50_fastq",
    "stdev_read_length": "STDEV read length_fastq",
    "total_bases": "Total bases_fastq",
}


def nanoplot_stats(multiqc_data, sample_name, stage):
    key = f"{sample_name}_nanoplot_{stage}"
    for section in multiqc_data.get("report_general_stats_data", []):
        if key in section:
            return section[key]
    return {}


def summary_stats(multiqc_path, sample_name):
    multiqc_data = load_multiqc_data(multiqc_path)
    unprocessed = nanoplot_stats(multiqc_data, sample_name, "unprocessed")
    processed = nanoplot_stats(multiqc_data, sample_name, "processed")
    stats = {"reads_before_filtering": unprocessed.get("Number of reads_fastq")}
    stats.update({name: processed.get(key) for name, key in NANOPLOT_PROCESSED_STATS.items()})
    return stats
