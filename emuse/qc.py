import json

import yaml


def load_multiqc_data(path):
    with open(path) as f:
        return json.load(f)


def trana_version(software_versions_path):
    with open(software_versions_path) as v:
        software_versions = yaml.safe_load(v)
    return software_versions["Workflow"]["genomic-medicine-sweden/TRANA"]
