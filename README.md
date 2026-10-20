# Emuse

This repository contains Emuse, a small reporting tool for generating a static
HTML report from taxonomic abundance results produced by
[Emu](https://github.com/treangenlab/emu), currently as part of the
[TRANA](https://github.com/genomic-medicine-sweden/TRANA) 16S rRNA taxonomic
profiling pipeline. TRANA uses [EMU](https://github.com/treangenlab/emu) for
species-level taxonomic abundance estimation from full-length 16S reads.
Support for generating reports directly from a standalone
`emu abundance` run (without TRANA) is planned for a future release; see
[Roadmap](#roadmap).


The script parses pipeline outputs and renders a human-readable summary report
using a Jinja2 HTML template and CSS styling.

## Repository Contents

| File                                  | Description                                                                                      |
| ------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `emuse/make_report.py`                | Python script that parses pipeline output and generates the report                               |
| `emuse/data/configs/config.toml`      | Default configuration file for customizing report generation                                     |
| `emuse/data/templates/report.html.j2` | Jinja2 template used to render the HTML report                                                   |
| `emuse/data/static/style.css`         | CSS styling for the report                                                                       |
| `emuse/data/taxonomy.tsv`             | TSV file with mapping of taxonomic ids to species                                                |
| `output/`                             | Directory where the generated HTML report (`report.html`) will be saved after running the script |
| `README.md`                           | Project documentation                                                                            |
| `pyproject.toml`                      | Configuration file including tool dependecies for the project.                                   |
| `tests/test_make_report.py`           | Minimal test python script for alignment metrics                                                 |
| `recipe/meta.yaml`                    | Draft bioconda recipe                                                                            |
| `THIRD_PARTY_LICENSES.txt`            | License and attribution for third-party tools whose output this project parses (Emu)             |

## Setup

This project uses dependencies defined in pyproject.toml and requires Python
≥3.11. You can set up the environment using Conda, Pixi, or pip.

### Dependencies

- python >=3.11,<=3.13
- jinja2
- pandas
- pysam
- pyyaml

### Installing

```bash
pip install .
```

This installs the `emuse` command-line tool.

## Usage

Run the report generator using the installed command:

```bash
emuse --input-dir <results_directory> \
                 --output-file <output_file> \
                 --sample-name <sample_name> \
                 --neg-control <negative_control_name> \
                 [--config config.toml] \
                 [--prob-score] \
                 [--alignment-metrics]
```

### Example
Generating a minimal report for 16S analysis:
```bash
emuse \
  --input-dir results/sample_01 \
  --output-file output/sample_01_report.html \
  --sample-name sample_01 \
  --neg-control neg_control \
  --config config.toml
```

Generating a minimal report for ITS analysis:
```bash
emuse \
  --input-dir results/sample_01 \
  --output-file output/sample_01_report.html \
  --sample-name sample_01 \
  --neg-control neg_control \
  --config config.toml \
  --analysis-type its
```

Generating a report with probability scores and alignment metrics:
```bash
emuse \
  --input-dir results/sample_01 \
  --output-file output/sample_01_report.html \
  --sample-name sample_01 \
  --neg-control neg_control \
  --config config.toml \
  --prob-score \
  --alignment-metrics
```

### Arguments

| Argument              | Short | Required | Description                                                                                                 |
| ---------------       | ----- | -------- | -------------------------------------------------------                                                     |
| `--input-dir`         | `-i`  | Yes      | Path to the directory containing TRANA pipeline results                                                     |
| `--output-file`       | `-o`  | Yes      | Name of the output html report to be generated                                                              |
| `--sample-name`       | `-s`  | Yes      | Name of the sample to generate the report for                                                               |
| `--neg-control`       | `-n`  | Yes      | Name of the negative control sample                                                                         |
| `--config`            | `-c`  | No       | Path to configuration file (default: the bundled default `emuse/data/configs/config.toml`)                  |
| `--prob-score`        | `-p`  | No       | Include the generation and addition of probability scores in the report                                     |
| `--alignment-metrics` | `-m`  | No       | Include metrics based on the raw alignment of reads to the database (percent identity and percent coverage) |
| `--version`           | `-v`  | No       | Print the installed emuse version and exit                                                                  |
| `--analysis-type`     | `-a`  | No       | Type of analysis: 16s or its. Default is 16s                                                                  |

### Roadmap

Emuse currently expects a TRANA pipeline run directory. A future release will
add support for generating reports directly from the output of a plain `emu
abundance` call (without requiring TRANA).

### Customization

The report can be customized by editing the bundled configuration file or
creating a new `config.toml` file and supplying it through `--config` (`-c`).
When `--config` is omitted, the bundled default config
(`emuse/data/configs/config.toml`) is used.

The following settings control spike species and abundance highlighting in the
report:
| Setting                     | Description                                                                                                                                               | 
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `spike_species`             | List of species treated as spike species and highlighted in purple in the report.                                                                         |
| `normalising_spike_species` | Spike species used to normalise sample abundances when comparing them with negative control.                                                              |
| `enrichment_fold_cutoff`    | Fold-change threshold for highlighting enriched species in green. A value of 25 means the normalised abundance in the sample must be greater than 25 times the abundance in the negative control to be highlighted green. Default: 25.                                                                                                               |
| `low_abundance_cutoff`      | Relative-abundance threshold below which species are displayed in a lighter grey text, reducing visual emphasis. Default: 0.005 (0.5% relative abundance) |

Define all spike species in the `spike_species` list. The species used for
normalisation is specified separately using `normalising_spike_species`, which
must correspond to one of the entries in the `spike_species` list. If no
`spike_species` or `normalising_spike_species` are defined in the config, the
functionality of colouring is removed from the report. If `enrichment_fold_cutoff`
or `low_abundance_cutoff` is omitted from the config, the default value is used.

### Output

The script generates a **static HTML report** (`report.html`) in the `output/`
directory summarizing results from the
[TRANA](https://github.com/genomic-medicine-sweden/TRANA) pipeline. By default,
the report specifies 16S results, but can be changed to specify ITS results using
`--analysis-type its`. The report contains:

#### Summary Statistics

Key sequencing metrics, including number of reads before and after filtering
(including downsampling), mean/median read length and read quality (Phred score),
read length N50, standard deviation (STDEV) of read lengths, and total bases.

#### Sample Abundance Table

A table summarizing the abundance and taxonomic composition of the sample, with
the following columns:

| Column                    | Description                                           |
| ------------------------- | ----------------------------------------------------- |
| Abundance                 | Relative abundance of the taxon                       |
| Species                   | Assigned species                                      |
| Genus                     | Assigned genus                                        |
| Family                    | Assigned family                                       |
| TaxID                     | NCBI Taxonomy ID                                      |
| Estimated read counts     | Number of reads assigned                              |
| Median probability*       | Median probability of the assigned taxon across reads |
| Mean probability*         | Mean probability of the assigned taxon across reads   |
| Median aligned identity** | Median aligned percent identity                       |
| Median aligned coverage** | Median aligned percent coverage                       |

*: The two columns for median and mean probability is by default not included in
the table. These can be added by supplying the argument --prob_score (-p).
**: The two columns for median aligned identity and coverage is by default not
included in the table. These can be added by supplying the argument
--alignment-metrics (-m).

**Color coding:**

- **Purple rows** indicate spike species
- **Green rows** indicate species absent in negative control or species whose
normalised abundance in the sample exceeds the configured enrichment fold cutoff
relative to the negative control. Normalised abundance is calculated through
comparing the abundance of species to the abundance of spike species.

#### Negative Control Table

A table summarizing the negative control sample, with the same first 6 columns as the
previous table.

**Color coding:**

- **Purple rows** indicate spike species

## Citations
- [EMU](https://github.com/treangenlab/emu)
  > Kristen D. Curry et al., “Emu: Species-Level Microbial Community Profiling
  > of Full-Length 16S RRNA Oxford Nanopore Sequencing Data,” Nature Methods,
  > June 30, 2022, 1–9, https://doi.org/10.1038/s41592-022-015>
