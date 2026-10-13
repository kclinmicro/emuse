# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [dev]

### Added

- [#28](https://github.com/kclinmicro/emuse/pull/28) Added `--version` (`-v`) to print the installed Emuse version and exit
  without requiring report arguments, enabling version capture in Nextflow
  workflows (by @rannick).

### Changed
- [#24](https://github.com/kclinmicro/emuse/pull/24) Changed rows in Summary statistics, adding corrected rows for Number of reads before and after qc and downsampling (by @AnnaNoren)
- [#26](https://github.com/kclinmicro/emuse/pull/26) Version is read from `emuse/__init__.py` (by @ryanjameskennedy)

### Fixed
- [#25](https://github.com/kclinmicro/emuse/pull/25) Fixed bundled config path handling (reported by @ryanjameskennedy, fixed by @AnnaNoren)
- [#26](https://github.com/kclinmicro/emuse/pull/26) Fixed version mismatch across 3 files (by @ryanjameskennedy)

## [1.1.0]

### Documentation

- [#15](https://github.com/kclinmicro/emuse/pull/15) Documented plans to
  support generating reports directly from raw Emu output, not just from a
  TRANA pipeline run (see the README "Roadmap" section) (by @samuell).

### Changed

- [#15](https://github.com/kclinmicro/emuse/pull/15)/[#17](https://github.com/kclinmicro/emuse/pull/17)
  Restructured the project into an installable `emuse` Python package to
  prepare for packaging as a bioconda recipe (by @samuell):
  - Moved `make_report.py`, `templates/`, `static/`, `configs/`, and
    `taxonomy.tsv` into `emuse/` (bundled under `emuse/data/`).
  - Moved `make_report_test.py` to `tests/test_make_report.py`.
  - Added an `emuse` console script entry point.
  - Renamed the distribution/CLI name from `16s-report`/`make_report.py` to
    `emuse`.
- [#17](https://github.com/kclinmicro/emuse/pull/17) `pyproject.toml` now
  declares a proper `[build-system]`, package data, and MIT license
  metadata (by @samuell).

### Added

- [#17](https://github.com/kclinmicro/emuse/pull/17) Added `LICENSE` (MIT)
  (by @samuell).
- [#17](https://github.com/kclinmicro/emuse/pull/17) Added `MANIFEST.in` for
  sdist packaging (by @samuell).
- [#17](https://github.com/kclinmicro/emuse/pull/17) Added a draft bioconda
  recipe at `recipe/meta.yaml` (by @samuell).
- [#15](https://github.com/kclinmicro/emuse/pull/15)/[#17](https://github.com/kclinmicro/emuse/pull/17)
  Added `THIRD_PARTY_LICENSES.txt` and a README "Acknowledgements" section
  crediting [Emu](https://github.com/treangenlab/emu) (MIT licensed), whose
  output this project parses and reports on (by @samuell).

## [1.0.0]

Initial release of the report generator for the TRANA 16S rRNA taxonomic
profiling pipeline.

### Added

- [#1](https://github.com/kclinmicro/emuse/pull/1) Added the HTML report
  generator that turns Emu/TRANA abundance and read-assignment output into
  a self-contained, styled report (by @AnnaNoren).
- [#2](https://github.com/kclinmicro/emuse/pull/2) Added a `--output_file`
  CLI parameter to control where the generated report is written
  (by @samuell).
- [#3](https://github.com/kclinmicro/emuse/pull/3) Added alignment-based
  metrics (percent identity and coverage) computed from Emu's alignment
  output (by @samuell).
- [#4](https://github.com/kclinmicro/emuse/pull/4)/[#5](https://github.com/kclinmicro/emuse/pull/5)
  Refined assignment metrics and added a customizable negative control
  table in the report (by @samuell, @AnnaNoren).
- [#8](https://github.com/kclinmicro/emuse/pull/8)/[#9](https://github.com/kclinmicro/emuse/pull/9)
  Added normalised abundance calculation relative to a spike species, with
  the spike species later made configurable via `config.toml` instead of
  being hardcoded (by @AnnaNoren).
- [#6](https://github.com/kclinmicro/emuse/pull/6) Added a test suite with
  fixtures based on realistic sequences (by @samuell).

### Fixed

- [#10](https://github.com/kclinmicro/emuse/issues/10)/[#11](https://github.com/kclinmicro/emuse/pull/11)
  Removed insertions (instead of deletions) from the
  `query_alignment_length` calculation, which had skewed coverage figures
  (by @samuell).
- [#13](https://github.com/kclinmicro/emuse/issues/13)/[#14](https://github.com/kclinmicro/emuse/pull/14)
  Report layout now adapts its dimensions to the browser window instead of
  using a fixed size (by @AnnaNoren).
