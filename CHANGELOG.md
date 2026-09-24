# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Documentation

- Documented plans to support generating reports directly from raw Emu
  output, not just from a TRANA pipeline run (see the README "Roadmap"
  section).
### Changed

- Restructured the project into an installable `emuse` Python package to
  prepare for packaging as a bioconda recipe:
  - Moved `make_report.py`, `templates/`, `static/`, `configs/`, and
    `taxonomy.tsv` into `emuse/` (bundled under `emuse/data/`).
  - Moved `make_report_test.py` to `tests/test_make_report.py`.
  - Added an `emuse` console script entry point.
  - Renamed the distribution/CLI name from `16s-report`/`make_report.py` to
    `emuse`.
- `pyproject.toml` now declares a proper `[build-system]`, package data, and
  MIT license metadata.

### Added

- `LICENSE` (MIT).
- `MANIFEST.in` for sdist packaging.
- Draft bioconda recipe at `recipe/meta.yaml`.
- `THIRD_PARTY_LICENSES.txt` and a README "Acknowledgements" section crediting
  [Emu](https://github.com/treangenlab/emu) (MIT licensed), whose output this
  project parses and reports on.

## [0.1.0]

- Initial release of the report generator for the TRANA 16S rRNA taxonomic
  profiling pipeline.
