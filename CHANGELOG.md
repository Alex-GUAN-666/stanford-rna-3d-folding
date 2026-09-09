# Changelog

## Documentation and package naming

- Use `rna-folding-toolkit` as the distribution name; the `rna-folding` command and `rna_folding` import stay the same.
- Consolidate historical inference settings and known issues in `docs/HISTORICAL_NOTES.md`.
- Update the technical report and installation guides.

## Unreleased

- Simplified the project overview, technical notes, and running instructions.
- Consolidated historical implementation details in `docs/HISTORICAL_NOTES.md`.

## September 9, 2026

- Added English and Chinese installation and verification guides.
- Added `scripts/verify_repo.py` to check a regular package installation, the
  unit suite, CPU commands, submission import and reassembly, and optional
  walkthrough notebook execution.
- Fixed console-command discovery on Windows and added installed-source
  consistency checks. Notebook kernels are closed after verification.
- Expanded CI to Ubuntu and Windows on Python 3.10 and 3.12, with Jupyter
  notebook execution on Ubuntu / Python 3.12 and downloadable run reports.
- Corrected the DRfold2 reference notebook's schema version to 4.5 to support
  its existing cell IDs.
- Added the final leaderboard record: Team GZSL, 76th, private score 0.42295.
  The separate later-test public result remains 0.409.
- Added a revised technical report, editable Markdown source, and PDF build
  script.

## 0.2.0

- Added two historical notebooks and two Protenix scripts with source hashes.
  Reference notebooks have cleared outputs; one workspace-deletion call is
  disabled.
- Documented the Protenix/trRNA concatenation failure, DRfold2 configuration,
  and external dependencies in the historical notes.
- Added a submission importer that preserves five coordinate slots by residue
  ID, validates input coverage and values, and protects existing output files.
- Added ten importer tests, including numeric roundtrips and reordered rows.

## 0.1.0

- Added sequence routing, candidate selection, overlapping-coordinate assembly,
  and submission validation.
- Added a synthetic walkthrough notebook and competition notes.
