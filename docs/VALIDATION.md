# Validation

## GitHub Actions

The [verification workflow](https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions/workflows/tests.yml) runs on every push. Open the latest run for `main` and check its commit and conclusion. Each job saves its command output and verification report.

| Environment | Checks run |
| --- | --- |
| Ubuntu, Python 3.10 | 45 tests and 11 command checks |
| Ubuntu, Python 3.12 | 45 tests and 11 command checks; walkthrough's 7 code cells execute in Jupyter |
| Windows, Python 3.10 | 45 tests and 11 command checks |
| Windows, Python 3.12 | 45 tests and 11 command checks |

The workflow installs the package, runs commands outside the checkout, compares
the installed modules with the source, and saves verification reports. The
submission checks include an import/reassembly roundtrip with numeric comparison.
A successful run verifies the commit shown on its page; a pending run is not a passing result.

## Local checks

The [fresh-install report](../evidence/reader_verification.json) records the
commands, outputs, environment, and source hashes from a separate virtual
environment. All 45 tests and 11 command checks passed.

| Demo property | Result |
| --- | --- |
| Synthetic targets | 3, with lengths 72, 160, and 520 residues |
| Available candidates | 21, or 7 per target |
| Selected candidates | 5 per target |
| Submission dimensions | 752 rows and 18 columns |

The [local notebook check](../evidence/notebook_verification.json)
executed all seven code cells sequentially in a shared namespace because kernel
sockets were unavailable there. Actual Jupyter execution is covered by the
Ubuntu / Python 3.12 CI job above.

## Coverage and limits

The tests cover routing boundaries, candidate IDs, coordinate validity, missing
or duplicate residues, path handling, energy comparison groups, rigid alignment,
window coverage, PDB formatting, and TM-score output parsing. Importer tests also
cover reordered rows, all-zero targets, numeric roundtrips, and output-directory
protection.

The historical GPU notebooks and scripts have been inspected but have not been
rerun in this environment. Their external dependencies and known issues are
recorded in [HISTORICAL_NOTES.md](HISTORICAL_NOTES.md). All three notebooks pass
schema validation; that check covers file format only.

The USalign wrapper is tested with a stub. Real USalign execution, neural-model
inference, and the historical competition evaluator are outside these checks.
The synthetic geometry plot and coordinate roundtrip measure numerical behavior,
not RNA prediction accuracy.

## Run locally

Follow [the installation guide](VERIFY.md) or [中文指南](VERIFY_ZH.md), then run:

```bash
python scripts/verify_repo.py --output outputs/verification.json
```

For the optional walkthrough check, install `requirements-notebooks.txt` and
add `--notebook`. The command writes the result and source hashes
to its report.
