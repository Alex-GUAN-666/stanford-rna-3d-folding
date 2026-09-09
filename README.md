# Stanford RNA 3D Folding

**Team GZSL · 76 / 1,516 teams · Final private score: 0.42295 · Bronze medal**

[![CPU verification](https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions/workflows/tests.yml/badge.svg)](https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions/workflows/tests.yml)

I'm Yuzhen Guan (Alex). This repository collects my team's inference experiments
from the 2025 Stanford RNA 3D Folding competition and the tools I've added since
then to inspect, combine, and validate saved predictions. My work on the team
focused on coordination, comparing experiments by sequence length, and integrating
model outputs.

[Technical report (PDF)](evidence/solution_notes.pdf) · [中文说明](START_HERE_ZH.md)
· [Run locally](docs/VERIFY.md) · [Competition](https://www.kaggle.com/competitions/stanford-rna-3d-folding)

## Task and approach

The task is to predict RNA structure from sequence. For each sequence of length
`L`, a submission contains five candidate structures, each represented by an
`(L, 3)` array of C1′ atom coordinates.

We explored pretrained Protenix, DRfold2, and trRNA components, with different
inference settings and candidate handling for different sequence lengths. The
saved code includes a Protenix + trRNA notebook and a separate DRfold2 notebook.
It captures several experiments rather than one final three-model pipeline.

I organized those backups under [`historical/`](historical/README.md). The current
[`rna_folding/`](rna_folding) package handles the next stage: checking residue
mapping, aligning overlapping segments, selecting five candidates, and exporting
a valid submission. It runs on saved coordinates and includes a synthetic CPU
demo. It does not run the pretrained models; the original weights, external model
code, and final scored submission are still needed to reproduce our competition run.

## Quick start

Clone the repository and use **Python 3.10 or later** in a virtual environment.
From the repository root:

```bash
python -m pip install .
python -m rna_folding demo --output-dir outputs/demo
python -m rna_folding validate --sequences outputs/demo/sequences.csv --submission outputs/demo/submission.csv
python scripts/verify_repo.py --output outputs/verification.json
```

The demo creates three synthetic targets, candidate coordinates, a routing plan,
`submission.csv`, and a selection audit. It needs only a CPU and NumPy. Choose a
new or empty demo output directory for each run. The verifier runs the tests and
checks the installed package, commands, and import/export round trip.

For environment setup on Windows, macOS, or Linux, see [local verification](docs/VERIFY.md).
The [walkthrough notebook](notebooks/01_workflow_walkthrough.ipynb) explains the
coordinate operations. These examples check software behavior; they do not measure
RNA prediction accuracy. [GitHub Actions](https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions)
also runs the CPU checks on Windows and Linux.

## Working with saved predictions

To import a five-candidate submission together with its sequence CSV:

```bash
python -m examples.import_submission --sequences data/test_sequences.csv --submission data/recovered_submission.csv --model historical_mixed --output-dir outputs/recovered
```

The importer checks target IDs, residue mapping, and coordinates, then writes
NumPy arrays and a candidate manifest. Use a fresh output directory.
`historical_mixed` labels a submission whose per-slot model sources are unknown.

To assemble predictions already described by a candidate manifest:

```bash
python -m rna_folding plan --sequences data/test_sequences.csv --output outputs/plan.json
python -m rna_folding assemble --sequences data/test_sequences.csv --manifest data/candidates/manifest.json --output outputs/submission.csv --audit outputs/selection_audit.json
python -m rna_folding validate --sequences data/test_sequences.csv --submission outputs/submission.csv
```

The default routing preference uses `trrosettarna` for 1–100 nt, `drfold2` for
101–480 nt, and `protenix` above 480 nt. I retained this policy from the written
solution as a configurable starting point; the historical notebooks use different
rules. `trrosettarna` is a working label, and the exact upstream identity of the
historical trRNA component remains unresolved. Candidate selection keeps five
separate structures and compares energy only within compatible score groups.
See [methods](docs/METHOD.md) and [input formats](docs/REPRODUCIBILITY.md) for details.

## Results

| Result | Score or placement |
| --- | --- |
| Final private leaderboard | **0.42295** |
| Team GZSL | **76 / 1,516 teams, bronze** |
| Later public leaderboard test | **0.409** |

Our final result is shown in the [leaderboard screenshot](evidence/private_leaderboard.png)
and [award certificate](evidence/certificate.jpg). At submission time, we could see
the public score; the private result became available when the final leaderboard
was released. The public and private numbers cover different evaluation sets,
so they are not a before-and-after improvement measurement.

I've kept the earlier development scores and their available context in
[results](docs/RESULTS.md). The filename `0.381.ipynb` has no matching score record.
The current CPU tools have not reproduced the leaderboard score.

## Repository guide

| Path | Contents |
| --- | --- |
| [`historical/`](historical/README.md) | Two inference notebooks and two Protenix scripts |
| [`docs/HISTORICAL_NOTES.md`](docs/HISTORICAL_NOTES.md) | Inference settings, saved-run observations, and known issues |
| [`rna_folding/`](rna_folding) | Sequence validation, geometry, candidate selection, and export |
| [`examples/`](examples) and [`notebooks/`](notebooks) | Import/export utilities and the CPU walkthrough |
| [`tests/`](tests) and [`scripts/verify_repo.py`](scripts/verify_repo.py) | Automated checks and local verification |
| [`evidence/`](evidence/README.md) | Technical report, competition results, and file inventories |
| [`docs/PROVENANCE.md`](docs/PROVENANCE.md) | Project history and the relationship between the backups and current tools |

## Acknowledgments

The competition result belongs to **Team GZSL**. Our inference work builds on
pretrained models and public implementations by their upstream authors. Model
architectures and weights are not my work; see [references](docs/REFERENCES.md)
and [NOTICE](NOTICE.md) for attribution. Competition data and model assets must
be obtained separately under their original terms.
