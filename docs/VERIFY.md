# Install and verify

[中文说明](VERIFY_ZH.md) · [Command and input reference](REPRODUCIBILITY.md)

The CPU toolkit runs without a GPU or model weights. This guide covers a fresh
installation, the synthetic demo, and the walkthrough notebook. For our
competition result and its source records, see [Results](RESULTS.md).

## 1. Clone and install

Install Git and **Python 3.10 or later**. Cloning and installing dependencies
requires an internet connection.

### Windows PowerShell

Use the virtual environment's Python directly:

```powershell
git clone https://github.com/Alex-GUAN-666/stanford-rna-3d-folding.git
cd stanford-rna-3d-folding
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install .
.\.venv\Scripts\python.exe scripts/verify_repo.py --output outputs/verification.json
```

If `py` is unavailable, check `python --version` and use
`python -m venv .venv` instead.

### macOS / Linux

```bash
git clone https://github.com/Alex-GUAN-666/stanford-rna-3d-folding.git
cd stanford-rna-3d-folding
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install .
.venv/bin/python scripts/verify_repo.py --output outputs/verification.json
```

If Linux reports that `venv` is unavailable, install the distribution's Python
virtual-environment package first. After changing package code or switching
commits, repeat `pip install .` in the same environment before verifying.

## 2. Check the result

The verifier runs the test suite and installed command-line tools outside the
checkout. It covers planning, sequence summaries, model-input export, candidate
assembly, submission validation, and an import/reassembly roundtrip.

| Check | Expected result |
| --- | --- |
| Unit suite | 45 tests pass |
| Command checks | 11 pass |
| Demo targets | 3 synthetic sequences: 72, 160, and 520 residues |
| Submission | 752 rows, 18 columns, 5 coordinate sets per target |
| Import/reassembly | Coordinate values preserved by target and residue ID |

Success returns exit code 0 and writes `"passed": true` to
`outputs/verification.json`. The report includes command output, the environment,
source hashes, and checks that the installed package matches the checkout.
Failures return a nonzero exit code and include an `error` in the report.

Each verification run uses temporary demo data. To retain a demo for inspection,
run these commands from the repository root:

```powershell
.\.venv\Scripts\python.exe -m rna_folding demo --output-dir outputs/my_demo
.\.venv\Scripts\python.exe -m rna_folding validate --sequences outputs/my_demo/sequences.csv --submission outputs/my_demo/submission.csv
```

On macOS/Linux, replace `.\.venv\Scripts\python.exe` with `.venv/bin/python`.
Open `submission.csv`, `plan.json`, and `audit.json` in `outputs/my_demo` to
inspect the coordinates, length routing, and candidate selections. Choose a new
or empty output directory for subsequent runs.

## 3. Run the notebook

Install the optional notebook dependencies and enable notebook execution:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-notebooks.txt
.\.venv\Scripts\python.exe scripts/verify_repo.py --notebook --output outputs/verification_notebook.json
```

On macOS/Linux, use `.venv/bin/python`. This executes the walkthrough's **seven
code cells** in a Jupyter kernel and saves `executed_walkthrough.ipynb` beside the
report, including outputs and the synthetic geometry plot.

For an interactive session:

```powershell
.\.venv\Scripts\python.exe -m jupyterlab
```

Open `notebooks/01_workflow_walkthrough.ipynb`, select the environment where the
package is installed, and run the cells in order.

## 4. Compare with CI

[GitHub Actions](https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions)
runs the CPU checks on Ubuntu and Windows with Python 3.10 and 3.12. The
Ubuntu / Python 3.12 job also executes the notebook. Select the run for the
commit you are using and open its logs or downloadable verification reports.

To include your source revision in a bug report:

```bash
git rev-parse HEAD
```

Attach that commit ID and your verification JSON. Recorded runs are summarized
in [VALIDATION.md](VALIDATION.md).

## Scope

These checks cover the CPU toolkit and synthetic geometry. They do not generate
neural-model predictions or recompute the competition score. The historical GPU
notebooks need separate model environments and assets; setup details and known
issues are in [HISTORICAL_NOTES.md](HISTORICAL_NOTES.md).

The optional `evaluate` command also needs a separately installed USalign binary
and real prediction/reference structures. Its automated tests check the wrapper
with a stub; they do not execute the competition evaluator.
