# Running the toolkit

## Environment

The CPU toolkit requires **Python 3.10 or later** and NumPy. No GPU, model checkpoint, Kaggle credential, or competition download is needed for the demo.

From the repository root:

```bash
python -m venv .venv
```

Activate the environment on macOS/Linux:

```bash
source .venv/bin/activate
```

Or in Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the package:

```bash
python -m pip install -e .
```

Folding models use separate upstream environments. This installation covers the CPU toolkit.

## Run the end-to-end demo

```bash
python -m rna_folding demo --output-dir outputs/demo
```

The demo uses synthetic coordinate traces for target lengths **72, 160, and 520 nt**, covering all three routing groups. The pipeline selects five of seven candidates per target. No neural model runs in this example.

Choose a new or empty output directory for each demo run. Existing demo files are protected from accidental overwrite; use a different directory name when rerunning the command.

| Output | What to inspect |
| --- | --- |
| `outputs/demo/sequences.csv` | Target identifiers and RNA sequences. |
| `outputs/demo/candidates.json` | Candidate metadata and relative coordinate-file paths. |
| `outputs/demo/coordinates/` | One numeric `.npy` coordinate array per candidate. |
| `outputs/demo/plan.json` | Routing choices and sequence-window plans. |
| `outputs/demo/submission.csv` | Five complete C1′ coordinate sets per target. |
| `outputs/demo/audit.json` | Which supplied candidates were selected. |
| `outputs/demo/demo_summary.json` | Demo output and validation summary. |
| `outputs/demo/synthetic_trace.pdb` | A C1′ trace exported from synthetic geometry, for inspecting PDB formatting. |

Re-run individual stages:

```bash
python -m rna_folding plan --sequences outputs/demo/sequences.csv --output outputs/demo/plan.json
python -m rna_folding assemble --sequences outputs/demo/sequences.csv --manifest outputs/demo/candidates.json --output outputs/demo/submission.csv --audit outputs/demo/audit.json
python -m rna_folding validate --sequences outputs/demo/sequences.csv --submission outputs/demo/submission.csv
```

The output contains 752 residue rows and 18 columns. Validation checks identifiers, residue coverage, and coordinate values. The synthetic demo has no TM-score.

## Inspect sequences and walk through the notebook

To summarize the demo sequences, or a sequence CSV you obtained separately:

```bash
python -m examples.inspect_sequences --sequences outputs/demo/sequences.csv --output outputs/sequence_summary.json
```

The JSON reports target count, total residues, length minimum/median/maximum, routing counts, and nucleotide counts/proportions for the input file.

The [walkthrough notebook](../notebooks/01_workflow_walkthrough.ipynb) covers the CPU demo, candidate audit, sequence summaries, routing, rigid window assembly, and an optional 3D plot. For an interactive session:

```bash
python -m pip install -r requirements-notebooks.txt
python -m jupyterlab
```

Open `notebooks/01_workflow_walkthrough.ipynb` and select the environment in which the toolkit dependencies are installed.

## Export upstream model inputs

After running the demo, create input files in a new output directory:

```bash
python -m examples.export_model_inputs --sequences outputs/demo/sequences.csv --output-dir outputs/model_inputs
```

This writes FASTA files, `protenix_jobs.json`, and `window_manifest.json`. The manifest maps each job to its original target and zero-based, half-open residue range. Targets longer than 480 nt are exported in overlapping windows; other targets retain their full sequence. FASTA files can serve as inputs to the relevant upstream tool, while the JSON follows the documented Protenix RNA job schema.

The exporter writes inputs for external tools; model installation and execution are separate steps. Choose a fresh output directory for each export because existing named files are protected. Follow the upstream model instructions and record its version and configuration.

## Import an existing submission

To process predictions saved by another workflow, supply a sequence file and a matching submission CSV with five C1′ coordinate slots:

```bash
python -m examples.import_submission --sequences data/test_sequences.csv --submission data/recovered_submission.csv --model historical_mixed --output-dir outputs/recovered
```

Use a **fresh output directory**. The importer validates target/residue IDs, sequence identity, full coverage, coordinate columns, and finite values. It writes the five slots by ID to `.npy` arrays and creates `candidates.json`. It rejects a target whose five coordinate sets are all zero; other zero-padded tails need separate inspection.

The required `--model` argument records the prediction source. Use a label such as `historical_mixed` for mixed predictions without a model-to-slot mapping. Optional `--prefix` sets candidate names. The importer does not infer source models or assign energy values.

To inspect the imported candidates through the assembler:

```bash
python -m rna_folding assemble --sequences data/test_sequences.csv --manifest outputs/recovered/candidates.json --output outputs/recovered/reassembled_submission.csv --audit outputs/recovered/selection_audit.json
python -m rna_folding validate --sequences data/test_sequences.csv --submission outputs/recovered/reassembled_submission.csv
```

The assembler applies its own selection and output ordering. Keep the original prediction CSV alongside the new export to retain the original slot order.

## Supply real external predictions

1. Obtain permitted sequence data and install a model from its [upstream project](REFERENCES.md).
2. Record the upstream commit, checkpoint, input data version, inference options, seeds, and hardware.
3. Generate complete candidate structures. Extract the C1′ atom by atom name, chain, and residue correspondence. Do not assume a universal all-atom tensor index for C1′.
4. Save each candidate as an unpickled NumPy array of shape `(L, 3)`, where rows match the full target sequence in order and coordinates are in Ångströms.
5. Create a versioned manifest and run `assemble`, then `validate`.

A minimal sequence CSV is:

```csv
target_id,sequence
example,ACGUACGU
```

A manifest entry looks like this:

```json
{
  "schema_version": 1,
  "candidates": [
    {
      "candidate_id": "example_protenix_seed0",
      "target_id": "example",
      "model": "protenix",
      "path": "coordinates/example_seed0.npy"
    }
  ]
}
```

This one-entry example illustrates the schema; **assembly requires at least five distinct candidate IDs per target**. Use the demo’s complete `candidates.json` as a runnable example. Candidate IDs are globally unique. The target set must match the sequence file, and coordinate files must stay beneath the manifest directory. Absolute paths and paths that escape it are rejected. Unknown manifest fields and pickled arrays are not accepted.

The optional fields `energy` and `energy_group` are for a finite numeric score and its comparison group. Omit them when a trustworthy comparable score is unavailable. The implementation orders lower energies first only within compatible groups of the same source-model pool; it never compares raw scores across models. Give the group a meaningful name, such as the evaluator version and units, and keep the exact scoring configuration in a separate run log. Do not put arbitrary unrecognized fields into the strict manifest.

For canonical input, use uppercase `A`, `C`, `G`, and `U`. The toolkit rejects unsupported symbols instead of silently replacing modified nucleotides. If preprocessing a real dataset requires mapping modifications, document the mapping and its scientific limitations before producing the canonical sequence file.

Python helpers support single-chain C1′ PDB import/export and segment stitching. Stitching requires compatible overlapping predictions and a complete sequence plan.

## Optional pairwise structural evaluation

Install USalign separately from its [official repository](https://github.com/pylelab/USalign), then compare PDB files you already have:

```bash
python -m rna_folding evaluate --predicted path/to/predicted.pdb --reference path/to/reference.pdb --usalign path/to/USalign
```

On Windows, point `--usalign` to the appropriate executable, such as `C:/tools/USalign.exe`. The default command requests RNA structural alignment using C1′ atoms and returns the score normalized by the second, reference structure.

To explicitly use residue-index correspondence instead of default structural alignment:

```bash
python -m rna_folding evaluate --predicted path/to/predicted.pdb --reference path/to/reference.pdb --usalign path/to/USalign --by-residue
```

The optional flag adds USalign’s `-TMscore 1`. Preserve original residue indices in input PDBs when correspondence matters, and record the selected mode with the result. `read_c1_pdb` returns sequence and coordinates in file order, without residue identifiers; passing those arrays to `write_c1_pdb` numbers residues from 1. Use the original PDBs for residue-index evaluation when their numbering matters. A reference-based fit is for evaluation only; the submission assembler never receives reference structures.

This helper scores **one predicted/reference pair**. Reproducing the complete competition evaluation requires the correct reference conformers, masks, target split, and aggregation from the [organizer scoring reference](https://github.com/DasLab/k1_tools/blob/main/ribonanza_tmscore.py). 

## Software checks

Run the tests from the repository root:

```bash
python -m unittest discover -s tests -v
```

Test coverage focuses on scientific data contracts and error cases: coordinate shape/validity, routing boundaries, proper rigid alignment, window coverage, compatible energy ranking, and submission integrity.

For experiments, preserve per-target predictions and a run manifest with the data split and cutoff, model/weight hashes, software versions, seeds, model settings, candidate score definitions, target failures, wall time, hardware, and evaluator flags.

## Historical inference backups

The four backup files are listed in [historical/README.md](../historical/README.md).
For configuration details and known issues, see
[HISTORICAL_NOTES.md](HISTORICAL_NOTES.md). Original hashes and changes to the
reference copies are recorded in
[historical/source_manifest.json](../historical/source_manifest.json).

The Protenix + trRNA notebook has cleared outputs and a disabled
workspace-deletion cell. Its saved run ended with a concatenation error. The
DRfold2 notebook has no saved execution outputs. Neither workflow is part of
the supported CPU installation.

Rerunning these workflows requires:

- Protenix checkpoints and chemistry resources, plus the trRNA/RNAformer,
  SPOT-RNA, and folding modules and weights used by the first notebook.
- The modified DRfold2 `cfg_97`, RNALM2, and PotentialFold sources; its RCLM
  checkpoint, 20 structure checkpoints, coordinate templates, and configuration.
- Competition inputs and offline dependencies, with handling for failed and
  truncated targets.
- The final submission and matching evaluation inputs and settings to reproduce
  the historical score.

The available backups retain useful implementation details, but the complete
GPU environments and final prediction CSV are still missing.
