# Historical inference notes

These notes describe the four source backups in [`historical/`](../historical/README.md). They are useful for understanding the experiments, but depend on Kaggle inputs, external model directories and checkpoints that are not included. Use the current [CPU verification guide](VERIFY.md) for a supported local starting point.

## Protenix + trRNA notebook

[`protenix_trrna_0381.ipynb`](../historical/reference/protenix_trrna_0381.ipynb) came from `0.381.ipynb`. Its pipeline is:

1. Generate five Protenix candidates with Protenix 0.4.6, `model_v0.2.0`, 10 cycles, 200 diffusion steps, seed 42 and `use_msa=False`.
2. Score the all-atom Protenix candidates with PyRosetta, order them by ascending score and extract coordinates using the notebook's atom-index mapping.
3. Run the trRNA/RNAformer pipeline, including SPOT-RNA and folding components, for targets of at most 300 nt.
4. Keep the first four Protenix candidates and replace slot 5 with the **first trRNA candidate**.

The notebook calculates a minimum trRNA score elsewhere, but does not use its candidate index for this replacement. It does not invoke DRfold2 or implement the current planner's short/medium/long routes. Exact versions of its external trRNA components are missing.

The original saved output contains 12 Protenix targets and nine trRNA replacements. The final cell raises `TypeError` because `pd.concat(results)` receives a mixture of dictionaries and DataFrames. Intermediate CSV writes occur earlier, but their resulting files are not in the archive. `VALIDATION=False`, and there is no matching score record for the filename `0.381`.

## DRfold2 notebook

| Setting | Archived behavior |
| --- | --- |
| Configuration | Modified `cfg_97` source, fp32 |
| Checkpoint count | ≤100 nt: 20; 101–200: 10; 201–480: 5; >480: 5 |
| Segmentation | 480 nt windows, stride 479: one shared residue |
| Maximum processed length | 2400 nt, followed by padding for omitted positions |
| Default ranking | `GET_CENTER=True`: distance from the mean geometric heuristic score |
| Energy/DR-specific scoring and optimization | Disabled by default |
| Output | Five sets of C1′ coordinates |

A single shared point cannot determine a unique rotation between two 3D windows. The default score is a custom geometric heuristic, not a calibrated molecular free energy; choosing candidates near its mean is also different from choosing the lowest scores. No saved outputs or executed cell counts remain in this notebook.

## Protenix scripts

The [Kaggle script](../historical/reference/protenix_kaggle_reference.py) and [Windows script](../historical/reference/protenix_windows_reference.py) both use 10 cycles, 400 diffusion steps, seed 0 and `use_msa=False`. They are different configurations from the notebook above. Both are inference scripts, without a training or optimizer loop and without accompanying run logs.

The Kaggle version references undefined `time0` in its validation path. Its exception handler also uses `target_id == ...` where assignment was intended. The Windows version corrects those two statements but still requires external files and model dependencies.

## Known issues

| Issue | Effect and current treatment |
| --- | --- |
| Whole-workspace cleanup | Original cell 20 of `0.381.ipynb` deletes `/kaggle/working`, including earlier outputs. The reference copy replaces that call with an explicit `RuntimeError`. |
| Mixed result types | The saved final concatenation fails. Disabling workspace deletion does not repair this separate bug. |
| Truncation and padding | The Protenix notebook truncates beyond 960 nt; DRfold2 beyond 2400 nt. Zero-filled positions are not model predictions. The Protenix path also stores a result before subsequently padding its coordinates. |
| Residue and atom mapping | Positional coordinate mixing, fixed atom indices, append-mode CSV writes and zero fallbacks need checks against the exact upstream formats. |
| Missing environment | External source snapshots, checkpoints and coordinate resources are needed before attempting a GPU rerun. Installing this repository's CPU dependencies does not supply them. |

The reference notebooks retain their source cells apart from the disabled cleanup call. Outputs and execution counts are cleared for publication; original cell numbers remain in metadata. The Python scripts retain their source bytes under clearer filenames. The [manifest](../historical/source_manifest.json) records hashes, transformations and the original saved error.

## Reusing an archived prediction

A recovered `submission.csv` can enter the current toolkit without another GPU run. The [archive README](../historical/README.md#reuse-a-saved-submission) contains the import commands. The importer checks target IDs, residue numbering and sequence agreement, preserves all five coordinate slots, and rejects all-zero targets. It does not infer which model produced a coordinate set or assign a leaderboard score.

The most useful remaining files are the final submission CSV, its Kaggle input list, exact source revisions and checkpoint names, and the corresponding score record. These would connect the archived settings to a reproducible run.
