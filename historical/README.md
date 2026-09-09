# Inference archive

I kept these four backups from the RNA competition work. They preserve the inference settings and integration logic, while the current CPU utilities live separately in `rna_folding/`.

| Original file | Archive copy | Notes |
| --- | --- | --- |
| `0.381.ipynb` | [Protenix + trRNA notebook](reference/protenix_trrna_0381.ipynb) | Source cells retained; whole-workspace deletion disabled |
| `drfold2-runnable.ipynb` | [DRfold2 notebook](reference/drfold2_reference.ipynb) | Source retained; original had no saved outputs |
| `protenix.py` | [Kaggle script](reference/protenix_kaggle_reference.py) | Original source bytes, renamed |
| `RNA代码.py` | [Windows script](reference/protenix_windows_reference.py) | Original source bytes, renamed |

These are historical references rather than supported entry points. They require external model code, checkpoints, specific Kaggle datasets and GPU environments. The current `requirements.txt` only installs dependencies for the CPU toolkit. Start with [Verification](../docs/VERIFY.md) to run that code locally.

## Archive preparation

The reference notebooks have cleared outputs and execution counts, reduced metadata and an archive notice. Original cell numbers remain in metadata. In the Protenix + trRNA notebook, original cell 20 calls `shutil.rmtree` on the whole `/kaggle/working` directory; the archive copy replaces this call with `RuntimeError`. Other historical defects remain documented rather than silently changed.

The [source manifest](source_manifest.json) records original and reference SHA-256 hashes, byte counts, transformations and saved-output findings. In particular:

- `0.381.ipynb` saved 12 Protenix target predictions and nine trRNA replacements, then failed in the final concatenation with a `TypeError`. The intermediate CSVs are missing.
- `0.381` is an experiment filename; there is no corresponding stored leaderboard score.
- The DRfold2 notebook has no saved execution output. Neither Python script has accompanying run logs.

[Historical notes](../docs/HISTORICAL_NOTES.md) covers exact settings, truncation, selection rules and known bugs. The archived models have not been rerun during the current toolkit verification.

## Reuse a saved submission

To inspect a saved prediction CSV without rerunning model inference:

```bash
python -m examples.import_submission --sequences data/test_sequences.csv --submission data/archived_submission.csv --model protenix-trrna --output-dir outputs/imported
python -m rna_folding assemble --sequences data/test_sequences.csv --manifest outputs/imported/candidates.json --output outputs/reassembled.csv --audit outputs/reassembled_audit.json
```

Use a new or empty output directory. The importer validates residue IDs and sequence agreement, preserves the five coordinate slots, rejects targets whose five slots are all zero, and records source hashes. `--model` is a label you supply; the importer cannot infer model identity from coordinates.

Upstream model code, scoring libraries and datasets retain their original authorship and terms. See [NOTICE](../NOTICE.md) for reuse information.
