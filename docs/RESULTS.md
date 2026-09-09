# Results

Our team, **Team GZSL**, finished Stanford RNA 3D Folding with a **0.42295 private leaderboard score**, placing **76th out of 1,516 teams** and earning a **bronze medal**. That is approximately the top 5% of the field (`76 / 1516 = 5.013%`).

| Record | Details |
| --- | --- |
| [Kaggle certificate](../evidence/certificate.jpg) | Alex GUAN; 76 / 1,516 teams; bronze; awarded September 24, 2025 |
| [Leaderboard screenshot](../evidence/private_leaderboard.png) | Team GZSL; rank 76; score 0.42295; bronze icon |
| [Official leaderboard](https://www.kaggle.com/competitions/stanford-rna-3d-folding/leaderboard) | Competition results |

## Score context

At submission time, we could see public leaderboard feedback; the private score was revealed later. **0.409 was a later-test public leaderboard result. Our final private result was 0.42295.** These use different evaluation subsets, so their difference is not a measured improvement on a common test set. I have not recovered the submission files needed to determine whether both scores came from the same CSV.

The earlier documents contain several other scores:

| Source | Experiment described | Score | Context |
| --- | --- | --- | --- |
| Technical report | Protenix baseline → improved pipeline | 0.334 → 0.368 | Development/validation comparison; the report alternates between a training subset and a 12-target validation set, so the exact split remains unresolved. |
| Original solution notes | Protenix baseline → adapted Protenix | 0.367 → 0.373 LB | The original notes do not identify the leaderboard phase. |
| Original solution notes | Adapted DRfold2 | 0.363 LB | Leaderboard phase unspecified. |
| Original solution notes | trRNA | 0.320 LB | Leaderboard phase unspecified. |
| Original solution notes, with my later clarification | Protenix + trRNA + DRfold2 | 0.409 public LB | Later-test public result; separate from the final private result above. |

The notes call the adapted models “fine-tuned”; their matching training records and weights are missing. The two Protenix comparisons should be kept separate rather than joined into one improvement curve. The report's runtime comparison of approximately one hour to 20 minutes also lacks enough hardware, target and configuration detail to serve as a reproducible benchmark.

The [solution PDF](../evidence/solution_notes.pdf) summarizes these records and the implementation. It is a current write-up, not an additional experiment. Original document hashes are in the [source manifest](../evidence/source_manifest.json), and the public/private clarification is recorded in [score context](../evidence/score_context.json).

## Saved inference runs

| Backup | What remains in the original saved output |
| --- | --- |
| `0.381.ipynb` | 12 Protenix targets, a 2,515-row × 18-column intermediate DataFrame and nine trRNA replacements; the final concatenation cell ends in a `TypeError`. Earlier cells write CSVs, but those files are missing. |
| `drfold2-runnable.ipynb` | No saved outputs or executed cell counts. |
| Two Protenix scripts | Source only; no accompanying run logs. |

`0.381` is a filename, not a verified score. That notebook has `VALIDATION=False`; its saved folding diagnostics do not provide a corresponding TM-score or leaderboard result. The [historical notes](HISTORICAL_NOTES.md) cover the configurations and known bugs, and the [archive manifest](../historical/source_manifest.json) records the saved-output findings.

## What runs in this repository

The current CPU toolkit has a synthetic end-to-end demo and automated software checks. It accepts structures generated elsewhere; it does not rerun the archived GPU models or reproduce the leaderboard score.

| Check | Scope |
| --- | --- |
| Synthetic demo | Produces a submission CSV from example coordinates without competition data or weights. |
| Routing tests | Check length boundaries and deterministic model preferences. |
| Geometry tests | Check known rigid transforms, overlap handling and complete residue coverage. |
| Input validation | Rejects missing, mismatched and invalid coordinate data. |
| Submission import/export | Preserves five supplied coordinate slots and residue mapping; rejects all-five-slots-zero targets. |
| Optional USalign comparison | Scores a supplied prediction/reference PDB pair; separate from the complete competition evaluator. |

For an actual run on your machine, follow [Verification](VERIFY.md). Current automated runs are available in [GitHub Actions](https://github.com/Alex-GUAN-666/stanford-rna-3d-folding/actions). Passing these checks verifies software behavior; RNA accuracy requires evaluation against experimental structures.

## Next experiments

My next priority is recovering the final `submission.csv` and its matching configuration. For a new benchmark, I would fix the target split, data cutoff and candidate budget, then compare a common-model baseline with length routing and candidate allocation. Per-target scores, length-group scores, failures, runtime and hardware should accompany the aggregate result. No such new accuracy benchmark is included here yet.
