# Project history and scope

I participated in Stanford RNA 3D Folding in 2025 as part of **Team GZSL**. We finished with **0.42295 on the private leaderboard**, **76 / 1,516 teams**, and a **bronze medal**. My work focused on coordination, experiment analysis and model integration. The result belongs to our team; the underlying predictors are the work of their upstream authors.

I put this repository together after the competition from the reports and code backups I still had. It contains an archive of inference experiments, a revised technical write-up and a CPU package for handling prediction files. The CPU package was added after the competition and is maintained separately from the archived notebooks.

## Repository contents

| Material | Location and purpose |
| --- | --- |
| Competition result | [Certificate](../evidence/certificate.jpg) and [leaderboard screenshot](../evidence/private_leaderboard.png) |
| Current write-up | [Solution PDF](../evidence/solution_notes.pdf), [text version](SOLUTION.md) and [method](METHOD.md) |
| Historical code | [Four inference backups](../historical/README.md), with configuration and defect notes in [Historical notes](HISTORICAL_NOTES.md) |
| Runnable CPU code | `rna_folding/`, `examples/`, the walkthrough notebook and tests |
| Source records | [Document manifest](../evidence/source_manifest.json) and [code manifest](../historical/source_manifest.json), including filenames, hashes and reference-copy changes |

The current solution PDF replaces the original two-page notes. It incorporates the code details and the final private score. The original notes and `RNA_IEEE_reportpdf.pdf` are summarized rather than included as current repository documents. The latter is a project report; its filename is not a claim of IEEE publication.

## Historical work and later additions

| Component | Historical backups | Current CPU package |
| --- | --- | --- |
| Prediction | Protenix + trRNA integration and a separate DRfold2 workflow | Reads externally generated coordinates; no neural inference |
| Candidate selection | PyRosetta ordering plus trRNA slot replacement; DRfold2 mean-score heuristic | Deterministic allocation across model pools; energy ordering within declared comparable groups |
| Length handling | Protenix truncation at 960 nt; DRfold2 480 nt windows with one-residue overlap and a 2400 nt cap | Configurable full-coverage windows, default 300 nt / 50 nt overlap |
| Coordinate assembly | Historical model-specific processing, including zero padding | Proper Kabsch alignment with at least three non-collinear overlap points |
| Validation | Partial saved notebook output and reports | Synthetic tests, submission checks and a repeatable CPU verification script |

The three-model strategy in the solution notes is broader than either surviving notebook. The Protenix + trRNA notebook does not invoke DRfold2, and the DRfold2 notebook does not run Protenix. The planner's length boundaries and geometric checks are later engineering choices, not recovered settings for the final submission.

## What is still missing

The archive is incomplete. I do not currently have the final submitted CSV, its exact run configuration, all upstream source snapshots or the matching model weights. The archived models have not been rerun as part of this repository's verification.

The older notes also mention fine-tuning, additional input features and MSA changes without their matching training artifacts. The surviving Protenix configurations have `use_msa=False`. Their 200-step notebook and 400-step scripts represent different inference settings. Those records are useful implementation references, but they do not resolve which exact setup produced every historical score.

[Results](RESULTS.md) keeps the older experiments separate from the final leaderboard result. [Historical notes](HISTORICAL_NOTES.md) documents the saved failure and other code issues in one place.

## Upstream credit

This project uses ideas and components from Protenix, DRfold2 and the RNAformer/SPOT-RNA/folding pipeline labeled trRNA in our backups. PyRosetta and USalign supply separate scoring and alignment functions. Their source, weights, datasets and licenses remain governed by their respective projects. See [References](REFERENCES.md) and [NOTICE](../NOTICE.md).
