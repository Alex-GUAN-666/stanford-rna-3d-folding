# Method

I worked on combining pretrained RNA structure predictors and handling their outputs under different sequence lengths and inference constraints. Our competition notes describe a three-model strategy using Protenix, DRfold2 and a pipeline called trRNA. The surviving code contains a Protenix + trRNA notebook and a separate DRfold2 notebook.

The current `rna_folding` package is a later addition: it handles coordinate files, candidate selection and submission validation on a CPU. It does not run the neural models. This page separates the historical settings from that package's behavior.

## Prediction format

For an RNA sequence of length `L`, each candidate is an `(L, 3)` array with one C1′ coordinate per nucleotide, in Ångströms. A submission has five complete candidates per target, with residue identity and one-based residue numbering. Rotating or translating a structure changes its coordinate frame, not its fold quality.

## Historical inference

| Setting | Protenix + trRNA notebook | DRfold2 notebook |
| --- | --- | --- |
| Model configuration | Protenix 0.4.6, `model_v0.2.0`, 10 cycles, 200 diffusion steps, five samples, seed 42, `use_msa=False` | Modified external `cfg_97` source; fp32 |
| Length-dependent work | trRNA runs for targets ≤300 nt | 20 / 10 / 5 / 5 checkpoints for lengths ≤100 / 101–200 / 201–480 / >480 nt |
| Candidate selection | Sort five Protenix candidates by PyRosetta score; replace slot 5 with the first trRNA candidate for targets ≤300 nt | With `GET_CENTER=True`, order candidates by distance from the mean geometric heuristic score and keep five |
| Long targets | Truncate Protenix inputs beyond 960 nt and pad omitted coordinates with zeros | 480 nt windows with one-residue overlap; truncate beyond 2400 nt and pad omitted frames |

The trRNA replacement does not compare scores across models: it takes the first trRNA candidate, even though another part of the notebook calculates a minimum score. DRfold2's default rule selects candidates near the mean heuristic score, rather than the five lowest-energy structures. Neither operation averages different folds together.

The two Protenix script backups use 10 cycles, **400** diffusion steps and seed **0**, with MSA disabled. These are separate configurations from the notebook. The backups contain inference code; the training artifacts mentioned in our older notes are missing. [Historical notes](HISTORICAL_NOTES.md) document the dependencies, saved run output and known defects.

## Current toolkit: routing

The planner uses these model preferences:

| Sequence length | Preferred model label |
| --- | --- |
| 1–100 nt | `trrosettarna` |
| 101–480 nt | `drfold2` |
| >480 nt | `protenix` |

I used the written solution's length groups as the starting point for this planner. Assigning 100 nt to the short route and 480 nt to the middle route resolves the overlapping boundaries in those notes. This is a planning policy, not the exact routing in either archived notebook or a measured claim that each model is best within its range.

`trrosettarna` is the toolkit label for the source called trRNA in our notes. The archived pipeline imports RNAformer, SPOT-RNA and folding components; its exact external source snapshots remain missing. Candidate generation happens outside the package.

## Current toolkit: window assembly

Default planning uses **300 nt windows with 50 nt overlap**. The window length follows the report; the overlap and boundary handling were added for this toolkit. They differ from the historical notebooks above. Ranges are zero-based and half-open. The last window shifts left to preserve its full length where possible, so its overlap can exceed 50 nt.

Windows can be predicted in unrelated coordinate frames. The assembler:

1. Matches residues in the shared overlap.
2. Fits a proper Kabsch rigid transform, excluding reflections.
3. Places the incoming window in the assembled frame.
4. Averages aligned contributions equally at shared residues and preserves complete coverage.

The first window anchors the frame. Each new window needs at least three non-collinear shared points. Gaps, ambiguous fits and incomplete coverage raise errors. The demo checks this behavior using known transforms of synthetic coordinates.

This operation assembles local geometry. It does not perform energy minimization, recover contacts between distant windows or guarantee a physically valid global fold. Reference coordinates are used only during evaluation, never to align inference windows.

## Current toolkit: five-candidate selection

Each candidate has a target ID, source model, coordinate file and optional energy metadata. Selection is deterministic:

- The routed model goes first. Remaining known models follow `trrosettarna`, `drfold2`, `protenix`, then other labels alphabetically.
- Round-robin allocation takes one candidate from each nonempty model pool until five are selected.
- A target needs at least five distinct candidate IDs; the package does not duplicate missing candidates.
- Within a model pool, finite energies can reorder candidates only within the same explicit `energy_group`. Other groups and unscored candidates retain their relative slots.

An energy group should identify the scoring method, units, direction and version. Raw scores from different models or force fields are not automatically comparable. The package supplies no physical energy function.

Round-robin allocation preserves multiple model sources. It does not measure geometric diversity, and distinct IDs can still describe very similar structures. Comparing allocation rules on the same saved candidate pool is a useful next experiment.

## Validation and evaluation

CSV validation checks target and residue mapping, coverage, five coordinate triplets and finite values. The submission importer also rejects targets whose five slots are all zero. These are input and output checks, not tests of biological plausibility.

The optional USalign wrapper compares one predicted PDB with one reference PDB and returns the second, reference-normalized TM-score. It requests C1′ structural alignment with `-mol RNA -atom " C1'"`; `--by-residue` adds `-TMscore 1` for residue-index correspondence. It is a pairwise diagnostic. The full competition metric additionally handles native conformers, candidate selection, invalid residues and aggregation; see [references](REFERENCES.md).

A real accuracy benchmark still needs fixed targets and a data cutoff, exact upstream versions and weights, saved predictions, and the official evaluator. The final competition result and the checks available in this repository are documented in [Results](RESULTS.md).
