# References and dependencies

These are the upstream projects used by the historical workflows or the current
CPU utilities. For new experiments, record the source commit and model checkpoint;
current releases may differ from the 2025 competition environment.

## Competition and evaluation

| Source | Use |
| --- | --- |
| [Stanford RNA 3D Folding](https://www.kaggle.com/competitions/stanford-rna-3d-folding) | Competition overview |
| [Competition data](https://www.kaggle.com/competitions/stanford-rna-3d-folding/data) | Data access and usage rules |
| [Official leaderboard](https://www.kaggle.com/competitions/stanford-rna-3d-folding/leaderboard) | Final standings |
| [DasLab Part 1 tools](https://github.com/DasLab/k1_tools) | C1′ coordinates, Ångström units, and five-conformation CSV format |
| [Organizer TM-score implementation](https://github.com/DasLab/k1_tools/blob/main/ribonanza_tmscore.py) | Native-conformer comparison, best-candidate scoring, and aggregation |
| [USalign](https://github.com/pylelab/USalign) | Structural alignment for the optional `evaluate` command |
| [USalign command-line source](https://github.com/pylelab/USalign/blob/master/USalign.cpp) | `-mol RNA`, `-atom`, and optional `-TMscore 1` arguments |

The organizer reference compares five predictions with available native
conformers, takes the best candidate/reference score for each target, and averages
across targets in a split. It uses the USalign score normalized by the reference
length. The current public script is a useful metric reference; reproducing the
historical evaluation also requires its original version, inputs, and split.
The local `evaluate` command scores one PDB pair.

## Folding models

| Project | Source | Notes |
| --- | --- | --- |
| Protenix | [bytedance/Protenix](https://github.com/bytedance/Protenix) | Used by the historical inference scripts |
| Protenix inputs | [Inference JSON format](https://github.com/bytedance/Protenix/blob/main/docs/infer_json_format.md) | RNA job schema used by the input exporter |
| Protenix setup | [Training and inference instructions](https://github.com/bytedance/Protenix/blob/main/docs/training_inference_instructions.md) | Upstream environment and execution instructions |
| DRfold2 | [leeyang/DRfold2](https://github.com/leeyang/DRfold2) | The backup also requires modified modules and external checkpoints |
| trRosettaRNA | [Official download page](https://yanglab.qd.sdu.edu.cn/trRosettaRNA/download/) | Reference for the toolkit's `trrosettarna` routing label |
| trRosettaRNA2 | [YangLab-SDU/trRosettaRNA2](https://github.com/YangLab-SDU/trRosettaRNA2/) | A separate project; not a drop-in identification of the historical trRNA branch |

My original notes call the third branch “trRNA”. The backup includes RNAformer,
SPOT-RNA, and folding components, but I do not have all source versions needed to
map it to a specific trRosettaRNA release. See
[historical implementation notes](HISTORICAL_NOTES.md) for the actual imports and
configuration. The earlier report also discusses MSA and Ribonanza-derived
features; the available backups do not contain a complete training pipeline.

## Upstream terms

The linked Protenix repository provides Apache-2.0 terms for its code and model
parameters; DRfold2 includes MIT permission text; trRosettaRNA2 provides
Apache-2.0 terms. USalign has its
[own license](https://github.com/pylelab/USalign/blob/master/LICENSE).
Check the terms of the exact code, model assets, and optional dependencies you
use. Terms for the older trRosettaRNA package are not recorded here. See
[NOTICE.md](../NOTICE.md) for attribution.

## Project records

- [My certificate](../evidence/certificate.jpg): bronze, 76th of 1,516 teams.
- [Final leaderboard screenshot](../evidence/private_leaderboard.png): Team GZSL,
  rank 76, private score 0.42295. [Score context](../evidence/score_context.json)
  also records the separate later-test public result of 0.409.
- `Stanford RNA 3D Folding-solution(1).pdf`: my original two-page solution notes.
- `RNA_IEEE_reportpdf.pdf`: the earlier project report, used as background for
  the current documentation.
- [Historical code backups](../historical/README.md): two notebooks and two
  Protenix scripts.

The original document names and hashes are retained in the
[source manifest](../evidence/source_manifest.json). The original PDFs are
summarized in the maintained [technical report](../evidence/solution_notes.pdf),
with [Markdown source](SOLUTION.md) and a
[PDF build script](../scripts/build_solution_pdf.py). The project's timeline is
in [PROVENANCE.md](PROVENANCE.md).
