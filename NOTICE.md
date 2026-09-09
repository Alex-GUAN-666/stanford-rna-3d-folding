# Attribution

This repository contains my notes and code backups from our Stanford RNA 3D
Folding project, alongside a CPU toolkit added after the competition.

## Upstream projects

Our workflows use Protenix, DRfold2, trRNA/RNAformer-related components,
SPOT-RNA, Rosetta/PyRosetta, and USalign. Credit for these models and tools belongs
to their original authors. Source links are in [docs/REFERENCES.md](docs/REFERENCES.md).
The historical code calls external implementations; model weights, source
snapshots, dependency wheels, and competition datasets are not bundled here.
Their respective licenses and data terms apply.

The trRNA backup includes RNAformer, SPOT-RNA, and folding components, but I do
not have a complete record of their versions. The CPU toolkit's `trrosettarna`
label is a routing label; it does not identify an exact upstream release.

## Project materials

The [evidence folder](evidence/README.md) contains my certificate, our final
leaderboard screenshot, a technical report, and validation records. The current
report summarizes the competition notes and the later software work.

[historical/reference/](historical/reference/) contains code backups, with the
original hashes and changes to the reference copies recorded in the
[historical manifest](historical/source_manifest.json). The supported
`rna_folding` package and examples are separate, post-competition utilities.

Competition names, logos, certificates, source documents, and historical code
retain their respective attribution and terms. This project is not affiliated
with or endorsed by Kaggle, Stanford, or the upstream authors. This notice does
not grant a repository-wide software license.
