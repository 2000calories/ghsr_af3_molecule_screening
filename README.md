# Finding Best Peptides and Natural Molecules for Human GHS-R Using Boltz-2

This repository provides a practical Boltz-2 screening workflow for **human GHS-R (GHSR1a, UniProt Q92847)** on Google Colab.

It builds Boltz YAML inputs from a curated ligand CSV, runs batch predictions, parses model outputs into one ranking table, and generates summary plots plus top-pose PDB exports.

## What this pipeline does

- Screens mixed ligand classes against human GHS-R:
  - peptides (modeled as a second protein chain),
  - small molecules (modeled as ligand chain with Boltz affinity head).
- Produces `results/ranking.csv` with unified metrics.
- Generates:
  - `results/plots/top20_affinity.png`
  - `results/plots/iptm_vs_affinity.png`
  - `results/plots/class_distribution.png`
- Copies top structures into `results/top_poses/` for quick 3D review.

## Important modeling caveats

- Boltz-2 affinity prediction is best suited for **small molecules**; peptide affinity values should not be trusted as primary ranking signals.
- Peptides are ranked from structural confidence proxies (primarily `iptm`, `interface_pae_mean`).
- Native ghrelin has an **octanoylated Ser3** PTM that is not directly represented by plain amino-acid sequence inputs. This repo includes a proxy row (`ghrelin_1_28_octanoyl_proxy`) and expects manual curation of an explicit chemistry-aware representation.
- Treat model rankings as a triage filter, not a final activity claim.

## Project layout

- `data/target_ghsr.fasta`: receptor sequence (Q92847).
- `data/ligand_library.csv`: curated ligands (editable).
- `scripts/build_inputs.py`: CSV -> Boltz YAML generator.
- `scripts/run_screen.py`: Boltz batch runner with resume/skip support.
- `scripts/parse_results.py`: output JSON parser -> ranking CSV.
- `scripts/visualize.py`: plots + top pose exports.
- `notebooks/boltz2_ghsr_screen.ipynb`: Colab workflow notebook.

## Quickstart (Colab Pro recommended)

1. Open `notebooks/boltz2_ghsr_screen.ipynb` in Colab.
2. Install dependencies:
   - `pip install -U boltz pandas matplotlib pyyaml py3Dmol`
3. Mount Drive and set:
   - `BOLTZ_CACHE=/content/drive/MyDrive/boltz_cache`
4. Ensure repo is available under `/content/ghsr_af3_molecule_screening`.
5. Run:
   - `python scripts/build_inputs.py`
   - `python scripts/run_screen.py --msa-mode server --recycling-steps 3 --diffusion-samples 1`
   - `python scripts/parse_results.py`
   - `python scripts/visualize.py`

## Editing the ligand library

`data/ligand_library.csv` columns:

- `name`: unique ligand identifier.
- `type`: `small_molecule` or `peptide`.
- `source`: free-text provenance (`natural`, `endogenous`, `positive_control`, etc.).
- `smiles`: required for `small_molecule`.
- `sequence`: required for `peptide`.
- `notes`: optional curation comments.

Rows missing required `smiles`/`sequence` are intentionally skipped by `build_inputs.py`.

## Ranking interpretation

- `affinity_pred_value_log10IC50uM`: lower is better (small molecules only).
- `affinity_pred_IC50_uM`: `10^(affinity_pred_value_log10IC50uM)`.
- `affinity_probability_binary`: binary-like binder confidence from Boltz affinity head.
- `iptm`, `ptm`, `interface_pae_mean`: structural confidence/proximity metrics for all inputs.

Suggested readout strategy:

- **Small molecules**: rank by affinity first, then check `iptm` and pose plausibility.
- **Peptides**: rank by `iptm` and `interface_pae_mean` only.

## Local execution (optional)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/build_inputs.py
python scripts/run_screen.py
python scripts/parse_results.py
python scripts/visualize.py
```

## Expected outputs

- `results/ranking.csv`
- `results/plots/*.png`
- `results/top_poses/*.pdb`

## Notes on curation quality

The default CSV intentionally includes a mixed-confidence starter set. For production screening, fill all missing SMILES/sequences with validated structures before launching long runs.
