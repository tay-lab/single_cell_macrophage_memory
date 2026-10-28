# single_cell_macrophage_memory

Code and data to reproduce the figures of Wang AG, Mercado-Vásquez JG, Son M, et al.,
"Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback" (Cell Systems).

## Requirements
- MATLAB R2024a (Update 8); all scripts were run and tested in this version.
- Tested on Windows 10 Pro (64-bit), Intel Core i5-13600K (14 cores), 64 GB RAM. No special hardware is needed; a multi-core CPU only shortens the optional recomputation.
- Toolboxes: Deep Learning, Statistics and Machine Learning, Signal Processing, Parallel Computing (`parfor`; without it the loops run serially and are ~10× slower).
- Third-party plotting functions are included in `functions/` (`Violin.m`, `violinplot.m`, `violin3.m` and variants; see `functions/license.txt`).
- R (Bioconductor: limma, edgeR, EGSEA, org.Mm.eg.db, gplots) only for `R/Fig4_RNAseq.R`; nichenetr 2.2.1.1 (R 4.5.2) for `R/Fig4BC_S4AB_NicheNet.R`, which downloads the mouse NicheNet v2 networks (Zenodo 7074291) on first run.

## How to run
1. Set MATLAB's current folder to this folder. All paths in the scripts are relative to it (each script adds `functions/`, `simulation_results/`, `classifier_results/` and `data/flow/` to the path), so the folder can be placed anywhere.
2. Run any script (`>> Fig5`), or run it section by section; figures open in MATLAB windows.

By default the scripts **load saved results**, so every figure is drawn in seconds and is identical to the paper:
- `simulation_results/`: 1000-cell model simulations for Figures 5–7 and the Jensen–Shannon scans for Figure S6 (only the protocols and variables each figure plots).
- `classifier_results/`: 5-fold cross-validation results of the LSTM, SVM, decision-tree and random classifiers (Figures 2, 3, S3).

To recompute instead, set the flag at the top of the script to `true` (`runSimulation` in `Fig5–7.m` and `FigS6.m`, `runTraining` in `Fig2.m` and `Fig3.m`); the saved files are then overwritten. Runtimes with a 14-worker pool: Fig5 ~1.5 min, Fig6 ~40 s, Fig7 ~1.2 min, FigS6 ~7 min, classifiers a few minutes per figure.

## Scripts and figures
| Script | Figures |
|---|---|
| `Fig1.m` | Figure 1 |
| `Fig2.m` | Figure 2 |
| `Fig3.m` | Figure 3 |
| `Fig4.m` | Figure 4D–I (also prints n and p-values for every comparison) |
| `Fig5.m` | Figure 5B–F (model) and 5C (TLR9 flow cytometry and Gamma fit) |
| `Fig6.m` | Figure 6C–H |
| `Fig7.m` | Figure 7B–G |
| `FigS2.m` | Figure S2 |
| `FigS3.m` | Figure S3 |
| `FigS4.m` | Figure S4C–F |
| `FigS6.m` | Figure S6B–C (Figure S6A was assembled from these results) |
| `R/Fig4_RNAseq.R` | Figure 4A, 4J, S4G (RNA-seq re-analysis) |
| `R/Fig4BC_S4AB_NicheNet.R` | Figure 4B–C, S4A–B (NicheNet ligand activity and regulatory potential) |

Figures S1 and S5 are schematics. R scripts are also run from this folder (the NicheNet script writes to `output_nichenet/`).

## Data
| Path | Content |
|---|---|
| `functions/genExtFeat.mat` | single-cell NF-κB traces and features (Figures 1–3, S2, S6) |
| `single_cell_features/` | traces and features of the cytokine-blocking experiments (Figure 4, S4); `genBlockFeat.mat` holds them as one single-cell table, `data4` (built by `functions/genBlockFeat.m`) |
| `data/flow/` | TLR9 flow cytometry: gated single-cell export used for Figure 5C, and the raw `.fcs` files (BD LSRFortessa, 03-Aug-2023) |
| `data/rnaseq/` | RNA-seq counts (GEO: GSE260996) and CellTalkDB ligand–receptor pairs |

The earlier analysis of the RNA-seq data (GEO: GSE260996) is at Zenodo (DOI: 10.5281/zenodo.14337874). This package (code and single-cell data) will be archived at Zenodo with its own DOI.
