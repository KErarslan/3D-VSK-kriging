# 3D-VSK: 3D Variogram Surface Kriging

Reproducible code accompanying:

**Erarslan, K.** "3D Variogram Surface Kriging: A PSD-Consistent Covariance
Framework for Anisotropic Ore Estimation." *Computers & Geosciences*
(submitted).

## Overview

This repository contains the Python implementation used to produce the
results, tables, and figures in the manuscript. The method (3D-VSK) replaces
discrete directional variogram model selection — which can produce an
indefinite (non-PSD) kriging covariance matrix under common anisotropy
conditions — with a continuous ellipsoidal covariance surface that is
PSD-consistent by construction.

Three covariance formulations are compared throughout:

1. **Classical Discrete** — nearest-direction model selection (standard practice)
2. **Classical Bilinear** — bilinear interpolation between directional models (Erarslan, 2001)
3. **3D-VSK (LMC-Ellipsoidal)** — the proposed method

on two case studies: a small, dense gold dataset (Kalgoorlie, Australia,
n=20) and a large, sparse lignite dataset (Seyitömer, Türkiye, n=191).

## Repository Structure

```
.
├── notebooks/
│   └── 3DSurfaceKriging.ipynb       ← main reproducible notebook (run this)
├── data/
│   └── synthetic_seyitomer_generator.py
├── requirements.txt
├── LICENSE
└── README.md
```

## Quick Start

The simplest way to run everything is via Google Colab:

1. Open `notebooks/3DSurfaceKriging.ipynb` in [Google Colab](https://colab.research.google.com/).
2. Run all cells top to bottom (Runtime → Run all).
3. Figures and result tables are written to `figures/`, `figures_sli/`,
   `results/`, and `params/` in the Colab session storage.

To run locally instead:

```bash
git clone https://github.com/<your-username>/3D-VSK-kriging.git
cd 3D-VSK-kriging
pip install -r requirements.txt
jupyter notebook notebooks/3DSurfaceKriging.ipynb
```

## Data Availability

**Kalgoorlie (Au) dataset.** The Kalgoorlie Northern Zone drill hole data is
included directly (hardcoded) in the notebook. It is derived from the
publicly available quarterly activities report of Riversgold Limited (2022),
with two modifications described in the manuscript (Section 3.1): the
analysis is limited to 20 RC drill holes, and spatial coordinates and grades
have been partially adjusted to preserve the spatial correlation structure
while obscuring the absolute field positions and exact assay values. These
modifications are consistent with standard practice for methodological
publications using proprietary exploration data.

**Seyitömer (lignite) dataset.** The real Seyitömer-Aslanlı dataset was
obtained during a site visit to a completed (closed) mining operation. It
is **not** the author's data to redistribute publicly, and is therefore not
included in this repository. It remains available from the corresponding
author upon reasonable request, consistent with the manuscript's Data
Availability Statement.

To preserve full reproducibility of the methodology without disclosing this
proprietary survey data, `data/synthetic_seyitomer_generator.py` generates a
**synthetic** replacement dataset with the same sample size (n=191), spatial
extent, directional variogram structure (nugget/sill/range per direction),
and marginal statistics as the real dataset (see the script's docstring for
full methodology). Running the notebook on this synthetic dataset reproduces
the same *qualitative* result reported in the manuscript — Classical
Discrete and Classical Bilinear formulations produce indefinite (non-PSD)
covariance matrices, while 3D-VSK remains PSD-consistent — but the exact
numerical values (R², RMSE, λ_min) will differ from the published tables,
since the underlying data is different.

## Requirements

See `requirements.txt`. Tested with Python 3.10+.

## Citation

If you use this code, please cite the manuscript (full citation will be
added upon publication) and, where relevant, the original method this work
builds on:

> Erarslan, K., 2001. Three dimensional variogram modeling and kriging. In:
> Xie, H., Wang, Y., Jiang, Y. (Eds.), Proceedings of the 29th International
> Symposium on Application of Computers and Operations Research in the
> Mineral Industries (APCOM 2001), Beijing. A.A. Balkema, Rotterdam, pp. 51-56.

## License

See `LICENSE`.

## AI Assistance Declaration

Portions of this codebase were developed with the assistance of Claude
(Anthropic), used as a coding and refactoring tool. All scientific content,
methodology, and validation of results are the sole responsibility of the
author.
