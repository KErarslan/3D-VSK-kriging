# 3D-VSK: Four PSD-Consistent Covariance Formulations for Anisotropic Kriging

Reproducible code accompanying:

**Erarslan, K.** "3D Variogram Surface Kriging: Four PSD-Consistent Covariance
Formulations for Anisotropic Ore and Coal Deposits." *Natural Resources
Research* (submitted).

## Overview

Discrete or interpolated directional-variogram selection — the common
geostatistical practice of assigning the nearest-direction model at each lag
and azimuth — can render the resulting kriging covariance matrix indefinite
(non-PSD), a structural weakness documented for four decades but still
common in practice. This repository implements and benchmarks **four new
covariance formulations**, each carrying a genuine, structurally provable
PSD guarantee, against a reference heuristic baseline:

1. **3D-VSK-Elliptical** — a corrected geometric-anisotropy transform (true
   elliptical coordinate transform, PSD by a linear-map argument).
2. **CZM (Continuous Zonal Mixture)** — an exponentiated-Fourier-series
   parametrization of the directional-mixture framework of Allard, Senoussi
   & Porcu (2016, *Mathematical Geosciences*), providing the first practical
   fitting procedure for that framework.
3. **Nested-LMR** — a reparametrized classical nested Linear Model of
   Regionalization (unconstrained fitting via bₖ = exp(uₖ); same structure
   as the classical framework).
4. **KernelSum** — Allard et al.'s (2016) own finite kernel-sum
   parametrization, implemented here for equal-terms comparison against CZM.
5. **3D-VSK-empirical** (reference) — the original heuristic ellipsoidal
   anisotropy formula; PSD-consistency is only empirically observed, not
   structurally guaranteed (it loses PSD under a Gaussian kernel at one of
   the two case-study sites).

All five formulations are tested across three base kernels (spherical,
exponential, Gaussian) and two independent case studies: a small, strongly
anisotropic gold deposit (Kalgoorlie, Australia, n=20) and a larger, weakly
anisotropic lignite field (Seyitömer, Türkiye, n=191, two variables). A
central finding is that no single formulation is universally optimal —
predictive performance tracks each deposit's own anisotropy signature
(Nested-LMR wins at Kalgoorlie, CZM wins at Seyitömer).

## Repository Structure

```
.
├── 3DSurfaceKriging.ipynb            ← main reproducible notebook (run this)
├── data/
│   └── synthetic_seyitomer_generator.py
├── requirements.txt
├── CITATION.cff
├── LICENSE
└── README.md
```

## Quick Start

The simplest way to run everything is via Google Colab:

1. Open `3DSurfaceKriging.ipynb` in [Google Colab](https://colab.research.google.com/).
2. Run all cells top to bottom (Runtime → Run all).
3. Figures, benchmark tables, and fitted-parameter files are written to
   `figures/`, `results/`, and `params/` in the Colab session storage.

To run locally instead:

```
git clone https://github.com/KErarslan/3D-VSK-kriging.git
cd 3D-VSK-kriging
pip install -r requirements.txt
jupyter notebook 3DSurfaceKriging.ipynb
```

## Data Availability

**Kalgoorlie (Au) dataset.** The Kalgoorlie Northern Zone drill hole data is
included directly (hardcoded) in the notebook. It is derived from the
publicly available quarterly activities report of Riversgold Limited (2022),
with modifications described in the manuscript (Section 3.1): the analysis
is limited to 20 RC drill holes, and spatial coordinates and grades have
been partially adjusted to preserve the spatial correlation structure while
obscuring the absolute field positions and exact assay values.

**Seyitömer (lignite) dataset.** The real Seyitömer-Aslanlı dataset was
obtained during a site visit to a completed (closed) mining operation. It
is **not** the author's data to redistribute publicly, and is therefore not
included in this repository. It remains available from the corresponding
author upon reasonable request, consistent with the manuscript's Data
Availability Statement.

To preserve full reproducibility of the methodology without disclosing this
proprietary survey data, `data/synthetic_seyitomer_generator.py` generates a
**synthetic** replacement dataset with the same sample size (n=191), spatial
extent, and directional variogram structure (nugget/sill/range per
direction) as the real dataset. Running the notebook on this synthetic
dataset reproduces the same *qualitative* results reported in the
manuscript, but the exact numerical values (R², RMSE, λ_min) will differ
from the published tables, since the underlying data is different.

## Requirements

See `requirements.txt`. Tested with Python 3.10+.

## Citation

If you use this code, please cite the manuscript (full citation will be
updated upon publication; see `CITATION.cff`) and, where relevant, the two
works this study directly builds on:

> Allard, D., Senoussi, R., Porcu, E., 2016. Anisotropy models for spatial
> data. *Mathematical Geosciences* 48, 305–328.
> https://doi.org/10.1007/s11004-015-9594-x

> Erarslan, K., 2001. Three dimensional variogram modeling and kriging. In:
> Xie, H., Wang, Y., Jiang, Y. (Eds.), Proceedings of the 29th International
> Symposium on Application of Computers and Operations Research in the
> Mineral Industries (APCOM 2001), Beijing. A.A. Balkema, Rotterdam, pp. 51-56.

## License

See `LICENSE` (MIT).

## AI Assistance Declaration

Portions of this codebase were developed with the assistance of Claude
(Anthropic), used as a coding and refactoring tool. All scientific content,
methodology, and validation of results are the sole responsibility of the
author.
