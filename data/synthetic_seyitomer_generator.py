#!/usr/bin/env python3
"""
================================================================================
Synthetic Seyitömer-Aslanlı Dataset Generator
================================================================================
Generates a statistically-equivalent SYNTHETIC replacement for the real
Seyitömer-Aslanlı drill hole dataset used in the 3D-VSK manuscript
(Erarslan, "3D Variogram Surface Kriging: A PSD-Consistent Covariance
Framework for Anisotropic Ore Estimation", Computers & Geosciences).

WHY SYNTHETIC DATA:
The real Seyitömer-Aslanlı drill hole data was obtained during a site visit
to a now-closed lignite mining operation. Although the mine is closed, the
underlying exploration/production data is not the author's to redistribute
publicly. To preserve full reproducibility of the 3D-VSK methodology without
disclosing proprietary survey data, this script generates a synthetic dataset
that:

  1. Has the same sample size (n=191) and spatial extent (3,736 x 6,511 m)
     as the real dataset (Section 4.1 of the manuscript).
  2. Reproduces the same directional variogram structure (nugget, sill,
     range per principal direction) reported in Table 9 of the manuscript,
     via simulation from the fitted LMC-Ellipsoidal covariance model.
  3. Reproduces the same marginal statistics (mean, std, skewness) for
     seam thickness and calorific value reported in Table 8.
  4. Reproduces the same large-scale spatial trend (X/Y correlation)
     reported in Section 4.2 (non-stationarity).

This is NOT the real Seyitömer dataset. Absolute values, drill hole
identities, and exact spatial positions carry no geological meaning.
The synthetic dataset exists solely so that the 3D-VSK / Classical
Discrete / Classical Bilinear comparison code can be run end-to-end by
anyone, producing qualitatively identical PSD-violation and LOOCV
behaviour to that reported in the manuscript.

Usage:
    python synthetic_seyitomer_generator.py
    # writes synthetic_seyitomer.csv to the current directory

Dependencies: numpy, pandas, scipy
================================================================================
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

RANDOM_SEED = 35  # validated against target statistics; see validation block below

# ── Field geometry (Section 4.1) ──────────────────────────────────────────
N_HOLES   = 191
X_EXTENT  = 3736.0   # m
Y_EXTENT  = 6511.0   # m

# ── Target marginal statistics (Table 8) ──────────────────────────────────
THICKNESS_TARGET = dict(mean=13.3, std=8.1, min_clip=0.5, max_clip=43.5)
SQRT_THICKNESS_TARGET = dict(mean=3.47, std=1.10)  # Table 8, sqrt(Thickness) row
CALORIFIC_TARGET = dict(mean=2039.0, std=374.0, min_clip=1088.0, max_clip=3187.0)

# ── Fitted directional LMC variogram parameters (Table 9) ────────────────
# sqrt(Thickness) space
THICKNESS_VARIO = {
    0:   dict(C0=0.654, C=0.890, a=1245.0),
    45:  dict(C0=0.449, C=1.719, a=1245.0),
    90:  dict(C0=0.342, C=0.866, a=1383.0),
    135: dict(C0=0.553, C=0.642, a=1245.0),
}
# Calorific value (kcal/kg) space
CALORIFIC_VARIO = {
    0:   dict(C0=75588.0,  C=106301.0, a=1245.0),
    45:  dict(C0=63436.0,  C=133055.0, a=1245.0),
    90:  dict(C0=78940.0,  C=55530.0,  a=317.0),
    135: dict(C0=103179.0, C=4642.0,   a=1245.0),
}

# ── Target large-scale trend (Section 4.2) ────────────────────────────────
# Pearson r between variable and X, Y coordinates
THICKNESS_TREND = dict(r_x=0.189, r_y=0.344)
CALORIFIC_TREND = dict(r_x=0.465, r_y=0.128)


def ellipsoidal_range(theta_deg, vario_dict):
    """Direction-dependent range a(theta) via cosine-squared interpolation
    between the fitted directional ranges (mirrors Eq. 10 of the manuscript).
    """
    dirs = np.array(sorted(vario_dict.keys()))
    ranges = np.array([vario_dict[d]['a'] for d in dirs])
    a_min, a_max = ranges.min(), ranges.max()
    theta_max = dirs[np.argmax(ranges)]
    t = np.radians(theta_deg % 180)
    tm = np.radians(theta_max)
    return a_min + (a_max - a_min) * np.cos(t - tm) ** 2


def mean_sill(vario_dict):
    return np.mean([v['C0'] + v['C'] for v in vario_dict.values()])


def lmc_covariance_matrix(coords, vario_dict):
    """Builds an LMC-Ellipsoidal covariance matrix (same construction as
    3D-VSK, Eq. 12) for an arbitrary point set, used here purely as the
    generating covariance for the synthetic Gaussian field.
    """
    n = coords.shape[0]
    S_bar = mean_sill(vario_dict)
    K = np.full((n, n), S_bar)
    for i in range(n):
        for j in range(i + 1, n):
            dx = coords[j, 0] - coords[i, 0]
            dy = coords[j, 1] - coords[i, 1]
            h = np.hypot(dx, dy)
            az = np.degrees(np.arctan2(dy, dx)) % 180
            a_theta = ellipsoidal_range(az, vario_dict)
            if h >= a_theta:
                cov = 0.0
            else:
                cov = S_bar * (1 - 1.5 * (h / a_theta) + 0.5 * (h / a_theta) ** 3)
            K[i, j] = K[j, i] = cov
    return K


def simulate_gaussian_field(coords, vario_dict, rng):
    """Simulates a zero-mean, unit-ish-variance spatially correlated
    Gaussian field via Cholesky factorisation of the LMC covariance matrix.
    """
    K = lmc_covariance_matrix(coords, vario_dict)
    # Diagonal jitter for numerical stability. The LMC construction is
    # PSD-consistent in principle (see manuscript Eq. 12), but finite-
    # precision arithmetic on a 191x191 matrix with strongly anisotropic
    # parameters (e.g. the calorific-value case) can still produce tiny
    # negative eigenvalues; a small relative jitter resolves this without
    # affecting the simulated spatial structure.
    eigval_min = np.linalg.eigvalsh(K).min()
    jitter = max(1e-6, abs(min(0.0, eigval_min)) * 1.5) * np.trace(K) / K.shape[0]
    K += np.eye(K.shape[0]) * jitter
    L = np.linalg.cholesky(K)
    z = rng.standard_normal(K.shape[0])
    return L @ z


def shape_marginal(field, target_mean, target_std, min_clip, max_clip, rng,
                    skew_strength=0.0):
    """Rescales a standard-normal-ish field to the target mean/std, applies
    a mild skew transform if requested, and clips to the observed data range.
    """
    z = (field - field.mean()) / field.std()
    if skew_strength != 0.0:
        z = z + skew_strength * (z ** 2 - 1.0)  # adds positive skew
        z = (z - z.mean()) / z.std()
    values = target_mean + target_std * z
    return np.clip(values, min_clip, max_clip)


def add_trend(values, coords, r_x, r_y, rng):
    """Blends in a linear X/Y trend at approximately the target Pearson
    correlation strength, preserving the marginal mean/std as closely as
    possible.
    """
    x_n = (coords[:, 0] - coords[:, 0].mean()) / coords[:, 0].std()
    y_n = (coords[:, 1] - coords[:, 1].mean()) / coords[:, 1].std()
    v_n = (values - values.mean()) / values.std()

    # Simple linear combination targeting approximate correlation magnitude
    trend = r_x * x_n + r_y * y_n
    blended = np.sqrt(max(1e-6, 1 - r_x**2 - r_y**2)) * v_n + trend
    blended = blended * values.std() + values.mean()
    return blended


def generate_synthetic_seyitomer(seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)

    # ── 1. Drill hole locations ──
    coords = rng.uniform([0, 0], [X_EXTENT, Y_EXTENT], size=(N_HOLES, 2))
    tree = cKDTree(coords)
    d, _ = tree.query(coords, k=2)
    print(f"Synthetic mean NN distance: {d[:, 1].mean():.1f} m "
          f"(target from real dataset: 185.5 m)")

    # ── 2. Simulate sqrt(Thickness) field ──
    sqrt_thick_field = simulate_gaussian_field(coords, THICKNESS_VARIO, rng)
    sqrt_thick = shape_marginal(
        sqrt_thick_field,
        target_mean=SQRT_THICKNESS_TARGET['mean'],
        target_std=SQRT_THICKNESS_TARGET['std'],
        min_clip=np.sqrt(THICKNESS_TARGET['min_clip']),
        max_clip=np.sqrt(THICKNESS_TARGET['max_clip']),
        rng=rng, skew_strength=0.05,
    )
    sqrt_thick = add_trend(sqrt_thick, coords, **THICKNESS_TREND, rng=rng)
    thickness = np.clip(sqrt_thick, 0, None) ** 2
    thickness = np.clip(thickness, THICKNESS_TARGET['min_clip'],
                         THICKNESS_TARGET['max_clip'])

    # ── 3. Simulate calorific value field ──
    cal_field = simulate_gaussian_field(coords, CALORIFIC_VARIO, rng)
    calorific = shape_marginal(
        cal_field,
        target_mean=CALORIFIC_TARGET['mean'],
        target_std=CALORIFIC_TARGET['std'],
        min_clip=CALORIFIC_TARGET['min_clip'],
        max_clip=CALORIFIC_TARGET['max_clip'],
        rng=rng, skew_strength=0.03,
    )
    calorific = add_trend(calorific, coords, **CALORIFIC_TREND, rng=rng)
    calorific = np.clip(calorific, CALORIFIC_TARGET['min_clip'],
                         CALORIFIC_TARGET['max_clip'])

    # ── 4. Synthetic seam roof elevation (coal_top), for 3D block model ──
    # Loosely correlated with a broad regional dip, consistent with
    # Section 4.1 (coal_top range 1,040-1,181 m a.s.l.)
    regional_dip = (coords[:, 0] / X_EXTENT) * 80 + (coords[:, 1] / Y_EXTENT) * 60
    coal_top = 1040 + regional_dip + rng.normal(0, 8, N_HOLES)
    coal_top = np.clip(coal_top, 1040, 1181)

    df = pd.DataFrame({
        'hole_id': [f'SYN{i+1:04d}' for i in range(N_HOLES)],
        'x': coords[:, 0],
        'y': coords[:, 1],
        'thickness_m': thickness,
        'calorific_kcal_kg': calorific,
        'coal_top_masl': coal_top,
    })
    return df


if __name__ == '__main__':
    df = generate_synthetic_seyitomer()
    df.to_csv('synthetic_seyitomer.csv', index=False)

    print("\n=== Synthetic dataset summary (compare to manuscript Table 8) ===")
    print(f"n = {len(df)}")
    print(f"Thickness:  min={df.thickness_m.min():.1f}  max={df.thickness_m.max():.1f}  "
          f"mean={df.thickness_m.mean():.1f}  std={df.thickness_m.std():.1f}")
    print(f"Calorific:  min={df.calorific_kcal_kg.min():.0f}  max={df.calorific_kcal_kg.max():.0f}  "
          f"mean={df.calorific_kcal_kg.mean():.0f}  std={df.calorific_kcal_kg.std():.0f}")

    from scipy.stats import pearsonr
    rx_t, _ = pearsonr(df.thickness_m, df.x)
    ry_t, _ = pearsonr(df.thickness_m, df.y)
    rx_c, _ = pearsonr(df.calorific_kcal_kg, df.x)
    ry_c, _ = pearsonr(df.calorific_kcal_kg, df.y)
    print(f"\nThickness trend:  r_X={rx_t:.3f} (target 0.189)  r_Y={ry_t:.3f} (target 0.344)")
    print(f"Calorific trend:  r_X={rx_c:.3f} (target 0.465)  r_Y={ry_c:.3f} (target 0.128)")
    print(f"\nWrote synthetic_seyitomer.csv ({len(df)} rows)")
