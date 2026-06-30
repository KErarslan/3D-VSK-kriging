# Data

## `synthetic_seyitomer_generator.py`

Generates `synthetic_seyitomer.csv`: a statistically-equivalent SYNTHETIC
replacement for the real Seyitömer-Aslanlı drill hole dataset used in the
manuscript. See the script's docstring for full methodology and rationale.

Run it directly to regenerate the CSV (deterministic, fixed random seed):

```bash
python synthetic_seyitomer_generator.py
```

## `synthetic_seyitomer.csv`

A pre-generated copy of the synthetic dataset (191 rows), included here for
convenience so the notebook can be run without first executing the generator
separately. Columns:

| Column | Description |
|---|---|
| `hole_id` | Synthetic drill hole identifier |
| `x`, `y` | Synthetic horizontal coordinates (m) |
| `thickness_m` | Synthetic seam thickness (m) |
| `calorific_kcal_kg` | Synthetic calorific value (kcal/kg) |
| `coal_top_masl` | Synthetic seam roof elevation (m above sea level) |

**This is not real survey data.** Absolute values, hole identities, and exact
spatial positions carry no geological meaning — see the main `README.md`
("Data Availability") for why this synthetic dataset exists and what it does
and does not reproduce from the published manuscript.

## Kalgoorlie dataset

The Kalgoorlie Au dataset is not provided as a separate file: it is included
directly (hardcoded) inside `notebooks/3DSurfaceKriging.ipynb`, Cell 1. See
the main `README.md` for provenance and modifications applied.
