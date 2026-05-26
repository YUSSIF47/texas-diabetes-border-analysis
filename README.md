# Texas Border Diabetes Analysis

## Overview

Analysis code and dataset for:

> Saha, P.R., Khan, S., Yahaya, Y., & Meia, M.A.A. (2025). Border-Region Status and Diagnosed Diabetes Prevalence in Texas: A Cross-Sectional Ecological Analysis. *PLOS ONE* (under review).

This study examined whether Texas–Mexico border-region county status is independently associated with diagnosed diabetes prevalence across 253 Texas counties, after controlling for physical inactivity and low food access.

---

## Key Findings

- Border-region counties had **33% higher** unadjusted mean diagnosed diabetes prevalence (16.1% vs. 12.1%)
- Border-region status remained significant after adjustment: β = 0.625 (95% CI [0.357, 0.893], p < 0.001)
- Physical inactivity was the strongest independent predictor (β = 0.404, p < 0.001)
- Global Moran's I confirmed strong spatial clustering (I = 0.5734, p = 0.001), reduced in OLS residuals (I = 0.1696, p = 0.001)
- Model R² = 0.960, N = 253 counties

---

## Repository Structure

```
texas-diabetes-border-analysis/
│
├── Texas_Diabetes_Analysis_v2.ipynb          ← main analysis notebook
├── requirements.txt                          ← Python dependencies
├── README.md
│
├── scripts/
│   ├── figure6_choropleth_publication_FINAL_v6.py  ← Figure 6 choropleth map
│   ├── morans_i_test.py                            ← Moran's I spatial analysis
│   └── slurm_morans_i_test.sh                      ← HPC job script (UTRGV CRADLE)
│
├── data/
│   ├── texas_analysis_dataset_v2.csv         ← analytical dataset (253 counties)
│   ├── PLACES_County_2025.csv                ← download separately (see below)
│   └── FoodEnvironmentAtlas_2025.csv         ← download separately (see below)
│
└── outputs/                                  ← figures and tables (auto-created)
```

---

## Data Downloads

The analytical dataset (`texas_analysis_dataset_v2.csv`) is included in this repository.
The two raw federal source files are **not included** due to size — download them directly:

| File | Source | URL |
|------|--------|-----|
| `PLACES_County_2025.csv` | CDC PLACES 2025 | https://data.cdc.gov (search: "PLACES County Data 2025 release") |
| `FoodEnvironmentAtlas_2025.csv` | USDA Food Environment Atlas 2025 | https://www.ers.usda.gov/data-products/food-environment-atlas |
| `cb_2023_us_county_5m.zip` | US Census Bureau shapefile | https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_county_5m.zip |

After downloading, rename files to match the names above and place them in the `data/` folder.
The choropleth and Moran's I scripts will also auto-download the Census shapefile if not provided.

---

## How to Run

### 1. Main Analysis Notebook

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/texas-diabetes-border-analysis.git
cd texas-diabetes-border-analysis

# Install dependencies
pip install -r requirements.txt

# Place raw data files in data/ (see Data Downloads above)

# Run the notebook top-to-bottom
jupyter notebook Texas_Diabetes_Analysis_v2.ipynb
```

All figures and result tables are saved to `outputs/`.

---

### 2. Figure 6 — Choropleth Map

Generates the publication-ready Texas county choropleth map (Figure 6) with La Paz border-region county outlines.

```bash
python scripts/figure6_choropleth_publication_FINAL_v6.py \
    --data_csv data/texas_analysis_dataset_v2.csv \
    --shapefile_zip data/cb_2023_us_county_5m.zip \
    --output_dir outputs/ \
    --dpi 600
```

**Key arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--data_csv` | required | Path to analysis CSV |
| `--shapefile_zip` | auto-download | Census county shapefile ZIP |
| `--output_dir` | required | Output folder |
| `--figure_name` | `Figure6_choropleth_map_publication_FINAL_v6.png` | Output filename |
| `--dpi` | 600 | Export resolution |
| `--vmin` / `--vmax` | 8.0 / 22.0 | Color scale range |
| `--label_counties` | `Dimmit,Zavala,Starr,Jim Hogg` | Counties to label |
| `--save_tiff` | off | Also save TIFF file |

**Outputs:** PNG, PDF (and optionally TIFF) + `Figure6_merge_quality_check.csv`

---

### 3. Moran's I Spatial Autocorrelation Test

Builds a Queen's contiguity spatial weights matrix and runs Global Moran's I on both
diagnosed diabetes prevalence and OLS model residuals.

```bash
python scripts/morans_i_test.py \
    --data_csv data/texas_analysis_dataset_v2.csv \
    --shapefile_zip data/cb_2023_us_county_5m.zip \
    --output_dir outputs/morans_i/ \
    --permutations 999
```

**Outputs:**

| File | Description |
|------|-------------|
| `morans_i_results.txt` | Full results + publication-ready text + LaTeX code |
| `FigureS2_morans_diabetes.png` | Moran scatter plot — diabetes prevalence |
| `FigureS2_morans_residuals.png` | Moran scatter plot — OLS residuals |
| `spatial_weights_summary.txt` | Weights matrix diagnostics |

---

### 4. Running on HPC (UTRGV CRADLE)

The Moran's I test was run on the UTRGV CRADLE high-performance computing cluster
using the provided SLURM job script:

```bash
# Submit job
sbatch scripts/slurm_morans_i_test.sh
```

Edit the `PROJECT_DIR` variable in the script to match your CRADLE home directory
before submitting. The script activates the `geofig` conda environment and runs
`morans_i_test.py` with 999 permutations.

> **Note for CRADLE users:** The script uses the `normal` CPU partition.
> Moran's I does not require GPU resources.

---

## Computational Environment

| Component | Version |
|-----------|---------|
| Python | 3.9 |
| pandas | ≥ 1.5 |
| numpy | ≥ 1.23 |
| matplotlib | ≥ 3.6 |
| seaborn | ≥ 0.12 |
| scipy | ≥ 1.9 |
| statsmodels | 0.14.6 |
| libpysal | 4.13.0 |
| esda | 2.7.0 |
| geopandas | ≥ 0.12 |

Primary regression analyses were run locally (Python 3.9, Jupyter Notebook).
Spatial autocorrelation diagnostics were run on the UTRGV CRADLE HPC cluster.

---

## Border County Definition

Border-region counties were defined using the **official La Paz Agreement 32-county
definition** — counties within 100 km (62.1 miles) of the US–Mexico international
boundary, as recognized by the US–Mexico Border Health Commission (2010).

The 32 counties are:
El Paso, Hudspeth, Culberson, Jeff Davis, Presidio, Brewster, Terrell, Val Verde,
Kinney, Maverick, Webb, Zapata, Starr, Hidalgo, Cameron, Willacy, Brooks, Jim Hogg,
Dimmit, Zavala, Uvalde, Edwards, Real, Bandera, Medina, Frio, La Salle, McMullen,
Duval, Jim Wells, Kleberg, Kenedy.

A sensitivity analysis using a narrower 20-county direct-boundary definition is
reported in Supplementary Table S2.

---

## License

This code is released under the [MIT License](LICENSE).
Data are from federal public domain sources (CDC and USDA).

---

## Citation
If you use this code or data, please cite:
> Saha, P.R., Khan, S., Yahaya, Y., & Meia, M.A.A. (2025). Border-Region Status
> and Diagnosed Diabetes Prevalence in Texas: A Cross-Sectional Ecological Analysis.
> *PLOS ONE* (under review). Code and data: https://doi.org/10.5281/zenodo.20390172

---

## Contact

**Yussif Yahaya** (Corresponding Author)  
PhD Student, Department of Mathematics and Statistics  
University of Texas Rio Grande Valley  
yussif.yahaya01@utrgv.edu
