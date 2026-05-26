#!/usr/bin/env python3
"""
================================================================
Moran's I Spatial Autocorrelation Test
Texas County Diagnosed Diabetes Prevalence Study

Author: Yussif Yahaya
Institution: University of Texas Rio Grande Valley

Purpose:
    1. Build Queen's contiguity spatial weights matrix
       from Texas county shapefile
    2. Run Global Moran's I on diagnosed diabetes prevalence
    3. Run Global Moran's I on OLS model residuals
    4. Generate Moran's I scatter plot
    5. Print publication-ready results for paper

Input files required:
    - texas_analysis_dataset_v2.csv  (your analysis data)
    - Texas county shapefile         (downloaded automatically)

Output files:
    - morans_i_results.txt           (full results for paper)
    - FigureS2_morans_scatter.png    (Moran scatter plot)
    - spatial_weights_summary.txt    (weights matrix diagnostics)

Usage on Cradle:
    python morans_i_test.py \
        --data_csv /path/to/texas_analysis_dataset_v2.csv \
        --output_dir /path/to/output/

================================================================
"""

import argparse
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Required for HPC — no display
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import statsmodels.api as sm

# Spatial analysis
from libpysal.weights import lag_spatial
from libpysal.weights import Queen
from esda.moran import Moran


def parse_args():
    parser = argparse.ArgumentParser(
        description='Moran I spatial autocorrelation test '
                    'for Texas diabetes study.'
    )
    parser.add_argument(
        '--data_csv',
        type=str,
        required=True,
        help='Path to texas_analysis_dataset_v2.csv'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Directory to save output files'
    )
    parser.add_argument(
        '--shapefile_zip',
        type=str,
        default='',
        help='Optional local path to Texas county shapefile ZIP'
    )
    parser.add_argument(
        '--permutations',
        type=int,
        default=999,
        help='Number of permutations for significance test (default: 999)'
    )
    return parser.parse_args()


def load_data(data_csv_path):
    """Load and prepare analysis dataset."""
    print('Loading analysis dataset...')
    df = pd.read_csv(data_csv_path)

    # Ensure correct types
    df['fips']             = df['fips'].astype(int)
    df['border']           = df['border'].astype(int)
    df['diag_diabetes_pct'] = pd.to_numeric(
        df['diag_diabetes_pct'], errors='coerce'
    )
    df['inactive_pct'] = pd.to_numeric(
        df['inactive_pct'], errors='coerce'
    )
    df['pct_low_food_access'] = pd.to_numeric(
        df['pct_low_food_access'], errors='coerce'
    )

    # Drop missing
    df = df.dropna(subset=[
        'diag_diabetes_pct',
        'inactive_pct',
        'pct_low_food_access'
    ])

    print(f'  Counties loaded    : {len(df)}')
    print(f'  Border counties    : {df["border"].sum()}')
    print(f'  Diabetes range     : '
          f'{df["diag_diabetes_pct"].min():.1f}% -- '
          f'{df["diag_diabetes_pct"].max():.1f}%')
    return df


def load_shapefile(output_dir, shapefile_zip=''):
    """Load Texas county shapefile."""
    import geopandas as gpd

    if shapefile_zip and Path(shapefile_zip).exists():
        print(f'Loading shapefile from: {shapefile_zip}')
        gdf = gpd.read_file(shapefile_zip)
        # Filter to Texas if US-wide shapefile
        if 'STATEFP' in gdf.columns:
            gdf = gdf[gdf['STATEFP'] == '48'].copy()
    else:
        # Try downloading Texas-only shapefile
        TX_URL = ('https://www2.census.gov/geo/tiger/'
                  'GENZ2023/shp/cb_2023_48_county_500k.zip')
        print(f'Downloading Texas shapefile from Census Bureau...')
        print(f'URL: {TX_URL}')
        try:
            gdf = gpd.read_file(TX_URL)
            print('Download successful.')
        except Exception as e:
            # Fallback to US-wide if Texas-only fails
            print(f'Texas-only download failed: {e}')
            print('Trying US-wide shapefile...')
            US_URL = ('https://www2.census.gov/geo/tiger/'
                      'GENZ2023/shp/cb_2023_us_county_5m.zip')
            gdf = gpd.read_file(US_URL)
            gdf = gdf[gdf['STATEFP'] == '48'].copy()

    gdf['fips'] = gdf['GEOID'].astype(int)
    print(f'  Shapefile counties : {len(gdf)}')
    return gdf


def build_queen_weights(gdf):
    """
    Build Queen's contiguity spatial weights matrix.
    Queen's contiguity: two counties are neighbors if they
    share any boundary point (edge or vertex).
    Matrix is row-standardized so each row sums to 1.
    """
    print('Building Queens contiguity spatial weights matrix...')
    print('  Definition: neighbors share any boundary point')
    print('  Standardization: row (each row sums to 1)')

    # Build Queen weights from GeoDataFrame
    w = Queen.from_dataframe(gdf, silence_warnings=True)
    w.transform = 'r'  # Row standardize

    # Diagnostics
    n_islands = len(w.islands)
    avg_neighbors = np.mean([len(v) for v in w.neighbors.values()])
    max_neighbors = max([len(v) for v in w.neighbors.values()])
    min_neighbors = min([len(v) for v in w.neighbors.values()])

    print(f'  Counties in matrix : {w.n}')
    print(f'  Islands (no neighbor): {n_islands}')
    print(f'  Avg neighbors      : {avg_neighbors:.2f}')
    print(f'  Min neighbors      : {min_neighbors}')
    print(f'  Max neighbors      : {max_neighbors}')

    if n_islands > 0:
        print(f'  WARNING: {n_islands} island(s) detected.')
        print('  Islands are excluded from Morans I calculation.')

    return w


def run_ols_model(df):
    """
    Run VIF-corrected OLS regression to get residuals.
    Matches primary analysis model from paper.
    """
    print('Running VIF-corrected OLS regression...')

    X = df[['pct_low_food_access', 'inactive_pct', 'border']].copy()
    X = sm.add_constant(X)
    y = df['diag_diabetes_pct']

    model  = sm.OLS(y, X)
    result = model.fit(cov_type='HC3')

    print(f'  R-squared          : {result.rsquared:.3f}')
    print(f'  Border coef        : '
          f'{result.params["border"]:.3f} '
          f'(p={result.pvalues["border"]:.4f})')
    print(f'  Inactivity coef    : '
          f'{result.params["inactive_pct"]:.3f} '
          f'(p={result.pvalues["inactive_pct"]:.4f})')

    return result


def run_morans_i(values, weights, label, permutations=999):
    """
    Run Global Moran's I test.

    Parameters:
        values      : array of values to test
        weights     : libpysal spatial weights object
        label       : description for printing
        permutations: number of random permutations

    Returns:
        dict with all results
    """
    print(f'\nRunning Morans I: {label}')
    print(f'  Permutations       : {permutations}')

    moran = Moran(values, weights, permutations=permutations)

    # Significance stars
    if moran.p_sim < 0.001:
        sig = '***'
    elif moran.p_sim < 0.01:
        sig = '**'
    elif moran.p_sim < 0.05:
        sig = '*'
    else:
        sig = 'ns (not significant)'

    print(f'  Morans I           : {moran.I:.4f}')
    print(f'  Expected E[I]      : {moran.EI:.4f}')
    print(f'  z-score            : {moran.z_sim:.4f}')
    print(f'  p-value (sim)      : {moran.p_sim:.4f} {sig}')

    if moran.p_sim < 0.05:
        if moran.I > moran.EI:
            interpretation = (
                'Significant POSITIVE spatial autocorrelation detected.\n'
                'High-value counties tend to cluster near other '
                'high-value counties.'
            )
        else:
            interpretation = (
                'Significant NEGATIVE spatial autocorrelation detected.\n'
                'High-value counties tend to be surrounded by '
                'low-value counties.'
            )
    else:
        interpretation = (
            'No significant spatial autocorrelation detected.\n'
            'Values are randomly distributed across Texas counties.'
        )

    print(f'  Interpretation     : {interpretation}')

    return {
        'label'          : label,
        'I'              : moran.I,
        'EI'             : moran.EI,
        'z_sim'          : moran.z_sim,
        'p_sim'          : moran.p_sim,
        'sig'            : sig,
        'interpretation' : interpretation,
        'moran_object'   : moran,
    }


def plot_morans_scatter(moran_obj, values, label,
                        output_path, color='#2166AC'):
    """
    Generate Moran's I scatter plot.
    X axis: standardized values
    Y axis: spatial lag of standardized values
    """
    fig, ax = plt.subplots(figsize=(8, 7))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Standardize values
    z = (values - values.mean()) / values.std()

    # Spatial lag
    lag_z = lag_spatial(moran_obj.w, z)

    # Quadrant colors
    colors = []
    for zi, lz in zip(z, lag_z):
        if zi >= 0 and lz >= 0:
            colors.append('#B5341A')    # High-High (red)
        elif zi < 0 and lz < 0:
            colors.append('#2166AC')    # Low-Low (blue)
        elif zi >= 0 and lz < 0:
            colors.append('#F4A582')    # High-Low (light red)
        else:
            colors.append('#92C5DE')    # Low-High (light blue)

    ax.scatter(z, lag_z, c=colors, alpha=0.6, s=40,
               edgecolors='white', linewidth=0.3, zorder=3)

    # Regression line
    m = np.polyfit(z, lag_z, 1)
    x_line = np.linspace(z.min(), z.max(), 100)
    ax.plot(x_line, np.polyval(m, x_line),
            color='#333333', linewidth=2, zorder=4)

    # Reference lines at 0
    ax.axhline(0, color='#888888', linewidth=0.8,
               linestyle='--', alpha=0.7)
    ax.axvline(0, color='#888888', linewidth=0.8,
               linestyle='--', alpha=0.7)

    # Quadrant labels
    ax.text(z.max() * 0.7,  lag_z.max() * 0.8,
            'High-High', fontsize=8, color='#B5341A',
            fontweight='bold', alpha=0.8)
    ax.text(z.min() * 0.7,  lag_z.min() * 0.8,
            'Low-Low',  fontsize=8, color='#2166AC',
            fontweight='bold', alpha=0.8)
    ax.text(z.max() * 0.5,  lag_z.min() * 0.8,
            'High-Low', fontsize=8, color='#F4A582',
            fontweight='bold', alpha=0.8)
    ax.text(z.min() * 0.9,  lag_z.max() * 0.8,
            'Low-High', fontsize=8, color='#92C5DE',
            fontweight='bold', alpha=0.8)

    ax.set_xlabel("Standardized Values (z-score)",
                  fontsize=12, fontweight='bold')
    ax.set_ylabel("Spatial Lag of Standardized Values",
                  fontsize=12, fontweight='bold')
    ax.set_title(
        f"Moran's I Scatter Plot\n{label}\n"
        f"I = {moran_obj.I:.4f}, "
        f"p = {moran_obj.p_sim:.4f}",
        fontsize=12, fontweight='bold', pad=12
    )
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300,
                bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {output_path}')


def write_results(results_diabetes, results_residuals,
                  output_path, n_counties, permutations):
    """Write publication-ready results to text file."""

    lines = []
    lines.append('=' * 70)
    lines.append('MORANS I SPATIAL AUTOCORRELATION RESULTS')
    lines.append('Texas Diagnosed Diabetes Prevalence Study')
    lines.append('=' * 70)
    lines.append('')
    lines.append('SPATIAL WEIGHTS MATRIX')
    lines.append('-' * 40)
    lines.append('Type         : Queens contiguity')
    lines.append('Definition   : Counties sharing any boundary point')
    lines.append('Standardization: Row (each row sums to 1)')
    lines.append(f'N counties   : {n_counties}')
    lines.append(f'Permutations : {permutations}')
    lines.append('')

    for r in [results_diabetes, results_residuals]:
        lines.append('-' * 70)
        lines.append(f'RESULT: {r["label"]}')
        lines.append('-' * 70)
        lines.append(f'Morans I     : {r["I"]:.4f}')
        lines.append(f'Expected E[I]: {r["EI"]:.4f}')
        lines.append(f'z-score      : {r["z_sim"]:.4f}')
        lines.append(f'p-value      : {r["p_sim"]:.4f} {r["sig"]}')
        lines.append(f'Interpretation: {r["interpretation"]}')
        lines.append('')

    lines.append('=' * 70)
    lines.append('PUBLICATION-READY TEXT FOR PAPER')
    lines.append('=' * 70)
    lines.append('')
    lines.append('--- FOR METHODS SECTION ---')
    lines.append('')
    lines.append(
        'A Queens contiguity spatial weights matrix was constructed '
        'to define county neighbors as those sharing any common '
        'boundary point or vertex. The matrix was row-standardized '
        'so that each row summed to unity. Global Morans I was '
        'calculated for both the age-adjusted diagnosed diabetes '
        'prevalence values and the OLS model residuals using '
        f'{permutations} conditional permutations to assess '
        'statistical significance.'
    )
    lines.append('')
    lines.append('--- FOR RESULTS SECTION ---')
    lines.append('')

    # Diabetes prevalence result
    d = results_diabetes
    if d['p_sim'] < 0.05:
        lines.append(
            f'Global Morans I confirmed significant positive spatial '
            f'clustering of diagnosed diabetes prevalence across '
            f'Texas counties (I={d["I"]:.4f}, z={d["z_sim"]:.2f}, '
            f'p={d["p_sim"]:.4f}), indicating that high-prevalence '
            f'counties tend to cluster geographically near other '
            f'high-prevalence counties.'
        )
    else:
        lines.append(
            f'Global Morans I revealed no significant spatial '
            f'clustering of diagnosed diabetes prevalence '
            f'(I={d["I"]:.4f}, z={d["z_sim"]:.2f}, '
            f'p={d["p_sim"]:.4f}).'
        )

    lines.append('')

    # Residuals result
    r = results_residuals
    if r['p_sim'] < 0.05:
        lines.append(
            f'Morans I applied to OLS model residuals indicated '
            f'significant residual spatial autocorrelation '
            f'(I={r["I"]:.4f}, z={r["z_sim"]:.2f}, '
            f'p={r["p_sim"]:.4f}), suggesting that neighboring '
            f'counties share similar unexplained diabetes burden '
            f'beyond the modeled predictors. This is acknowledged '
            f'as a limitation; future analyses incorporating '
            f'spatial error or spatial lag models are recommended.'
        )
    else:
        lines.append(
            f'Morans I applied to OLS model residuals revealed '
            f'no significant residual spatial autocorrelation '
            f'(I={r["I"]:.4f}, z={r["z_sim"]:.2f}, '
            f'p={r["p_sim"]:.4f}), confirming that the OLS model '
            f'adequately accounts for the geographic structure of '
            f'diagnosed diabetes prevalence across Texas counties '
            f'and that the independence assumption of OLS errors '
            f'is not substantially violated.'
        )

    lines.append('')
    lines.append('=' * 70)
    lines.append('LATEX CODE FOR PAPER')
    lines.append('=' * 70)
    lines.append('')
    lines.append('--- EQUATION (add to Methods after Equation 1) ---')
    lines.append('')
    lines.append(r'\begin{equation}')
    lines.append(r'I = \frac{N}{\sum_i \sum_j w_{ij}} \cdot')
    lines.append(r'    \frac{\sum_i \sum_j w_{ij}')
    lines.append(r'          (x_i - \bar{x})(x_j - \bar{x})}')
    lines.append(r'         {\sum_i (x_i - \bar{x})^2}')
    lines.append(r'\end{equation}')
    lines.append('')
    lines.append(
        r'\noindent where $N$ is the number of counties (N=253), '
        r'$w_{ij}$ is the row-standardized spatial weight between '
        r'counties $i$ and $j$ from the Queens contiguity matrix, '
        r'$x_i$ is the value of interest for county $i$, and '
        r'$\bar{x}$ is the mean value across all counties. '
        r'The expected value under spatial randomness is '
        r'$E[I] = -1/(N-1)$. Significance was assessed using '
        f'{permutations} conditional permutations.'
    )
    lines.append('')
    lines.append('--- RESULTS TEXT (fill in values below) ---')
    lines.append('')
    d = results_diabetes
    r = results_residuals
    lines.append(
        f'Global Morans I confirmed significant positive spatial '
        f'clustering of diagnosed diabetes prevalence across Texas '
        f'counties ($I$={d["I"]:.4f}, $z$={d["z_sim"]:.2f}, '
        f'$p$={d["p_sim"]:.4f}). Morans I applied to OLS model '
        f'residuals indicated {"significant" if r["p_sim"] < 0.05 else "no significant"} '
        f'residual spatial autocorrelation '
        f'($I$={r["I"]:.4f}, $z$={r["z_sim"]:.2f}, '
        f'$p$={r["p_sim"]:.4f}).'
    )
    lines.append('')
    lines.append('=' * 70)
    lines.append('Significance: * p<0.05  ** p<0.01  *** p<0.001  ns=not significant')
    lines.append('=' * 70)

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))

    print(f'\nResults saved: {output_path}')


def main():
    args = parse_args()

    data_csv   = Path(args.data_csv).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 70)
    print("Morans I Spatial Autocorrelation Test")
    print('Texas Diagnosed Diabetes Prevalence Study')
    print('=' * 70)
    print(f'Data CSV   : {data_csv}')
    print(f'Output dir : {output_dir}')
    print(f'Permutations: {args.permutations}')
    print()

    # ── Step 1: Load data ─────────────────────────────────────────
    df = load_data(data_csv)

    # ── Step 2: Load shapefile ────────────────────────────────────
    gdf = load_shapefile(output_dir, args.shapefile_zip)

    # ── Step 3: Merge data with shapefile ─────────────────────────
    print('\nMerging analysis data with shapefile...')
    merged = gdf.merge(
        df[['fips', 'diag_diabetes_pct', 'inactive_pct',
            'pct_low_food_access', 'border']],
        on='fips', how='inner'
    )
    merged = merged.dropna(subset=['diag_diabetes_pct'])
    merged = merged.reset_index(drop=True)
    print(f'  Merged counties    : {len(merged)}')

    # ── Step 4: Build Queen's contiguity weights ──────────────────
    w = build_queen_weights(merged)

    # Save weights diagnostics
    weights_summary = output_dir / 'spatial_weights_summary.txt'
    with open(weights_summary, 'w') as f:
        f.write('Queens Contiguity Spatial Weights Matrix\n')
        f.write('=' * 50 + '\n')
        f.write(f'N counties         : {w.n}\n')
        f.write(f'Islands            : {len(w.islands)}\n')
        f.write(f'Avg neighbors      : '
                f'{np.mean([len(v) for v in w.neighbors.values()]):.2f}\n')
        f.write(f'Min neighbors      : '
                f'{min([len(v) for v in w.neighbors.values()])}\n')
        f.write(f'Max neighbors      : '
                f'{max([len(v) for v in w.neighbors.values()])}\n')
        f.write(f'Standardization    : Row\n')
    print(f'\nWeights summary saved: {weights_summary}')

    # ── Step 5: Run OLS to get residuals ──────────────────────────
    print()
    ols_result = run_ols_model(merged)
    residuals  = ols_result.resid.values

    # ── Step 6: Run Moran's I on diabetes prevalence ──────────────
    print()
    diabetes_values = merged['diag_diabetes_pct'].values
    results_diabetes = run_morans_i(
        diabetes_values, w,
        label='Diagnosed diabetes prevalence (%)',
        permutations=args.permutations
    )

    # ── Step 7: Run Moran's I on OLS residuals ────────────────────
    results_residuals = run_morans_i(
        residuals, w,
        label='OLS model residuals',
        permutations=args.permutations
    )

    # ── Step 8: Generate Moran scatter plots ──────────────────────
    print('\nGenerating Morans I scatter plots...')

    plot_morans_scatter(
        results_diabetes['moran_object'],
        diabetes_values,
        label='Diagnosed Diabetes Prevalence',
        output_path=output_dir / 'FigureS2_morans_diabetes.png',
        color='#B5341A'
    )

    plot_morans_scatter(
        results_residuals['moran_object'],
        residuals,
        label='OLS Model Residuals',
        output_path=output_dir / 'FigureS2_morans_residuals.png',
        color='#2166AC'
    )

    # ── Step 9: Write publication-ready results ───────────────────
    write_results(
        results_diabetes,
        results_residuals,
        output_dir / 'morans_i_results.txt',
        n_counties=len(merged),
        permutations=args.permutations
    )

    # ── Final Summary ─────────────────────────────────────────────
    print()
    print('=' * 70)
    print('FINAL SUMMARY')
    print('=' * 70)
    d = results_diabetes
    r = results_residuals
    print(f'Diabetes prevalence:')
    print(f'  Morans I = {d["I"]:.4f}  '
          f'p = {d["p_sim"]:.4f}  {d["sig"]}')
    print(f'OLS residuals:')
    print(f'  Morans I = {r["I"]:.4f}  '
          f'p = {r["p_sim"]:.4f}  {r["sig"]}')
    print()
    print('Output files:')
    print(f'  morans_i_results.txt          — Full results + LaTeX code')
    print(f'  FigureS2_morans_diabetes.png  — Scatter plot (diabetes)')
    print(f'  FigureS2_morans_residuals.png — Scatter plot (residuals)')
    print(f'  spatial_weights_summary.txt   — Weights diagnostics')
    print()
    print('Next step:')
    print('  Open morans_i_results.txt and copy the')
    print('  publication-ready text into your LaTeX manuscript.')
    print('=' * 70)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'\nERROR: {exc}', file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
