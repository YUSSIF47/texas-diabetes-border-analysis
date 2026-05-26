#!/usr/bin/env python3
"""
Publication-ready Figure 6 choropleth map — FINAL_v6
Texas county-level age-adjusted diagnosed diabetes prevalence
with La Paz Agreement border-region counties outlined.

Main fixes included:
  1. County callout labels are placed OUTSIDE the Texas polygon.
  2. Labels are no longer grouped in one block and are not allowed to touch the map.
  3. Dimmit/Zavala and Jim Hogg/Starr are separated vertically in southwest white space.
  4. Short elbow connectors are used instead of long diagonal leader lines.
  5. Only small county markers remain inside the county polygons.
  6. Border-region county outlines remain visible, including counties with suppressed data.
  7. Legend is placed below the map, not over any counties.
  8. PNG and PDF are saved automatically; TIFF is optional.

Example use on UTRGV Cradle:

python figure6_choropleth_publication_FINAL_v6.py \
    --data_csv /home/USER/texas_diabetes_map/data/texas_analysis_dataset_v3.csv \
    --shapefile_zip /home/USER/texas_diabetes_map/data/cb_2023_us_county_5m.zip \
    --output_dir /home/USER/texas_diabetes_map/results \
    --figure_name Figure6_choropleth_map_publication_FINAL_v6.png

If you do not provide --shapefile_zip, the script first checks the output folder
for cb_2023_us_county_5m.zip and then tries to download it from the U.S. Census.
On HPC systems without internet access, provide --shapefile_zip directly.
"""

import argparse
import sys
from pathlib import Path
from urllib.request import urlretrieve

import matplotlib
matplotlib.use("Agg")

import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Point


DEFAULT_SHAPEFILE_URL = (
    "https://www2.census.gov/geo/tiger/"
    "GENZ2023/shp/cb_2023_us_county_5m.zip"
)


# -----------------------------------------------------------------------------
# Arguments and input utilities
# -----------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create a publication-ready Texas county choropleth map for "
            "age-adjusted diagnosed diabetes prevalence."
        )
    )
    parser.add_argument("--data_csv", required=True, help="Path to analysis CSV file.")
    parser.add_argument(
        "--shapefile_zip",
        default="",
        help="Path to Census county shapefile ZIP, e.g., cb_2023_us_county_5m.zip.",
    )
    parser.add_argument("--output_dir", required=True, help="Directory for saved outputs.")
    parser.add_argument(
        "--figure_name",
        default="Figure6_choropleth_map_publication_FINAL_v6.png",
        help="Output PNG filename.",
    )
    parser.add_argument("--dpi", type=int, default=600, help="PNG/TIFF export DPI.")

    # Data column options. Defaults match your current dataset.
    parser.add_argument("--fips_col", default="fips", help="FIPS column in CSV.")
    parser.add_argument("--county_col", default="county", help="County name column in CSV.")
    parser.add_argument(
        "--value_col",
        default="diag_diabetes_pct",
        help="Diabetes prevalence percentage column in CSV.",
    )
    parser.add_argument("--border_col", default="border", help="Border-region indicator column.")

    # Visual controls.
    parser.add_argument("--vmin", type=float, default=8.0, help="Color scale minimum.")
    parser.add_argument("--vmax", type=float, default=22.0, help="Color scale maximum.")
    parser.add_argument("--us_avg", type=float, default=12.7, help="U.S. average marker.")
    parser.add_argument(
        "--label_counties",
        default="Dimmit,Zavala,Starr,Jim Hogg",
        help="Comma-separated counties to label outside the map.",
    )
    parser.add_argument(
        "--save_tiff",
        action="store_true",
        help="Also save a high-resolution TIFF file.",
    )
    parser.add_argument(
        "--title",
        default="Diagnosed Diabetes Prevalence Across Texas Counties",
        help="Main figure title.",
    )
    parser.add_argument(
        "--subtitle",
        default="La Paz Border-Region Counties Outlined in Navy Blue",
        help="Figure subtitle.",
    )
    return parser.parse_args()


def normalize_fips(series: pd.Series) -> pd.Series:
    """Return 5-character FIPS strings, robust to int/float/string CSV values."""
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
        .str.extract(r"(\d+)", expand=False)
        .fillna("")
        .str.zfill(5)
    )


def require_columns(df: pd.DataFrame, required_cols):
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError("Missing required columns in data CSV: " + ", ".join(missing))


def get_shapefile_path(output_dir: Path, shapefile_zip: str) -> Path:
    if shapefile_zip:
        path = Path(shapefile_zip).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Provided shapefile does not exist: {path}")
        return path

    local_zip = output_dir / "cb_2023_us_county_5m.zip"
    if local_zip.exists():
        return local_zip

    print("No --shapefile_zip provided and local ZIP not found.")
    print("Trying to download Census county shapefile...")
    try:
        urlretrieve(DEFAULT_SHAPEFILE_URL, local_zip)
        return local_zip
    except Exception as exc:
        raise RuntimeError(
            "Could not download the shapefile. On Cradle/HPC, download it manually:\n"
            "wget https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_county_5m.zip\n"
            "Then rerun with:\n"
            "--shapefile_zip /path/to/cb_2023_us_county_5m.zip"
        ) from exc


def load_and_merge(args):
    data_csv = Path(args.data_csv).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not data_csv.exists():
        raise FileNotFoundError(f"Data CSV not found: {data_csv}")

    print("=" * 78)
    print("Figure 6 — Publication-ready Texas diabetes choropleth map FINAL_v6")
    print("=" * 78)
    print(f"Data CSV      : {data_csv}")
    print(f"Output folder : {output_dir}")

    df = pd.read_csv(data_csv)
    required = [args.fips_col, args.county_col, args.value_col, args.border_col]
    require_columns(df, required)

    df = df.copy()
    df["_fips5"] = normalize_fips(df[args.fips_col])
    df[args.value_col] = pd.to_numeric(df[args.value_col], errors="coerce")
    df[args.border_col] = pd.to_numeric(df[args.border_col], errors="coerce").fillna(0).astype(int)
    df["_county_csv"] = df[args.county_col].astype(str).str.replace(" County", "", regex=False).str.strip()

    shp = get_shapefile_path(output_dir, args.shapefile_zip)
    print(f"Shapefile     : {shp}")

    counties = gpd.read_file(str(shp))
    if "STATEFP" in counties.columns:
        tx = counties[counties["STATEFP"].astype(str) == "48"].copy()
    elif "GEOID" in counties.columns:
        tx = counties[counties["GEOID"].astype(str).str.startswith("48")].copy()
    else:
        raise ValueError("Shapefile must contain STATEFP or GEOID column.")

    if tx.empty:
        raise ValueError("No Texas counties found in shapefile. Check shapefile source.")

    tx["_fips5"] = tx["GEOID"].astype(str).str.zfill(5)
    tx["_county_shape"] = tx.get("NAME", "").astype(str).str.strip()

    merged = tx.merge(
        df[["_fips5", "_county_csv", args.value_col, args.border_col]],
        on="_fips5",
        how="left",
    )

    merged["county_label"] = merged["_county_csv"].where(
        merged["_county_csv"].notna() & (merged["_county_csv"].astype(str).str.len() > 0),
        merged["_county_shape"],
    )
    merged[args.border_col] = merged[args.border_col].fillna(0).astype(int)
    merged[args.value_col] = pd.to_numeric(merged[args.value_col], errors="coerce")

    # Projection requested in the figure footer.
    merged = merged.to_crs("EPSG:3083")

    # Basic QC summary.
    n_total = len(merged)
    n_with_data = int(merged[args.value_col].notna().sum())
    n_no_data = int(merged[args.value_col].isna().sum())
    n_border = int((merged[args.border_col] == 1).sum())

    print("-" * 78)
    print(f"Texas counties in shapefile : {n_total}")
    print(f"Counties with data          : {n_with_data}")
    print(f"Data suppressed/unavailable : {n_no_data}")
    print(f"La Paz border counties      : {n_border}")
    if n_total != 254:
        print("WARNING: Texas should have 254 counties. Check shapefile/merge.")
    if n_border != 32:
        print("WARNING: Expected 32 La Paz border-region counties. Check border column.")

    # Save merge QC file to identify any county that is missing/suppressed.
    qc = pd.DataFrame(
        {
            "fips": merged["_fips5"],
            "county": merged["county_label"],
            "diabetes_pct": merged[args.value_col],
            "border": merged[args.border_col],
            "has_diabetes_data": merged[args.value_col].notna(),
        }
    ).sort_values(["border", "county"], ascending=[False, True])
    qc_path = output_dir / "Figure6_merge_quality_check.csv"
    qc.to_csv(qc_path, index=False)
    print(f"Merge QC CSV saved          : {qc_path}")
    print("=" * 78)

    return merged, output_dir


# -----------------------------------------------------------------------------
# Label placement utilities
# -----------------------------------------------------------------------------

def _state_union(map_data: gpd.GeoDataFrame):
    """Return a single Texas polygon/multipolygon. Works across Shapely versions."""
    try:
        return map_data.geometry.union_all()
    except AttributeError:
        return map_data.geometry.unary_union


def find_southwest_label_point(row, state_poly, dx, dy, county_key):
    """
    FINAL_v6 placement rule:
      - Put the label anchor in the southwest white space, outside Texas.
      - Move the right edge of the label box far enough from the polygon so the
        box does not touch/cover South Texas.
      - Stagger labels so Dimmit/Zavala and Jim Hogg/Starr do not overlap.

    Important: for southwest labels we use ha="right" later, so the data point
    returned here is the RIGHT EDGE/CENTER of the label box. The full label box
    extends to the LEFT of this x-coordinate. Therefore, this x-coordinate must
    remain outside the map polygon with a generous westward margin.
    """
    c = row.geometry.centroid

    # FINAL_v6: force labels clearly into the southwest/Mexico-side white space.
    # The earlier v5 offsets were still too close to the county polygons after
    # the final axis scaling. These stronger offsets keep every label box away
    # from the Texas polygon while keeping arrows much shorter than the right-
    # side layout. Values are fractions of full Texas width/height in EPSG:3083.
    offsets = {
        "dimmit":   (-0.335 * dx, -0.065 * dy),
        "zavala":   (-0.335 * dx, -0.112 * dy),
        "jim hogg": (-0.270 * dx, -0.070 * dy),
        "starr":    (-0.270 * dx, -0.122 * dy),
    }
    ox, oy = offsets.get(county_key, (-0.300 * dx, -0.095 * dy))
    candidate = Point(c.x + ox, c.y + oy)

    # If the anchor still falls inside Texas, push it farther west/southwest.
    # This is a safeguard for future data/shapefile variations.
    push_vectors = {
        "dimmit":   (-1.00, -0.18),
        "zavala":   (-1.00, -0.22),
        "jim hogg": (-0.94, -0.34),
        "starr":    (-0.88, -0.42),
    }
    vx, vy = push_vectors.get(county_key, (-1.0, -0.25))
    length = (vx * vx + vy * vy) ** 0.5
    ux, uy = vx / length, vy / length

    step = 0.012 * max(dx, dy)
    retry = 0
    while state_poly.contains(candidate) and retry < 70:
        candidate = Point(candidate.x + ux * step, candidate.y + uy * step)
        retry += 1

    # Extra padding after the point is outside. This makes the right edge of each
    # label box stay clear of the Texas border even after PDF cropping.
    extra_pad = {
        "dimmit":   (-0.055 * dx,  0.000 * dy),
        "zavala":   (-0.055 * dx, -0.006 * dy),
        "jim hogg": (-0.046 * dx,  0.000 * dy),
        "starr":    (-0.046 * dx, -0.006 * dy),
    }
    px, py = extra_pad.get(county_key, (-0.050 * dx, 0.0))
    candidate = Point(candidate.x + px, candidate.y + py)

    return candidate

def add_nearby_callout(ax, row, text_x, text_y, value_col, dx, ha=None, box_edge="#8B2B20"):
    """
    Add a clean southwest-side callout:
      - small marker on the county centroid,
      - short elbow leader line,
      - label box outside Texas and separated from the other boxes.
    """
    centroid = row.geometry.centroid
    county = str(row["county_label"])
    value = row[value_col]

    ax.plot(
        centroid.x,
        centroid.y,
        marker="o",
        markersize=4.2,
        markerfacecolor="#2B2B2B",
        markeredgecolor="white",
        markeredgewidth=0.75,
        linestyle="None",
        zorder=9,
        clip_on=False,
    )

    if ha is None:
        ha = "right" if text_x < centroid.x else "left"

    # End the connector at the edge of the label box.
    label_gap = 0.010 * dx
    line_end_x = text_x + label_gap if ha == "right" else text_x - label_gap

    # Compact elbow connector. Keep the elbow close to the label so the line is
    # visible but does not run across the interior of the map.
    elbow_x = line_end_x + 0.014 * dx if ha == "right" else line_end_x - 0.014 * dx
    elbow_y = text_y

    ax.plot(
        [centroid.x, elbow_x, line_end_x],
        [centroid.y, elbow_y, elbow_y],
        color="#4A4A4A",
        linewidth=0.65,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=8,
        clip_on=False,
    )

    label = f"{county}: {value:.1f}%"
    ax.text(
        text_x,
        text_y,
        label,
        ha=ha,
        va="center",
        fontsize=8.5,
        fontweight="bold",
        color="#1A1A1A",
        bbox=dict(
            boxstyle="round,pad=0.23",
            facecolor="white",
            edgecolor=box_edge,
            linewidth=0.75,
            alpha=0.98,
        ),
        zorder=10,
        clip_on=False,
    )


# -----------------------------------------------------------------------------
# Plotting function
# -----------------------------------------------------------------------------

def make_figure(map_data: gpd.GeoDataFrame, output_dir: Path, args):
    value_col = args.value_col
    border_col = args.border_col

    output_png = output_dir / args.figure_name
    output_pdf = output_png.with_suffix(".pdf")
    output_tiff = output_png.with_suffix(".tif")

    cmap = plt.cm.YlOrRd
    norm = mcolors.Normalize(vmin=args.vmin, vmax=args.vmax)
    navy = "#001F8F"

    # Figure size balances map area, right callout strip, colorbar, legend, and footer.
    fig = plt.figure(figsize=(12.8, 8.5), facecolor="white")

    # Main map axes. Right side inside this axes is reserved for outside labels.
    ax = fig.add_axes([0.025, 0.155, 0.805, 0.755])
    ax.set_facecolor("white")

    # Choropleth fill.
    map_data.plot(
        column=value_col,
        ax=ax,
        cmap=cmap,
        norm=norm,
        edgecolor="#F7F7F7",
        linewidth=0.28,
        legend=False,
        missing_kwds={
            "color": "#D9D9D9",
            "edgecolor": "#AFAFAF",
            "linewidth": 0.35,
        },
        zorder=1,
    )

    # County boundary overlay for readability.
    map_data.boundary.plot(ax=ax, color="#C8C8C8", linewidth=0.22, zorder=2)

    # State outer boundary.
    try:
        map_data.dissolve().boundary.plot(ax=ax, color="#666666", linewidth=0.55, zorder=3)
    except Exception:
        pass

    # Border-region outline. Slightly thinner than v6 for cleaner print output.
    border_geo = map_data[map_data[border_col] == 1]
    border_geo.boundary.plot(ax=ax, color=navy, linewidth=1.50, zorder=6)

    # Compute map extent and dynamic label positions before setting limits.
    xmin, ymin, xmax, ymax = map_data.total_bounds
    dx = xmax - xmin
    dy = ymax - ymin

    requested = [c.strip() for c in args.label_counties.split(",") if c.strip()]
    selected_rows = []
    for cname in requested:
        rows = map_data[
            map_data["county_label"].astype(str).str.strip().str.lower() == cname.lower()
        ]
        if rows.empty:
            print(f"WARNING: Label county not found: {cname}")
            continue
        row = rows.iloc[0]
        if pd.isna(row[value_col]):
            print(f"WARNING: Label county has missing diabetes value and was skipped: {cname}")
            continue
        selected_rows.append(row)

    texas_poly = _state_union(map_data)
    label_positions = []
    for row in selected_rows:
        cname_key = str(row["county_label"]).strip().lower()
        p = find_southwest_label_point(row, texas_poly, dx, dy, cname_key)
        centroid = row.geometry.centroid
        ha = "right" if p.x < centroid.x else "left"
        label_positions.append((row, p.x, p.y, ha))

    # Axis extent includes nearby outside labels but avoids the excessive right-side whitespace.
    label_xs = [x for _, x, _, _ in label_positions]
    label_ys = [y for _, _, y, _ in label_positions]
    plot_xmin = min([xmin] + label_xs) - 0.030 * dx
    plot_xmax = max([xmax] + label_xs) + 0.055 * dx
    plot_ymin = min([ymin] + label_ys) - 0.070 * dy
    plot_ymax = ymax + 0.035 * dy
    ax.set_xlim(plot_xmin, plot_xmax)
    ax.set_ylim(plot_ymin, plot_ymax)
    ax.set_axis_off()

    # Nearby outside callout labels: separated and closer to their counties.
    for row, text_x, text_y, ha in label_positions:
        cname = str(row["county_label"]).strip().lower()
        edge = "#7B0000" if cname == "dimmit" else "#8B2B20"
        add_nearby_callout(
            ax,
            row,
            text_x=text_x,
            text_y=text_y,
            value_col=value_col,
            dx=dx,
            ha=ha,
            box_edge=edge,
        )

    # Colorbar with compact two-line title.
    cax = fig.add_axes([0.860, 0.315, 0.022, 0.465])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.ax.tick_params(labelsize=9.4, width=0.6, length=3)
    cbar.outline.set_linewidth(0.7)
    cbar.ax.set_title(
        "Age-adjusted diabetes\nprevalence (%)",
        fontsize=10.2,
        fontweight="bold",
        pad=10,
    )

    # U.S. average line. Use colorbar data coordinates, not normalized axes coords.
    if args.vmin <= args.us_avg <= args.vmax:
        cbar.ax.axhline(args.us_avg, color="#0A7A2A", linewidth=1.45, linestyle="--")
        cbar.ax.text(
            1.35,
            args.us_avg,
            f"US avg\n{args.us_avg:.1f}%",
            transform=cbar.ax.get_yaxis_transform(),
            fontsize=7.7,
            color="#0A7A2A",
            va="center",
            ha="left",
            fontweight="bold",
        )

    # Title and subtitle.
    fig.text(
        0.50,
        0.964,
        args.title,
        ha="center",
        va="top",
        fontsize=15.4,
        fontweight="bold",
        color="#111111",
    )
    fig.text(
        0.50,
        0.930,
        args.subtitle,
        ha="center",
        va="top",
        fontsize=12.8,
        fontweight="bold",
        color="#111111",
    )

    # Legend below the map, never over counties.
    n_border = int((map_data[border_col] == 1).sum())
    n_no_data = int(map_data[value_col].isna().sum())
    border_handle = mpatches.Patch(
        facecolor="white",
        edgecolor=navy,
        linewidth=1.50,
        label=f"La Paz border-region counties (n={n_border})",
    )
    nodata_handle = mpatches.Patch(
        facecolor="#D9D9D9",
        edgecolor="#AFAFAF",
        linewidth=0.8,
        label=f"Data suppressed/unavailable (n={n_no_data})",
    )
    legend = fig.legend(
        handles=[border_handle, nodata_handle],
        loc="lower center",
        bbox_to_anchor=(0.420, 0.070),
        bbox_transform=fig.transFigure,
        ncol=2,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="#CCCCCC",
        fontsize=9.5,
        handlelength=2.2,
        handleheight=1.05,
        borderpad=0.65,
        columnspacing=1.8,
        handletextpad=0.55,
    )
    legend.get_frame().set_linewidth(0.75)

    # Footer. Use a typographic en dash between US and Mexico for publication.
    footer = (
        "Data: CDC PLACES 2025   |   "
        "Border definition: La Paz Agreement — 32 counties within 100 km of the US–Mexico boundary   |   "
        "Projection: Texas Centric Albers Equal Area (EPSG:3083)"
    )
    fig.text(
        0.016,
        0.023,
        footer,
        ha="left",
        va="bottom",
        fontsize=7.8,
        color="#4D4D4D",
        style="italic",
    )

    # Save outputs.
    plt.savefig(output_png, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    if args.save_tiff:
        plt.savefig(output_tiff, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("=" * 78)
    print(f"Saved PNG : {output_png}")
    print(f"Saved PDF : {output_pdf}")
    if args.save_tiff:
        print(f"Saved TIFF: {output_tiff}")
    print("Done.")
    print("=" * 78)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    try:
        map_data, output_dir = load_and_merge(args)
        make_figure(map_data, output_dir, args)
    except Exception as exc:
        print("\nERROR:", str(exc), file=sys.stderr)
        raise


if __name__ == "__main__":
    main()
