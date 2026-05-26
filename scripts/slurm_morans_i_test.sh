#!/bin/bash
#SBATCH -p normal          # CPU partition — Moran's I does not use GPU
#SBATCH --cpus-per-task=4  # 4 CPUs sufficient for 999 permutations on 253 counties
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J morans_i
#SBATCH -o /home/yussifyahaya01/texas_diabetes_map/logs/morans_i_%j.out
#SBATCH -e /home/yussifyahaya01/texas_diabetes_map/logs/morans_i_%j.err

set -euo pipefail

echo "============================================================"
echo "Moran's I Spatial Autocorrelation Test"
echo "Job started on: $(hostname)"
echo "Start time    : $(date)"
echo "============================================================"

PROJECT_DIR=/home/yussifyahaya01/texas_diabetes_map
SCRIPT_PATH=${PROJECT_DIR}/scripts/morans_i_test.py
DATA_CSV=${PROJECT_DIR}/data/texas_analysis_dataset_v2.csv
SHAPEFILE_ZIP=${PROJECT_DIR}/data/cb_2023_us_county_5m.zip
OUTPUT_DIR=${PROJECT_DIR}/results/morans_i_test

cd ${PROJECT_DIR} || exit 1
mkdir -p ${OUTPUT_DIR}
mkdir -p ${PROJECT_DIR}/logs

source /home/yussifyahaya01/miniconda3/etc/profile.d/conda.sh || exit 1
conda activate geofig || exit 1

# Match thread count to SLURM CPU allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

echo "Project dir   : ${PROJECT_DIR}"
echo "Python path   : $(which python)"
python --version

echo "Checking input files..."
test -f "${SCRIPT_PATH}" || { echo "ERROR: Missing script: ${SCRIPT_PATH}"; exit 1; }
test -f "${DATA_CSV}" || { echo "ERROR: Missing CSV: ${DATA_CSV}"; exit 1; }
test -f "${SHAPEFILE_ZIP}" || { echo "ERROR: Missing shapefile ZIP: ${SHAPEFILE_ZIP}"; exit 1; }

echo "Testing required Python packages..."
python -c "import geopandas, libpysal, esda, statsmodels, pandas, numpy, matplotlib; print('All required packages are available.')"

echo "Running Moran's I test..."
python ${SCRIPT_PATH} \
    --data_csv ${DATA_CSV} \
    --shapefile_zip ${SHAPEFILE_ZIP} \
    --output_dir ${OUTPUT_DIR} \
    --permutations 999

echo "============================================================"
echo "Output files:"
ls -lh ${OUTPUT_DIR}
echo "End time: $(date)"
echo "Job finished successfully."
echo "============================================================"
