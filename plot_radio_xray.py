import argparse
import os
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import linmix

# ----------------------------
# 0. Parse command-line arguments
# ----------------------------

parser = argparse.ArgumentParser(description="Process X-ray and radio surface brightness data.")

parser.add_argument(
    "-c", "--cluster",
    type=int,
    required=True,
    help="Which cluster (e.g., 697)?"
)

parser.add_argument(
    "-n", "--cens",
    type=int,
    required=True,
    help="Censorship on = 1; exclude <2RMS = 0"
)

args = parser.parse_args()

# ----------------------------
# 1. Load data from FITS files
# ----------------------------

# TODO: Specify the correct file paths before running the script.
# Replace "<filename_RH>" and "<filename_MH>" with the actual FITS file paths.

file_RH = "<filename_RH>"  # Replace with actual path to RH FITS file
file_MH = "<filename_MH>"  # Replace with actual path to MH FITS file

print(f"Cluster ID: {args.cluster}")
print(f"RH file: {file_RH}")
print(f"MH file: {file_MH}")

# Check if files exist before proceeding
if not os.path.exists(file_RH) or not os.path.exists(file_MH):
    raise FileNotFoundError(f"Error: FITS file(s) not found: {file_RH} or {file_MH}")

# Read data from FITS files
RH_table = Table.read(file_RH, format='fits')
MH_table = Table.read(file_MH, format='fits')

# ----------------------------
# 2. Extract relevant columns
# ----------------------------

# TODO: Verify that these column names match the FITS file structure.
# If the column names differ, update them accordingly.
# Use `print(RH_table.colnames)` or `print(MH_table.colnames)` to check available columns.

print("Available columns in RH file:", RH_table.colnames)
print("Available columns in MH file:", MH_table.colnames)

# Update these column names if necessary
Ix_RH, Ir_RH = RH_table['xray_sb'], RH_table['radio1_sb']
Ix_RH_err, Ir_RH_err = RH_table['xray_sb_err'], RH_table['radio1_sb_err']

Ix_MH, Ir_MH = MH_table['xray_sb'], MH_table['radio1_sb']
Ix_MH_err, Ir_MH_err = MH_table['xray_sb_err'], MH_table['radio1_sb_err']

# ----------------------------
# 3. Filter data (remove invalid values)
# ----------------------------

# Define threshold based on the cluster ID
threshold = {697: 1.96e-7, 2219: 4.42e-7}.get(args.cluster, 0)
filter_r = threshold if args.cens == 0 else 0

valid_MH = (Ix_MH > 0) & (Ir_MH > filter_r) & np.isfinite(Ix_MH) & np.isfinite(Ir_MH) & np.isfinite(Ix_MH_err) & np.isfinite(Ir_MH_err)
valid_RH = (Ix_RH > 0) & (Ir_RH > 0) & np.isfinite(Ix_RH) & np.isfinite(Ir_RH) & np.isfinite(Ix_RH_err) & np.isfinite(Ir_RH_err)

Ix_MH, Ir_MH, Ix_MH_err, Ir_MH_err = Ix_MH[valid_MH], Ir_MH[valid_MH], Ix_MH_err[valid_MH], Ir_MH_err[valid_MH]
Ix_RH, Ir_RH, Ix_RH_err, Ir_RH_err = Ix_RH[valid_RH], Ir_RH[valid_RH], Ix_RH_err[valid_RH], Ir_RH_err[valid_RH]

# ----------------------------
# 4. Logarithmic transformation
# ----------------------------

log_Ix_RH, log_Ir_RH = np.log10(Ix_RH), np.log10(Ir_RH)
log_Ix_MH, log_Ir_MH = np.log10(Ix_MH), np.log10(Ir_MH)

log_Ix_RH_err = Ix_RH_err / (Ix_RH * np.log(10))
log_Ir_RH_err = Ir_RH_err / (Ir_RH * np.log(10))
log_Ix_MH_err = Ix_MH_err / (Ix_MH * np.log(10))
log_Ir_MH_err = Ir_MH_err / (Ir_MH * np.log(10))

# ----------------------------
# 5. Create censorship array
# ----------------------------

censorship_r = np.zeros_like(Ir_MH, dtype=int) if args.cens == 0 else (Ir_MH > threshold).astype(int)

# ----------------------------
# 6. Bayesian linear regression with LinMix
# ----------------------------

linmix_RH = linmix.LinMix(log_Ix_RH, log_Ir_RH, log_Ix_RH_err, log_Ir_RH_err)
linmix_RH.run_mcmc(miniter=2500, silent=True)

linmix_MH = linmix.LinMix(log_Ix_MH, log_Ir_MH, log_Ix_MH_err, log_Ir_MH_err, delta=censorship_r)
linmix_MH.run_mcmc(miniter=2500, silent=True)

# ----------------------------
# 7. Extract best-fit parameters
# ----------------------------

def extract_fit_params(linmix_result):
    alpha = np.mean(linmix_result.chain['alpha'])
    beta = np.mean(linmix_result.chain['beta'])
    alpha_std = np.std(linmix_result.chain['alpha'])
    beta_std = np.std(linmix_result.chain['beta'])
    return alpha, beta, alpha_std, beta_std

alpha_RH, beta_RH, alpha_RH_std, beta_RH_std = extract_fit_params(linmix_RH)
alpha_MH, beta_MH, alpha_MH_std, beta_MH_std = extract_fit_params(linmix_MH)

# ----------------------------
# 8. Plot results
# ----------------------------

plt.figure(figsize=(10, 8))

plt.errorbar(Ix_RH, Ir_RH, xerr=Ix_RH_err, yerr=Ir_RH_err, fmt='o', markersize=5, capsize=2, color='blue', alpha=0.5, label='RH')
plt.errorbar(Ix_MH, Ir_MH, xerr=Ix_MH_err, yerr=Ir_MH_err, fmt='o', markersize=5, capsize=2, color='red', alpha=0.5, label='MH')

plt.xscale('log')
plt.yscale('log')

# Updated labels with LaTeX formatting
plt.xlabel(r'$\log(\text{Surface brightness in X-band}) \quad [\text{count}/(\text{arcsec}^{2} \, \text{s})]$')
plt.ylabel(r'$\log(\text{Surface brightness in radio-band}) \quad [\text{Jy}/\text{arcsec}^{2}]$')

plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()

