## Overview
This script processes X-ray and radio surface brightness data from FITS files and performs Bayesian linear regression using the **LinMix** package. It applies data filtering, logarithmic transformations, and censorship handling before plotting the results.

## Requirements
Ensure you have the following dependencies installed:

```bash
pip install numpy matplotlib astropy linmix
```

## Usage
Run the script with the following command-line arguments:

```bash
python plot_radio_xray.py -c <cluster_id> -n <censorship>
```

### Arguments:
- `-c`, `--cluster`: Integer, specifies the cluster ID (e.g., `697` or `2219`). It is only used for labeling and threshold selection.
- `-n`, `--cens`: Integer, determines censorship mode:
  - `1`: Applies censorship based on threshold.
  - `0`: Excludes data below `2σ` threshold.

## Data Preparation
Before running the script, specify the correct paths for the FITS files:

```python
file_RH = "<filename_RH>"  # Replace with actual path to RH FITS file
file_MH = "<filename_MH>"  # Replace with actual path to MH FITS file
```

Additionally, verify that the column names match those in your FITS files. Use:
```python
print(RH_table.colnames)
print(MH_table.colnames)
```
Update the column names accordingly in:
```python
Ix_RH, Ir_RH = RH_table['xray_sb'], RH_table['radio1_sb']
Ix_RH_err, Ir_RH_err = RH_table['xray_sb_err'], RH_table['radio1_sb_err']
```

## Cluster-Specific Information
The script contains hardcoded threshold values for specific clusters (e.g., `697`, `2219`). If you are working with different clusters, ensure that you modify these values accordingly in the script.

## Workflow
1. **Load Data**: Reads X-ray and radio surface brightness from FITS files.
2. **Filter Data**: Removes invalid and non-finite values.
3. **Logarithmic Transformation**: Converts brightness values to log-scale.
4. **Apply Censorship**: Defines a censorship array when required.
5. **Bayesian Regression**: Uses `linmix.LinMix` for linear fitting.
6. **Plot Results**: Displays fitted regression lines and uncertainty regions.

## Output
- A **scatter plot** with error bars showing surface brightness correlations.
- Best-fit regression lines with shaded uncertainty regions.

## Example Command
```bash
python plot_radio_xray.py -c 697 -n 1
```

## Notes
- Ensure that the FITS files exist and their paths are correctly set.
- Adjust `plt.xlim` and `plt.ylim` manually if necessary for better visualization.

## License
MIT License. Feel free to modify and use the script as needed.

