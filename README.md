# The Role of International Trade in The Great Convergence

Replication code for the Master Thesis by Florian Gerzner (University of Zurich, 2026).

## Overview
This repository contains all MATLAB code to replicate the results in the thesis.
The analysis uses a static multi-sector Armington model with input-output linkages,
calibrated on data from the World Input-Output Database (WIOD, 2013 release).

## Repository Structure

| File | Description |
|---|---|
| `create_DATA_files.m` | Cleans and aggregates raw WIOD tables into `DATA_YY.mat` files |
| `parameters.m` | Computes structural parameters and saves `PARAMETERS_YY_AGG.mat` |
| `baseline_model.m` | One-sector, no-intermediates Armington model |
| `multisector_model.m` | Multi-sector, no-intermediates Armington model |
| `linkage_model.m` | Full model with sectors and intermediate input linkages |
| `plot_trade_shares.m` | Gives the code for the trade share plot in the paper |
| `plot_trade_volumes.m` | Gives the code for the trade volumes plot in the paper |

## Data
The raw WIOD data (`wiotYY_row_apr12.csv`) is not included in this repository
and must be obtained from (https://www.rug.nl/ggdc/valuechain/wiod/wiod-2013-release).
The aggregation mappings and elasticity files (`country_aggregation.csv`,
`sector_aggregation.csv`, `SIGMA.csv`, `country_list.csv`, `sector_list.csv`)
are included.
It is not necessary to load the files as all the aggregated trade matrices and parameters are included.

## How to Run

**If the `DATA_YY.mat` and `PARAMETERS_YY_AGG.mat` files are not included:**
1. Obtain the raw WIOD CSV files from [www.wiod.org](http://www.wiod.org) and place them in the same folder as the scripts
2. Run `create_DATA_files.m` to produce `DATA_YY.mat` files
3. Run `clean_aggregate_WIOD.m` to produce `PARAMETERS_YY_AGG.mat` files

**To reproduce the results:**

4. Run `baseline_model.m` — one sector, no intermediates
5. Run `multisector_model.m` — multiple sectors, no intermediates
6. Run `linkage_model.m` — full model with intermediate input linkages

## Reference
Costinot, A., & Rodríguez-Clare, A. (2014). Trade theory with numbers: Quantifying the consequences of globalization. In Handbook of international economics (Vol. 4, pp. 197-261). Elsevier.
Timmer, M. P., Dietzenbacher, E., Los, B., Stehrer, R., & De Vries, G. J. (2015). An illustrated user guide to the world input–output database: the case of global automotive production. Review of international economics, 23(3), 575-605.
