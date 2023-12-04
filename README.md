# HEI_Australian_wildfires
This repository contains the scripts for data processing procedures and data analysis
1. calculate_background_conc_new_fast.m & calculate_background_conc_new_fast_se_aug2jan.m: Calculate non-fire (background) PM2.5 concentrations at each receptors.
2. calculate_ratio_above_pblh.m: Calculate plume injection fractions based on MISR plume height data.
3. calculate_smoke_exposure_full_fast.m & calculate_smoke_exposure_full_fast_202001.m: Calculate daily smoke exposure based on INJ-CLIM and CTL experiments.
4. calculate_smoke_exposure_rf.m: Calculate smoke exposure based on INJ-RF experiment.
5. create_gfed_files.m & create_gfed_files_1997_2002.m: Calculate fire emission fluxes from GFED v4.1s original files and generate input files for HEMCO.
6. create_input_factors_for_rf_prediction.m: Create input variables for random forest predictions.
7. prepare_io_variables_for_RF.m: Create input/output variables for training the random forest model.
8. random_forest_regression_training.m: Train the random forest regression model.
9. read_merra2_var.m: Read, crop and save meteorological variables from MERRA2 files (Variables: PBLH, T2, U10M, V10M, RH, and PREP).
10. read_save_mcd14ml.m: Read and save MODIS fire data (MCD14ML).
11. read_save_plume_data.m: Read and save MISR plume height data.
12. demo_calculate_smoke_exposure.mlx: Demo for calculating smoke exposure in Nov 2019.
13. demo_train_random_forest_regression.mlx: Demo for training the random forest regression model.
14. demo_rf_prediction_daily_injection_fractions_au.mlx: Demo for predicting daily plume injection fractions using resulting random forest model.
