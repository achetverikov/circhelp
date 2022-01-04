## code to prepare `Pascucci_et_al_2019_data` dataset goes here

Pascucci_et_al_2019_data <- fread('https://zenodo.org/record/2544946/files/Experiment2_rawdata.csv?download=1')
Pascucci_et_al_2019_data[, err:=angle_diff_180(reported, orientation)]
usethis::use_data(Pascucci_et_al_2019_data, overwrite = TRUE)
