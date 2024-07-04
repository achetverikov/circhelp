## code to prepare `Bae_Luck_2018_data` dataset goes here

Bae_Luck_2018_data <- fread("https://osf.io/4m2kb/download")
Bae_Luck_2018_data[, Duration := NULL]
Bae_Luck_2018_data[, AngularError := NULL]
Bae_Luck_2018_data[, TarBinID := NULL]
Bae_Luck_2018_data[, BinCenter := NULL]
Bae_Luck_2018_data[, err := angle_diff_360(RespAngle, TargetDirection)]
usethis::use_data(Bae_Luck_2018_data, overwrite = TRUE)
