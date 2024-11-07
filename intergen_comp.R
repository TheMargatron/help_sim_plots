#!/usr/bin/env -S Rscript --vanilla

# Script for use on server/from terminal instead of RStudio
# Written by Margaret Bolton 17/06/24

# libraries
library(here)
source(here::here("Base_functions.R"))

# import data
intergen_set_data <- read_set_data(here::here(dirname(here::here()), "src", "files_folder", "Intergen"),
                                   pattern = "intergen3_[[:graph:]]+_f_10$")

intergen_set_params <- read_set_params(here::here(dirname(here::here()), "src", "files_folder", "Intergen"),
                                       pattern = "intergen3_[[:graph:]]+_f_10$")

# prep data
intergen_contour_dat <- contour_data(intergen_set_data, intergen_set_params,
                                     cost_var = "survival_cost_of_surv_help",
                                     keep_var = c("d", "baseline_survival", 
                                                  "juvenile_survival_weight"))

write.csv(intergen_contour_dat, here::here("Data", "intergen3_f_10.csv"), row.names = FALSE)
