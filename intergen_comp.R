#!/usr/bin/env -S Rscript --vanilla

# Script for use on server/from terminal instead of RStudio
# Written by Margaret Bolton 17/06/24

# libraries
library(here)
source(here::here("Base_functions.R"))

file_number <- 10
file_pattern <- paste0("set_intergen3_[[:graph:]]+_s_", file_number, "$")
file_out <- paste0("intergen3_s_", file_number, ".csv")

# import data
intergen_set_data <- read_set_data(here::here(dirname(here::here()), "src", "files_folder", "Intergen"),
                                   pattern = file_pattern)

intergen_set_params <- read_set_params(here::here(dirname(here::here()), "src", "files_folder", "Intergen"),
                                       pattern = file_pattern)

# prep data
intergen_contour_dat <- contour_data(intergen_set_data, intergen_set_params,
                                     cost_var = "survival_cost_of_surv_help",
                                     keep_var = c("d", "baseline_survival", 
                                                  "juvenile_survival_weight"))

write.csv(intergen_contour_dat, here::here("Data", file_out), row.names = FALSE)
