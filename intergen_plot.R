#!/usr/bin/env -S Rscript --vanilla

# attempt at a faster version on intergen_comp.R

# libraries
library(here)
source(here::here("Base_functions.R"))

file_pattern <- "set_intergen3_[[:graph:]]+_s_"
file_path    <- here(dirname(here()), "src", "files_folder", "Intergen")
file_list    <- list.files(file_path, pattern = file_pattern)
file_out     <- here("Data", "intergen3_s.csv")

run_header     <- read_run(file_path, file_list[[1]])
contour_header <- contour_data(run_header$data, run_header$params,
				cost_var = "survival_cost_of_surv_help",
				keep_var = c("d", "baseline_survival", "juvenile_survival_weight"))
writeLines(paste(names(contour_header), collapse = ","), file_out)
# print(names(contour_header))

CON <- file(file_out, "a")

lapply(file_list, function(file_current){
	run_current     <- read_run(file_path, file_current)
	contour_current <- contour_data(run_current$data, run_current$params,
					cost_var = "survival_cost_of_surv_help",
					keep_var = c("d", "baseline_survival", "juvenile_survival_weight"))
	writeLines(apply(contour_current, 1, paste, collapse = ","), CON)
})
close(CON)

