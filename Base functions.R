# Functions for reading in and making basic plots of simulation data
# Written by Margaret Bolton (mb804@exeter.ac.uk) 23/04/24

# library ####
library(tidyverse)
library(here)

# ______________________________________________________________________________
# reading data ####
# Find the first line with no written character
find_split <- function(file_name){
  file_temp <- readLines(file_name)
  split_index <- match("", file_temp)
  return(split_index)
}

# reads single run returning data and parameters in list, for tests rather than full sets
read_run <- function(folder_name, file_name){
  # find line in file where it changes from data to params, otherwise throw an error
  split_index <- find_split(paste(folder_name, file_name, sep = "/"))
  if(is.na(split_index)){
    stop("No split index found")
  }
  
  # read parameter data from file by skipping to split_index, then rename vars and add file info
  params <- read.table(paste(folder_name, file_name, sep = "/"), header = FALSE, sep = ";", skip = split_index)
  names(params) <- c("param", "val")
  params$folder_name <- folder_name
  params$file_name <- file_name
  
  # read data from file by stopping at split_index-2, remove empty variable, add file info
  dat <- read.table(paste(folder_name, file_name, sep = "/"), header = TRUE, sep = ";", nrows = (split_index - 2)) # subtract one for the header row, another for the first blank line
  dat$X <- NULL
  dat$foldername <- folder_name
  dat$filename <- file_name
  
  return(list("params" = params, "data" = dat))
}

# reads data from full set and binds into dataframe
read_set_data <- function(folder_name, pattern = NULL){
  # list files in folder, can include name pattern to isolate sets
  file_list <- list.files(folder_name, pattern = pattern, full.names = FALSE)
  
  # find line in each file where it changes from data to params, otherwise throw an error
  split_index <- unlist(lapply(paste(folder_name, file_list, sep = "/"), find_split))
  if(anyNA(split_index)) {
    print(paste0("No split index found for files ", which(is.na(split_index))))
    na.fail(split_index)
  }
  
  # read data from all files, remove stray empty var, add file info
  dat_list <- lapply(file_list, function(fn) {
    dat <- read.table(paste(folder_name, fn, sep = "/"), 
                      header = TRUE, sep = ";", 
                      nrows = (find_split(paste(folder_name, fn, sep = "/")) - 2))
    dat$X <- NULL
    dat$foldername <- folder_name
    dat$filename <- fn
    return(dat)
  })
  
  dat <- do.call(plyr::rbind.fill, dat_list)
  dat$filename <- as.factor(dat$filename)
  
  return(dat)
}

# Reads parameter values from full set and binds into dataframe
read_set_params <- function(folder_name, pattern = NULL){
  # list files in folder, can include name pattern to isolate sets
  file_list <- list.files(folder_name, pattern = pattern, full.names = FALSE)
  
  # find line in each file where it changes from data to params, otherwise throw an error
  split_index <- unlist(lapply(paste(folder_name, file_list, sep = "/"), find_split))
  if(anyNA(split_index)) {
    print(paste0("No split index found for files ", which(is.na(split_index))))
    na.fail(split_index)
  }
  
  # read parameter data from all files by skipping to split_index, then rename vars and add file info
  param_list <- lapply(file_list, function(fn) {
    dat <- read.table(paste(folder_name, fn, sep = "/"), 
                      header = FALSE, sep = ";", 
                      skip = find_split(paste(folder_name, fn, sep = "/")))
    names(dat) <- c("param", "val")
    dat$foldername <- folder_name
    dat$filename <- fn
    return(dat)
  })
  
  params <- do.call(rbind, param_list)
  return(params)
  
}

# Adds missing rows to parameter data
bind_params <- function(missing_vars, params_dat, filenams, foldernam = "Data/sets"){
  # Provide matrix with columns named "param" and "val" 
  missing_vars <- missing_vars %>% 
    as.data.frame() %>% 
    mutate(val = as.numeric(val)) %>% 
    bind_cols(as.data.frame(matrix(rep(unique(filenams), each = nrow(.)), 
                                   byrow  = FALSE, 
                                   ncol   = length(unique(filenams))))) %>% 
    pivot_longer(cols = paste0("V", 1:length(unique(filenams))), names_to = NULL, values_to = "filename") %>% 
    mutate(foldername = foldernam)
  
  # Join and return
  params_dat <- bind_rows(params_dat, missing_vars)
  
  return(params_dat)
  
}

# gets final value for each run within a set 
extract_final <- function(dat, extinct = FALSE){
  # alternative method would be to provide time step but I doubt I would use that
  if(extinct){
    ## adjust for when population goes extinct
    min_time <- dat %>% 
      group_by(filename) %>% 
      summarise(max_time = max(time_step)) %>% 
      {min(.$max_time)}
    
    final <- dat %>% 
      group_by(filename) %>% 
      dplyr::filter(time_step == min_time)
    
  } else {
    ## do things normally
    final <- dat %>%
      group_by(filename) %>%
      dplyr::filter(time_step == max(time_step))
    
  }
  
  return(as.data.frame(final))
}

# ______________________________________________________________________________
# plotting ####

# plot timeseries 
timeseries_plot <- function(ts_dat, 
                            y, 
                            title_text = "longnam", 
                            label_ylab = TRUE,
                            aes_col = "filename",
                            legend_title = "Run file"){
  plot_basic <- ggplot(ts_dat, aes(x = !! sym("time_step"), y = !! sym(y), colour = !! sym(aes_col))) +
    geom_point(size = 1) + 
    labs(title = plot_names[plot_names$varnam == y, title_text],
         colour = legend_title) +
    xlab(plot_names[plot_names$varnam == "time_step", "shortnam"]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  if(label_ylab){
    plot_basic <- plot_basic + ylab(plot_names[plot_names$varnam == y, "shortnam"]) 
  } else{
    plot_basic <- plot_basic + ylab(NULL)
  }
  return(plot_basic)
}

## contour plots ####
# TODO switch to ggplot throughout
# https://stackoverflow.com/questions/73949067/control-label-of-contour-lines-in-contour

contour_data <- function(dat, params_dat, 
                         cost_var = "fecundity_cost_of_fec_help",
                         x_var = "baseline_survival",
                         y_var = "(?!x)x"){
  
  contour_dat <- dat %>% 
    select(-contains("disp")) %>% 
    filter(time_step == 50000) %>%
    select(time_step, ends_with("1"), foldername, filename) %>% 
    rename_with(~ str_remove_all(.x, "[:digit:]"))
  
  contour_dat <- dat %>% 
    select(-contains("disp")) %>% 
    filter(time_step == 50000) %>% 
    select(time_step, ends_with("2"), foldername, filename) %>% 
    rename_with(~ str_remove_all(.x, "[:digit:]"))  %>% 
    bind_rows(contour_dat)
  
  contour_dat <- params_dat %>% 
    mutate(param = str_remove_all(param, "[:digit:]")) %>% 
    distinct() %>% 
    filter(str_detect(param, cost_var) | 
             str_detect(param, x_var) |
             str_detect(param, y_var)) %>% 
    pivot_wider(names_from = param,
                values_from = val) %>% 
    right_join(contour_dat, join_by(filename, foldername)) %>% 
    mutate(b_over_c = 1/get(cost_var)) 
  
  return(contour_dat)
}

contour_matrix <- function(contour_dat,
                           x_var = "baseline_survival",
                           help_var = "mean_fec_h",
                           summary_fun = mean){
  
  contour_mat <- contour_dat  %>% 
    select(!!sym(x_var), !! sym(help_var),
           b_over_c) %>% 
    arrange(b_over_c) %>% 
    pivot_wider(names_from = b_over_c,
                values_from = !! sym(help_var),
                values_fn = summary_fun) %>% 
    select(!! sym(x_var), matches("[[:digit:]]")) %>% 
    arrange(!!x_var) 
  
  return(contour_mat)

}

contour_points <- function(contour_dat,
                           x_var = "baseline_survival",
                           y_var = "b_over_c"){
  
  contour_pnt <- contour_dat %>% 
    select(!! sym(x_var), !! sym(y_var)) %>% 
    distinct()
  
  return(contour_pnt)
  
}

basic_contour <- function(contour_dat, 
                          x_var = "baseline_survival",
                          help_var = "mean_fec_h",
                          main = help_var,
                          SD_lines = FALSE){
  
  contour_mat <- contour_matrix(contour_dat,
                                x_var = x_var,
                                help_var = help_var)
  
  contour_pnt <- contour_points(contour_dat, x_var = x_var, y_var = "b_over_c")
  
  if(!SD_lines){
    contour_out <- filled.contour(x = contour_mat[[x_var]], 
                                  y = as.numeric(names(contour_mat)[2:length(names(contour_mat))]),
                                  z = as.matrix(contour_mat[,2:ncol(contour_mat)]), plot.axes = {
                                    axis(1)
                                    axis(2)
                                    contour(x = contour_mat[[x_var]],
                                            y = as.numeric(names(contour_mat)[2:length(names(contour_mat))]),
                                            z = as.matrix(contour_mat[,2:ncol(contour_mat)]),
                                            add = TRUE,
                                            levels = 0
                                    )
                                    points(contour_pnt)
                                  },
                                  xlab = x_var,
                                  ylab = "b_over_c",
                                  main = main
    )
    
  } else {
    plus_sd <- contour_matrix(contour_dat,
                              x_var = x_var,
                              help_var = help_var,
                              summary_fun = function(x){mean(x) + sd(x)})
    
    minus_sd <- contour_matrix(contour_dat,
                               x_var = x_var,
                               help_var = help_var,
                               summary_fun = function(x){mean(x) - sd(x)})
    
    contour_out <- filled.contour(x = contour_mat[[x_var]], 
                                  y = as.numeric(names(contour_mat)[2:length(names(contour_mat))]),
                                  z = as.matrix(contour_mat[,2:ncol(contour_mat)]), plot.axes = {
                                    axis(1)
                                    axis(2)
                                    contour(x = contour_mat[[x_var]],
                                            y = as.numeric(names(contour_mat)[2:length(names(contour_mat))]),
                                            z = as.matrix(contour_mat[,2:ncol(contour_mat)]),
                                            add = TRUE,
                                            levels = 0
                                    )
                                    contour(x = plus_sd[[x_var]],
                                            y = as.numeric(names(plus_sd)[2:length(names(plus_sd))]),
                                            z = as.matrix(plus_sd[,2:ncol(plus_sd)]),
                                            add = TRUE,
                                            levels = 0
                                    )
                                    contour(x = minus_sd[[x_var]],
                                            y = as.numeric(names(minus_sd)[2:length(names(minus_sd))]),
                                            z = as.matrix(minus_sd[,2:ncol(minus_sd)]),
                                            add = TRUE,
                                            levels = 0
                                    )
                                    points(contour_pnt)
                                  },
                                  xlab = x_var,
                                  ylab = "b_over_c",
                                  main = main
    )
    
  }
  
}

meta_contour <- function(contour_dat,
                         x_var = "baseline_survival",
                         y_var = "juvenile_survival_weight", 
                         help_var = "mean_fec_h",
                         threshold_lines = c(0),
                         SD_lines = FALSE){
  
  contour_list <- contour_dat %>% 
    arrange(!!x_var) %>% 
    mutate(split_var = as.factor(!! sym(x_var))) %>% 
    split(f = .$split_var)
  
  contour_list <- lapply(contour_list, function(x){
    x_out <- extract_contour(x, y_var = y_var, help_var = help_var)
    x_out$split_var <- unique(x$split_var)
    
    return(x_out)
  })
  
  contour_meta <- do.call(rbind, contour_list) %>% 
    dplyr::filter(x %in% contour_dat[[y_var]]) %>% 
    rename(!! sym(y_var) := x,
           b_over_c = y) %>% 
    select(-level) %>% 
    pivot_wider(names_from = !! sym(y_var),
                values_from = b_over_c)
  
  contour_pnt <- contour_points(contour_dat, x_var = x_var, y_var = y_var)
  
  if(class(threshold_lines) != "numeric") {
    threshold_lines <- mean(as.matrix(contour_meta[,2:ncol(contour_meta)]))
    warning("using default threshold line")
  }
  
  contour_out <- filled.contour(x = as.numeric(as.character(contour_meta$split_var)), 
                                y = as.numeric(names(contour_meta)[2:length(names(contour_meta))]),
                                z = as.matrix(contour_meta[,2:ncol(contour_meta)]), plot.axes = {
                                  axis(1)
                                  axis(2)
                                  contour(x = as.numeric(as.character(contour_meta$split_var)),
                                          y = as.numeric(names(contour_meta)[2:length(names(contour_meta))]),
                                          z = as.matrix(contour_meta[,2:ncol(contour_meta)]),
                                          add = TRUE,
                                          levels = threshold_lines
                                          )
                                  points(contour_pnt)
                                },
                                xlab = x_var,
                                ylab = y_var,
                                main = "threshold b/c"
  )
  
  # TODO add sd lines
  
  return(contour_out)
  
}


extract_contour <- function(contour_dat,
                            y_var = "juvenile_survival_weight",
                            help_var = "mean_fec_h"){
  
  contour_mat <- contour_matrix(contour_dat = contour_dat, 
                                help_var = help_var,
                                x_var = y_var)
  
  lines_dat <- contourLines(x = contour_mat[, y_var][[1]], 
                            y = as.numeric(names(contour_mat)[2:length(names(contour_mat))]),
                            z = as.matrix(contour_mat[,2:ncol(contour_mat)]),
                            levels = 0) %>% 
    as.data.frame()
  
  return(lines_dat)
    
}


speedy_gam_plot <- function(dat, y = "mean_fec_h", var1 = "fec_b_over_fec_c", var2 = "baseline_survival") {
  b <- mgcv::gam(get(y) ~ s(get(var2), get(var1)), data = dat)
  b <- mgcViz::getViz(b)
  plot(b)
}

# ______________________________________________________________________________
# essential data ####
plot_names <- c(# output variables
  "mean_disp", 
  "mean_fec_h", 
  "mean_surv_h", 
  "mean_given_fec_h", 
  "mean_given_surv_h",
  
  "var_disp", 
  "var_fec_h", 
  "var_surv_h", 
  "var_given_fec_h", 
  "var_given_surv_h",
  
  "mean_surv_prob", 
  "mean_surv_help_per_ind", 
  "patch_occupancy", 
  "nsurvivors", 
  "mean_offspring",
  
  # params
  "d", 
  "npp", 
  "baseline_survival", 
  "baseline_fecundity",
  
  "strength_survival",
  "fecundity_help",
  "survival_help",
  
  "survival_cost_of_surv_help",
  "survival_cost_of_fec_help", 
  "fecundity_cost_of_surv_help",
  "fecundity_cost_of_fec_help",
  
  "mu_fec_h",
  "mu_surv_h") %>% 
  paste0(rep(c("",1,2), each = length(.))) %>% 
  {data.frame(varnam = c(# output variables
    "time_step",
    
    # params
    "mu_fec_h", 
    "mu_surv_h",
    "mu_d",
    "sdmu",
    
    "between_species",
    "death_birth", 
    "partner_mechanism",
    "fidelity_prob",
    "negotiate_once",
    "sd_pcerr",
    .),
    longnam = NA,
    shortnam = NA)} %>% 
  
  mutate(longnam = case_when(# output variables
    varnam ==          "time_step"                     ~ "Time step",
    str_detect(varnam, "mean_disp")                    ~ "Mean dispersal rate",
    str_detect(varnam, "mean_fec_h")                   ~ "Mean fecundity help given",
    str_detect(varnam, "mean_surv_h")                  ~ "Mean survival help given",
    str_detect(varnam, "mean_given_fec_h")             ~ "Mean given fecundity help",
    str_detect(varnam, "mean_given_surv_h")            ~ "Mean given survival help",
    
    str_detect(varnam, "var_disp")                     ~ "Variance in dispersal rate",
    str_detect(varnam, "var_fec_h")                    ~ "Variance in fecundity help given",
    str_detect(varnam, "var_surv_h")                   ~ "Variance in survival help given",
    str_detect(varnam, "var_given_fec_h")              ~ "Variance in given fecundity help",
    str_detect(varnam, "var_given_surv_h")             ~ "Variance in given survival help",
    
    str_detect(varnam, "mean_surv_prob")               ~ "Mean survival probability",
    str_detect(varnam, "mean_surv_help_per_ind")       ~ "Mean survival help received per individual",
    str_detect(varnam, "nsurvivors")                   ~ "Number of adult survivors",
    str_detect(varnam, "mean_offspring")               ~ "Mean number of offspring",
    str_detect(varnam, "patch_occupancy")              ~ "Patch occupancy",
    
    # params
    str_detect(varnam, "mu_fec_h")                     ~ "Fecundity mutation rate",
    str_detect(varnam, "mu_surv_h")                    ~ "Survival mutation rate",
    varnam ==          "mu_d"                          ~ "Dispersal mutation rate",
    varnam ==          "sdmu"                          ~ "Standard deviation of mutation",
    
    varnam ==          "between_species"               ~ "Between species interaction",
    varnam ==          "death_birth"                   ~ "Death birth updating",
    varnam ==          "partner_mechanism"             ~ "Partner mechanism",
    varnam ==          "fidelity_prob"                 ~ "Probability of fidelity",
    varnam ==          "negotiate_once"                ~ "Negotiate once",
    varnam ==          "sd_pcerr"                      ~ "Standard deviation of error in partner choice",
    
    str_detect(varnam, "fecundity_cost_of_fec_help")   ~ "Fecundity cost of fecundity help",
    str_detect(varnam, "fecundity_cost_of_surv_help")  ~ "Fecundity cost of survival help",
    str_detect(varnam, "survival_cost_of_fec_help")    ~ "Survival cost of fecundity help",
    str_detect(varnam, "survival_cost_of_surv_help")   ~ "Survival cost of survival help",
    
    str_detect(varnam, "survival_help")                ~ "Initial survival help",
    str_detect(varnam, "fecundity_help")               ~ "Initial fecundity help",
    str_detect(varnam, "strength_survival")            ~ "Strength of survival help",
    
    str_detect(varnam, "baseline_fecundity")           ~ "Baseline fecundity rate",
    str_detect(varnam, "baseline_survival")            ~ "Baseline survival rate",
    str_detect(varnam, "npp")                          ~ "Number of individuals per patch",
    str_detect(varnam, "d")                            ~ "(Initial) dispersal rate"
  )) %>%
  mutate(shortnam = longnam) %>%
  mutate(longnam = case_when(str_detect(varnam, "1") ~ paste0(longnam, " (sp. 1)"),
                             str_detect(varnam, "2") ~ paste0(longnam, " (sp. 2)"),
                             TRUE ~ longnam))
