#!/bin/sh
#SBATCH --job-name=bias_coverage
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1

####Running over all 1000 simulations 100 time for each 
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(slurmR)
params_all <- read.csv("~/sima/CopyOfSims_calibrate_paper/params_all.csv")
# Parameters
n_runs        <- 100
ndays         <- 60
n_param_sets  <- max(params_all$sim_id)
cat("Starting bias and coverage analysis for", n_param_sets, "parameter sets...\n")
cat("Each parameter set will be run", n_runs, "times for", ndays, "days\n")
set.seed(123)

# ── Core simulation function ─────────────────────────────────────────────────
run_parameter_set_method <- function(sim_idx, method_name) {
  
  library(epiworldR)
  library(dplyr)
  
  method_params <- params_all %>%
    filter(sim_id == sim_idx, param_type == method_name)
  
  if (nrow(method_params) == 0) {
    cat("Warning: No parameters found for sim_id", sim_idx, "method", method_name, "\n")
    return(tibble())
  }
  
  model <- ModelSIRCONN(
    name              = paste0("sim_", sim_idx, "_", method_name),
    n                 = method_params$true_n,
    prevalence        = method_params$prevalence,
    contact_rate      = method_params$contact_rate,
    transmission_rate = method_params$transmission_rate,
    recovery_rate     = method_params$recovery_rate
  )
  
  # ── Save both total_hist AND transition ──────────────────────────────────
  saver <- epiworldR::make_saver("total_hist", "transition")
  epiworldR::run_multiple(model, ndays = ndays, nsims = n_runs,
                          saver = saver, nthreads = 18)
  ans <- epiworldR::run_multiple_get_results(model)
  
  # ── Infected prevalence stats (total_hist) ────────────────────────────────
  infected_stats <- ans$total_hist %>%
    filter(state == "Infected") %>%
    group_by(date) %>%
    summarize(
      mean_infected = mean(counts),
      std_infected  = sd(counts),
      q025          = quantile(counts, 0.025),
      q975          = quantile(counts, 0.975),
      .groups       = "drop"
    ) %>%
    mutate(sim_id     = sim_idx,
           param_type = method_name)
  
  # ── Daily incidence stats (S → I transitions) ─────────────────────────────
  incidence_stats <- ans$transition %>%
    filter(from == "Susceptible", to == "Infected", date > 0) %>%
    group_by(date) %>%
    summarize(
      mean_incidence = mean(counts),
      std_incidence  = sd(counts),
      q025_incidence = quantile(counts, 0.025),
      q975_incidence = quantile(counts, 0.975),
      .groups        = "drop"
    ) %>%
    mutate(sim_id     = sim_idx,
           param_type = method_name)
  
  gc()
  
  # Return both as a named list
  return(list(
    infected_stats  = infected_stats,
    incidence_stats = incidence_stats
  ))
}

# ── Method wrappers ───────────────────────────────────────────────────────────
run_true_method <- function(sim_idx) run_parameter_set_method(sim_idx, "true")
run_abc_method  <- function(sim_idx) run_parameter_set_method(sim_idx, "abc")
run_lstm_method <- function(sim_idx) run_parameter_set_method(sim_idx, "lstm")

# ── SLURM submission ──────────────────────────────────────────────────────────
slurm_opts <- list(
  partition      = "vegayon-np",
  account        = "vegayon-np",
  time           = "02:00:00",
  `mem-per-cpu`  = "8G",
  `cpus-per-task`= 1
)

cat("Submitting SLURM jobs for method: true\n")
ans <- Slurm_lapply(
  X        = 1:n_param_sets,
  FUN      = run_true_method,
  job_name = "bias_coverage_true",
  njobs    = 100,
  overwrite= TRUE,
  plan     = "submit",
  sbatch_opt = slurm_opts,
  export   = c("run_parameter_set_method", "run_true_method",
               "params_all", "n_runs", "ndays")
)

cat("Submitting SLURM jobs for method: abc\n")
ans1 <- Slurm_lapply(
  X        = 1:n_param_sets,
  FUN      = run_abc_method,
  job_name = "bias_coverage_abc",
  njobs    = 100,
  overwrite= TRUE,
  plan     = "submit",
  sbatch_opt = slurm_opts,
  export   = c("run_parameter_set_method", "run_abc_method",
               "params_all", "n_runs", "ndays")
)

cat("Submitting SLURM jobs for method: lstm\n")
ans2 <- Slurm_lapply(
  X        = 1:n_param_sets,
  FUN      = run_lstm_method,
  job_name = "bias_coverage_lstm",
  njobs    = 100,
  overwrite= TRUE,
  plan     = "submit",
  sbatch_opt = slurm_opts,
  export   = c("run_parameter_set_method", "run_lstm_method",
               "params_all", "n_runs", "ndays")
)

library(slurmR)
library(dplyr)

# Collect all three methods
results_true <- Slurm_collect(read_slurm_job("~/sima/CopyOfSims_calibrate_paper/bias_coverage_true"), any. = TRUE)
results_abc  <- Slurm_collect(read_slurm_job("~/sima/CopyOfSims_calibrate_paper/bias_coverage_abc"),  any. = TRUE)
results_lstm <- Slurm_collect(read_slurm_job("~/sima/CopyOfSims_calibrate_paper/bias_coverage_lstm"), any. = TRUE)

# Extract infected stats
infected_true <- do.call(rbind, lapply(results_true, function(x) x$infected_stats))
infected_abc  <- do.call(rbind, lapply(results_abc,  function(x) x$infected_stats))
infected_lstm <- do.call(rbind, lapply(results_lstm, function(x) x$infected_stats))

# Extract incidence stats
incidence_true <- do.call(rbind, lapply(results_true, function(x) x$incidence_stats))
incidence_abc  <- do.call(rbind, lapply(results_abc,  function(x) x$incidence_stats))
incidence_lstm <- do.call(rbind, lapply(results_lstm, function(x) x$incidence_stats))

# Save to CSV
write.csv(infected_true,  "infected_true.csv",  row.names = FALSE)
write.csv(infected_abc,   "infected_abc.csv",   row.names = FALSE)
write.csv(infected_lstm,  "infected_lstm.csv",  row.names = FALSE)
write.csv(incidence_true, "incidence_true.csv", row.names = FALSE)
write.csv(incidence_abc,  "incidence_abc.csv",  row.names = FALSE)
write.csv(incidence_lstm, "incidence_lstm.csv", row.names = FALSE)

cat("Done! Results saved.\n")
