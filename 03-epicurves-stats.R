# ───────────────────────────────────
# Complete code to compare epidemic curves for ONE parameter set
# ───────────────────────────────────

# 1) libraries
library(tidyverse)
library(glue)
source("~/sima/Calibrate_ABM/01a-bilstm.R")
source("~/sima/Calibrate_ABM/00-params.R")
ndays <- 60
N_SIMS <- nrow(theta_use)  
library(epiworldR)
# ──────────────────────────────────────────────────────────────────
# 2) Load ABC predictions and run LSTM predictions
# ──────────────────────────────────────────────────────────────────
params_all=read.csv("~/sima/Calibrate_ABM/params_all.csv")
# ──────────────────────────────────────────────────────────────────
# 4) Pick one parameter set for comparison
# ──────────────────────────────────────────────────────────────────
# Identify the sim_ids you want to keep (those with abc values not both over 30)
set.seed(121)
chosen_sim_id <- sample.int(N_SIMS, size = 1)
single_param_set <- params_all %>%
  filter(sim_id == chosen_sim_id)
cat(sprintf("Comparing epidemic curves for sim_id: %d\n", chosen_sim_id))
print(single_param_set)

# ── Setup ────────────────────────────────────────────────────────────────────
n_reps <- 100
single_param_set$run_id <- seq_len(nrow(single_param_set))
single_param_set[2, 7] <- single_param_set[1, 7]

# Containers
all_avg_series  <- tibble()
all_daily_cases <- tibble()

# ── Main Loop ────────────────────────────────────────────────────────────────
for (i in seq_len(nrow(single_param_set))) {
  row <- single_param_set[i, ]
  
  cat(sprintf("Running %s parameters (row %d/%d)...\n",
              row$param_type, i, nrow(single_param_set)))
  
  set.seed(123 + i)
  
  calibrated_model <- ModelSIRCONN(
    name              = paste0("single_", row$param_type, "_", row$sim_id),
    n                 = row$true_n,
    prevalence        = row$prevalence,
    contact_rate      = row$contact_rate,
    transmission_rate = row$transmission_rate,
    recovery_rate     = row$recovery_rate
  )
  
  # Save both total_hist and transition in one go
  saver <- make_saver("total_hist", "transition")
  run_multiple(calibrated_model, ndays = ndays, nsims = n_reps,
               saver = saver, nthreads = 18)
  ans <- run_multiple_get_results(calibrated_model, nthreads = 18)
  
  # ── Prevalence (S/I/R counts) ──────────────────────────────────────────────
  avg_states <- ans$total_hist %>%
    group_by(date, state) %>%
    summarize(
      mean_count = mean(counts),
      ci_lower   = quantile(counts, 0.025),
      ci_upper   = quantile(counts, 0.975),
      .groups    = "drop"
    ) %>%
    mutate(param_type = row$param_type,
           sim_id     = row$sim_id)
  
  # ── Daily Incidence (S → I transitions) ───────────────────────────────────
  incidence <- ans$transition %>%
    filter(from == "Susceptible", to == "Infected") %>%
    filter(date > 0) %>%
    group_by(date) %>%
    summarize(
      mean_cases = mean(counts),
      ci_lower   = quantile(counts, 0.025),
      ci_upper   = quantile(counts, 0.975),
      median     = quantile(counts, 0.5),
      .groups    = "drop"
    ) %>%
    mutate(param_type = row$param_type,
           sim_id     = row$sim_id,
           transition = "New Infections (S→I)")
  
  # ── Daily Recoveries (I → R transitions) ──────────────────────────────────
  recoveries <- ans$transition %>%
    filter(from == "Infected", to == "Recovered") %>%
    filter(date > 0) %>%
    group_by(date) %>%
    summarize(
      mean_cases = mean(counts),
      ci_lower   = quantile(counts, 0.025),
      ci_upper   = quantile(counts, 0.975),
      median     = quantile(counts, 0.5),
      .groups    = "drop"
    ) %>%
    mutate(param_type = row$param_type,
           sim_id     = row$sim_id,
           transition = "Recoveries (I→R)")
  
  # ── Susceptible Remaining (S → S transitions) ─────────────────────────────
  susceptible_remaining <- ans$transition %>%
    filter(from == "Susceptible", to == "Susceptible") %>%
    filter(date > 0) %>%
    group_by(date) %>%
    summarize(
      mean_cases = mean(counts),
      ci_lower   = quantile(counts, 0.025),
      ci_upper   = quantile(counts, 0.975),
      median     = quantile(counts, 0.5),
      .groups    = "drop"
    ) %>%
    mutate(param_type = row$param_type,
           sim_id     = row$sim_id,
           transition = "Susceptible Remaining (S→S)")
  
  # ── Combine daily cases ────────────────────────────────────────────────────
  daily_cases <- bind_rows(incidence, recoveries, susceptible_remaining)
  
  # ── Collect ───────────────────────────────────────────────────────────────
  all_avg_series  <- bind_rows(all_avg_series,  avg_states)
  all_daily_cases <- bind_rows(all_daily_cases, daily_cases)
}

# ── Plot Helpers ─────────────────────────────────────────────────────────────
library(ggplot2)
library(patchwork)

create_daily_cases_plot <- function(data, title = "SIR Model Comparison (Days 0–40)") {
  
  library(dplyr)
  library(ggplot2)
  
  data <- data %>%
    filter(date <= 40) %>%
    mutate(
      transition = factor(
        transition,
        levels = c(
          "Susceptible Remaining (S→S)",
          "New Infections (S→I)",
          "Recoveries (I→R)"
        )
      ),
      param_type = factor(param_type, levels = c("abc", "lstm", "true"))
    )
  
  ggplot(data, aes(x = date, color = param_type)) +
    
    # Confidence ribbons (not for true)
    geom_ribbon(
      data = data %>% filter(param_type != "true"),
      aes(ymin = ci_lower, ymax = ci_upper, fill = param_type),
      alpha = 0.2,
      color = NA
    ) +
    
    # Mean lines
    geom_line(aes(y = mean_cases), linewidth = 1.2) +
    
    facet_wrap(~ transition, scales = "free_y", ncol = 1) +
    
    scale_color_manual(values = c(
      "abc"  = "#E76F51",
      "lstm" = "#2A9D8F",
      "true" = "#3A86FF"
    )) +
    
    scale_fill_manual(values = c(
      "abc"  = "#E76F51",
      "lstm" = "#2A9D8F"
    )) +
    
    labs(
      title = title,
      subtitle = "Daily Counts Over Time by Method",
      x = "Day",
      y = "Average Daily Count (95% CI)",
      color = "Method",
      fill  = "Method"
    ) +
    
    theme_minimal(base_size = 13) +
    
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      strip.text = element_text(size = 13, face = "bold"),
      panel.spacing = unit(1.2, "lines")
    )
}
# ── Generate & Save Plots ─────────────────────────────────────────────────────

p_daily       <- create_daily_cases_plot(all_daily_cases)

ggsave("~/sima/Calibrate_ABM/figures/p_daily_cases.png", p_daily,       width = 12, height = 14, dpi = 300)

# Update the plotting function to increase y-axis breaks
library(ggplot2)
library(patchwork)
library(dplyr)


# -- Plotting Function with CI and axis breaks --
create_sir_plot_detailed <- function(data, state_name, subtitle) {
  plot_data <- data %>% filter(state == state_name)
  
  y_max <- max(plot_data$ci_upper, na.rm = TRUE)
  y_breaks <- pretty(c(0, y_max), n = 6)
  
  ggplot(plot_data, aes(x = date, color = param_type)) +
    geom_ribbon(
      data = plot_data %>% filter(param_type != "true"),
      aes(ymin = ci_lower, ymax = ci_upper, fill = param_type),
      alpha = 0.2, color = NA
    ) +
    geom_line(aes(y = mean_count), size = 1) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      title = subtitle,
      x = "Day",
      y = paste(state_name, "(mean ± 95% CI)"),
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12)
    )
}

# -- Create 3 SIR plots --
p1 <- create_sir_plot_detailed(all_avg_series, "Infected", "Infected Over Time")
p2 <- create_sir_plot_detailed(all_avg_series, "Susceptible", "Susceptible Over Time")
p3 <- create_sir_plot_detailed(all_avg_series, "Recovered", "Recovered Over Time")

# -- Combine Horizontally with Shared Legend and Caption --
final_plot <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +  # <-- Collect legends into one
  plot_annotation(
    title = paste("Epidemic Curves for Randomly Chosen sim_id:", chosen_sim_id),
    caption = "Legend:\nBlue Line = true (no CI)\nRed Line + Light Red = abc (95% CI)\nGreen Line + Light Cyan = lstm (95% CI)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 11, hjust = 0, margin = margin(t = 10))
    )
  )

# -- Save and Show --
ggsave("~/sima/Calibrate_ABM/figures/epidemic_curves_annotated.png", final_plot, width = 12, height = 7, dpi = 300)
print(final_plot)

create_sir_plot <- function(data, state_name, title) {
  # Filter data to only include days 0-40 and the specified state
  plot_data <- data %>% 
    filter(state == state_name, date <= 40)
  
  p <- ggplot(plot_data, aes(x = date, color = param_type)) +
    # Confidence interval ribbon only for non-"true" methods
    geom_ribbon(
      data = plot_data %>% filter(param_type != "true"),
      aes(ymin = ci_lower, ymax = ci_upper, fill = param_type), 
      alpha = 0.2, color = NA
    ) +
    # Mean lines for all methods
    geom_line(aes(y = mean_count), size = 1) +
    # Set x-axis limits to ensure consistent range across plots
    scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
    labs(
      title = title,
      x = "Day", 
      y = paste0("Average ", state_name, " (95% CI)"), 
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, hjust = 0.5)
    )
  
  return(p)
}

# Create the plots with the modified function
p_infected <- create_sir_plot(all_avg_series, "Infected", "Infected Over Time by Method")
p_susceptible <- create_sir_plot(all_avg_series, "Susceptible", "Susceptible Over Time by Method") 
p_recovered <- create_sir_plot(all_avg_series, "Recovered", "Recovered Over Time by Method")
final_plot <- (p_infected | p_susceptible| p_recovered) +
  plot_layout(guides = "collect") +  # <-- Collect legends into one
  plot_annotation(
    title = paste("Epidemic Curves for Randomly Chosen sim_id:", chosen_sim_id),
    caption = "Legend:\nBlue Line = true (no CI)\nRed Line + Light Red = abc (95% CI)\nGreen Line + Light Cyan = lstm (95% CI)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 11, hjust = 0, margin = margin(t = 10))
    )
  )

####Running over all 1000 simulations 100 time for each 
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(slurmR)

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
#here please run 03-saving-slurm.R file 
source("~/sima/Calibrate_ABM/03-saving-slurm.R")

# ── Collect & save results ────────────────────────────────────────────────────
collect_and_save <- function(slurm_obj, method_label) {
  cat("Collecting results for method:", method_label, "\n")
  raw <- Slurm_collect(slurm_obj)
  
  # Drop errored jobs
  clean <- raw[!sapply(raw, function(x) inherits(x, "error") || inherits(x$res, "error"))]
  
  infected_df  <- bind_rows(lapply(clean, function(x) x$infected_stats))
  incidence_df <- bind_rows(lapply(clean, function(x) x$incidence_stats))
  
  saveRDS(infected_df,  paste0("results_infected_",  method_label, ".rds"))
  saveRDS(incidence_df, paste0("results_incidence_", method_label, ".rds"))
  
  cat("Saved results for method:", method_label, "\n")
  list(infected = infected_df, incidence = incidence_df)
}

res_true <- collect_and_save(ans,  "true")
res_abc  <- collect_and_save(ans1, "abc")
res_lstm <- collect_and_save(ans2, "lstm")

# ── Combine ───────────────────────────────────────────────────────────────────
all_infected  <- bind_rows(res_true$infected,  res_abc$infected,  res_lstm$infected)
all_incidence <- bind_rows(res_true$incidence, res_abc$incidence, res_lstm$incidence)

saveRDS(all_infected,  "~/sima/Calibrate_ABM/data-result/all_simulation_results_infected.rds")
saveRDS(all_incidence, "~/sima/Calibrate_ABM/data-result/all_simulation_results_incidence.rds")

cat("Total infected rows:",  nrow(all_infected),  "\n")
cat("Total incidence rows:", nrow(all_incidence), "\n")

# ── (Optional) reload from disk ───────────────────────────────────────────────
# all_infected  <- readRDS("~/sima/Calibrate_ABM/data-result/all_simulation_results_infected.rds")
# all_incidence <- readRDS("~/sima/Calibrate_ABM/data-result/all_simulation_results_incidence.rds")

# ============================================================
# BIAS & COVERAGE HELPER — works for either dataset
# ============================================================
compute_bias_coverage <- function(all_results, value_col = "mean_infected",
                                  q025_col = "q025", q975_col = "q975") {
  
  wide <- all_results %>%
    select(sim_id, date, param_type,
           value = all_of(value_col),
           q025  = all_of(q025_col),
           q975  = all_of(q975_col)) %>%
    pivot_wider(
      names_from  = param_type,
      values_from = c(value, q025, q975),
      names_sep   = "_"
    )
  
  bias_cov <- wide %>%
    mutate(
      bias_abc       = value_abc  - value_true,
      bias_lstm      = value_lstm - value_true,
      abs_bias_abc   = abs(bias_abc),
      abs_bias_lstm  = abs(bias_lstm),
      pct_bias_abc   = (bias_abc  / value_abc)  * 100,
      pct_bias_lstm  = (bias_lstm / value_lstm) * 100,
      coverage_abc   = (value_true >= q025_abc)  & (value_true <= q975_abc),
      coverage_lstm  = (value_true >= q025_lstm) & (value_true <= q975_lstm)
    )
  
  # Coverage by day (quantile-based)
  coverage_by_day <- all_results %>%
    select(sim_id, date, param_type, value = all_of(value_col)) %>%
    pivot_wider(names_from = param_type, values_from = value) %>%
    group_by(date) %>%
    summarize(
      abc_q025      = quantile(abc,  0.025, na.rm = TRUE),
      abc_q975      = quantile(abc,  0.975, na.rm = TRUE),
      lstm_q025     = quantile(lstm, 0.025, na.rm = TRUE),
      lstm_q975     = quantile(lstm, 0.975, na.rm = TRUE),
      coverage_abc  = mean((true >= abc_q025)  & (true <= abc_q975),  na.rm = TRUE) * 100,
      coverage_lstm = mean((true >= lstm_q025) & (true <= lstm_q975), na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  # Daily stats
  daily <- bias_cov %>%
    group_by(date) %>%
    summarize(
      mean_bias_abc      = mean(bias_abc,      na.rm = TRUE),
      mean_bias_lstm     = mean(bias_lstm,     na.rm = TRUE),
      mean_abs_bias_abc  = mean(abs_bias_abc,  na.rm = TRUE),
      mean_abs_bias_lstm = mean(abs_bias_lstm, na.rm = TRUE),
      mean_pct_bias_abc  = mean(pct_bias_abc,  na.rm = TRUE),
      mean_pct_bias_lstm = mean(pct_bias_lstm, na.rm = TRUE),
      se_bias_abc        = sd(bias_abc,      na.rm = TRUE) / sqrt(n()),
      se_bias_lstm       = sd(bias_lstm,     na.rm = TRUE) / sqrt(n()),
      se_pct_bias_abc    = sd(pct_bias_abc,  na.rm = TRUE) / sqrt(n()),
      se_pct_bias_lstm   = sd(pct_bias_lstm, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    left_join(coverage_by_day, by = "date")
  
  # Overall stats
  overall <- bias_cov %>%
    summarize(
      mean_bias_abc      = mean(bias_abc,      na.rm = TRUE),
      median_bias_abc    = median(bias_abc,    na.rm = TRUE),
      rmse_abc           = sqrt(mean(bias_abc^2,  na.rm = TRUE)),
      mean_abs_bias_abc  = mean(abs_bias_abc,  na.rm = TRUE),
      mean_bias_lstm     = mean(bias_lstm,     na.rm = TRUE),
      median_bias_lstm   = median(bias_lstm,   na.rm = TRUE),
      rmse_lstm          = sqrt(mean(bias_lstm^2, na.rm = TRUE)),
      mean_abs_bias_lstm = mean(abs_bias_lstm, na.rm = TRUE),
      n_observations     = n()
    )
  
  # Overall coverage
  overall_cov <- all_results %>%
    select(sim_id, date, param_type, value = all_of(value_col)) %>%
    pivot_wider(names_from = param_type, values_from = value) %>%
    summarize(
      abc_q025      = quantile(abc,  0.025, na.rm = TRUE),
      abc_q975      = quantile(abc,  0.975, na.rm = TRUE),
      lstm_q025     = quantile(lstm, 0.025, na.rm = TRUE),
      lstm_q975     = quantile(lstm, 0.975, na.rm = TRUE),
      coverage_abc  = mean((true >= abc_q025)  & (true <= abc_q975),  na.rm = TRUE) * 100,
      coverage_lstm = mean((true >= lstm_q025) & (true <= lstm_q975), na.rm = TRUE) * 100
    )
  
  overall <- bind_cols(overall,
                       overall_cov %>% select(coverage_abc, coverage_lstm))
  
  list(bias_cov = bias_cov, daily = daily, overall = overall)
}

# ── Run for both infected prevalence and incidence ────────────────────────────
bc_infected  <- compute_bias_coverage(all_infected,
                                      value_col = "mean_infected",
                                      q025_col  = "q025",
                                      q975_col  = "q975")

bc_incidence <- compute_bias_coverage(all_incidence,
                                      value_col = "mean_incidence",
                                      q025_col  = "q025_incidence",
                                      q975_col  = "q975_incidence")

cat("\n=== INFECTED PREVALENCE SUMMARY ===\n"); print(bc_infected$overall)
cat("\n=== DAILY INCIDENCE (S→I) SUMMARY ===\n"); print(bc_incidence$overall)

# ============================================================
# PLOTTING HELPER — reusable for both datasets
# ============================================================
make_plots <- function(daily_stats, overall_stats, label = "Infected") {
  
  # ── Plot 1: Mean Bias ───────────────────────────────────────────────────────
  bias_long <- daily_stats %>%
    select(date, mean_bias_abc, mean_bias_lstm, se_bias_abc, se_bias_lstm) %>%
    pivot_longer(cols = c(mean_bias_abc, mean_bias_lstm),
                 names_to = "method", values_to = "bias") %>%
    mutate(
      method = case_when(method == "mean_bias_abc"  ~ "ABC",
                         method == "mean_bias_lstm" ~ "LSTM"),
      se = case_when(method == "ABC"  ~ se_bias_abc,
                     method == "LSTM" ~ se_bias_lstm)
    ) %>%
    select(-se_bias_abc, -se_bias_lstm)
  
  p1 <- ggplot(bias_long, aes(x = date, y = bias, color = method)) +
    geom_ribbon(aes(ymin = bias - 1.96*se, ymax = bias + 1.96*se, fill = method),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title    = paste("Mean Bias Over Time —", label),
         subtitle = "Bias = Predicted - True, with 95% confidence bands",
         x = "Day", y = "Mean Bias", color = "Method", fill = "Method") +
    theme_minimal() + theme(legend.position = "bottom") +
    scale_color_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
    scale_fill_manual(values  = c("ABC" = "#E74C3C", "LSTM" = "#3498DB"))
  
  # ── Plot 2: Coverage ────────────────────────────────────────────────────────
  cov_long <- daily_stats %>%
    select(date, coverage_abc, coverage_lstm) %>%
    pivot_longer(cols = c(coverage_abc, coverage_lstm),
                 names_to = "method", values_to = "coverage") %>%
    mutate(method = case_when(method == "coverage_abc"  ~ "ABC",
                              method == "coverage_lstm" ~ "LSTM"))
  
  p2 <- ggplot(cov_long, aes(x = date, y = coverage, color = method)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "gray50") +
    labs(title    = paste("Coverage Over Time —", label),
         subtitle = "% of true values within predicted 95% CI",
         x = "Day", y = "Coverage (%)", color = "Method") +
    theme_minimal() + theme(legend.position = "bottom") +
    scale_color_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
    ylim(0, 100)
  
  # ── Plot 3: Percentage Bias ─────────────────────────────────────────────────
  pct_long <- daily_stats %>%
    select(date, mean_pct_bias_abc, mean_pct_bias_lstm,
           se_pct_bias_abc, se_pct_bias_lstm) %>%
    pivot_longer(cols = c(mean_pct_bias_abc, mean_pct_bias_lstm),
                 names_to = "method", values_to = "pct_bias") %>%
    mutate(
      method = case_when(method == "mean_pct_bias_abc"  ~ "ABC",
                         method == "mean_pct_bias_lstm" ~ "LSTM"),
      se = case_when(method == "ABC"  ~ se_pct_bias_abc,
                     method == "LSTM" ~ se_pct_bias_lstm)
    ) %>%
    select(-se_pct_bias_abc, -se_pct_bias_lstm)
  
  p3 <- ggplot(pct_long, aes(x = date, y = pct_bias, color = method)) +
    geom_ribbon(aes(ymin = pct_bias - 1.96*se, ymax = pct_bias + 1.96*se, fill = method),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = paste("Percentage Bias Over Time —", label),
         x = "Day", y = "Mean Percentage Bias (%)",
         color = "Method", fill = "Method") +
    theme_minimal() + theme(legend.position = "bottom") +
    scale_color_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
    scale_fill_manual(values  = c("ABC" = "#E74C3C", "LSTM" = "#3498DB"))
  
  # ── Plot 4: Summary bar charts ──────────────────────────────────────────────
  summary_table <- data.frame(
    Method        = c("ABC", "LSTM"),
    Mean_Bias     = c(overall_stats$mean_bias_abc,      overall_stats$mean_bias_lstm),
    RMSE          = c(overall_stats$rmse_abc,            overall_stats$rmse_lstm),
    Mean_Abs_Bias = c(overall_stats$mean_abs_bias_abc,   overall_stats$mean_abs_bias_lstm),
    Coverage      = c(overall_stats$coverage_abc,        overall_stats$coverage_lstm)
  )
  
  p4 <- summary_table %>%
    select(-Coverage) %>%
    pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value") %>%
    ggplot(aes(x = Method, y = Value, fill = Method)) +
    geom_col(alpha = 0.7) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(title = paste("Summary Statistics —", label), x = "Method", y = "Value") +
    theme_minimal() + theme(legend.position = "none") +
    scale_fill_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB"))
  
  p5 <- ggplot(summary_table, aes(x = Method, y = Coverage, fill = Method)) +
    geom_col(alpha = 0.7) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "gray50") +
    labs(title    = paste("Overall Coverage —", label),
         subtitle = "Target: 95%",
         x = "Method", y = "Coverage (%)") +
    theme_minimal() + theme(legend.position = "none") +
    scale_fill_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
    ylim(0, 100)
  
  list(p_bias = p1, p_coverage = p2, p_pct_bias = p3,
       p_summary = p4, p_coverage_summary = p5,
       summary_table = summary_table)
}

# ── Generate plots for both outcomes ─────────────────────────────────────────
plots_infected  <- make_plots(bc_infected$daily,  bc_infected$overall,  "Infected Prevalence")
plots_incidence <- make_plots(bc_incidence$daily, bc_incidence$overall, "Daily Incidence (S→I)")

# ── Combined final figures ────────────────────────────────────────────────────
final_infected <- (plots_infected$p_pct_bias  | plots_infected$p_coverage) /
  (plots_infected$p_summary   | plots_infected$p_coverage_summary) +
  plot_annotation(
    title    = "Bias & Coverage — Infected Prevalence: ABC vs LSTM",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days",
    theme    = theme(plot.title = element_text(size = 16, face = "bold"))
  )

final_incidence <- (plots_incidence$p_pct_bias  | plots_incidence$p_coverage) /
  (plots_incidence$p_summary   | plots_incidence$p_coverage_summary) +
  plot_annotation(
    title    = "Bias & Coverage — Daily Incidence (S→I): ABC vs LSTM",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days",
    theme    = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(final_infected)
print(final_incidence)

ggsave("~/sima/Calibrate_ABM/figures/final_plot_infected.png",  final_infected,  width = 12, height = 10, dpi = 300)
ggsave("~/sima/Calibrate_ABM/figures/final_plot_incidence.png", final_incidence, width = 12, height = 10, dpi = 300)

cat("\n=== INFECTED PREVALENCE SUMMARY ===\n")
print(plots_infected$summary_table)
cat("\n=== DAILY INCIDENCE (S→I) SUMMARY ===\n")
print(plots_incidence$summary_table)

