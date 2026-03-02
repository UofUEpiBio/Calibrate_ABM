# ───────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(epiworldR)
library(glue)
source("~/sima/Calibrate_ABM/01a-bilstm.R")

# ── 1) LOAD ABC PREDICTIONS ───────────────────────────────────────────────────

filename <- paste0("~/sima/Calibrate_ABM/data-result/abc_parameters.csv")
ndays    <- 60
nsims    <- 1
nthreads <- 10
abc_pred <- read_csv(filename)

# ── 2) RUN LSTM CALIBRATION LOOP ─────────────────────────────────────────────

lstm_results <- vector("list", nrow(abc_pred))

for (i in 1:nrow(abc_pred)) {
  sim_id <- abc_pred$sim_id[i]
  cat(sprintf("Processing sim_id %d (%d/%d)\n", sim_id, i, nrow(abc_pred)))
  
  tryCatch({
    m0 <- ModelSIRCONN(
      name              = paste0("true_sim_", sim_id),
      n                 = abc_pred$true_n[i],
      prevalence        = abc_pred$true_preval[i],
      contact_rate      = abc_pred$true_crate[i],
      transmission_rate = abc_pred$true_ptran[i],
      recovery_rate     = abc_pred$true_recov[i]
    )
    run(m0, ndays = 60)
    
    incidence  <- plot_incidence(m0, plot = FALSE)
    all_counts <- incidence[, 1]
    ts         <- all_counts
    
    if (length(ts) != ndays + 1) {
      stop(glue("After dropping day0 and subsetting, length(ts) = {length(ts)}, not {ndays}!"))
    }
    
    lstm_out   <- calibrate_sir(ts, abc_pred$true_n[i], abc_pred$true_recov[i])
    ptran_lstm <- lstm_out[1]
    crate_lstm <- lstm_out[2]
    R0_lstm    <- lstm_out[3]
    
    lstm_results[[i]] <- data.frame(
      sim_id                 = sim_id,
      prevalence_lstm        = abc_pred$true_preval[i],
      contact_rate_lstm      = crate_lstm,
      transmission_rate_lstm = ptran_lstm,
      R0_lstm                = R0_lstm,
      stringsAsFactors       = FALSE
    )
    cat(sprintf("  ✓ Success for sim_id %d\n", sim_id))
    
  }, error = function(e) {
    cat(sprintf("  ✗ Error for sim_id %d: %s\n", sim_id, e$message))
    lstm_results[[i]] <<- data.frame(
      sim_id                 = sim_id,
      prevalence_lstm        = NA,
      contact_rate_lstm      = NA,
      transmission_rate_lstm = NA,
      R0_lstm                = NA,
      stringsAsFactors       = FALSE
    )
  })
}

# ── 3) COMBINE & FILTER LSTM RESULTS ─────────────────────────────────────────

lstm_preds <- bind_rows(lstm_results) %>%
  filter(!is.na(prevalence_lstm))

cat(sprintf("\nSummary:\n"))
cat(sprintf("Total simulations:          %d\n", nrow(abc_pred)))
cat(sprintf("Successful LSTM predictions:%d\n", nrow(lstm_preds)))
cat(sprintf("Success rate:               %.2f%%\n", 100 * nrow(lstm_preds) / nrow(abc_pred)))

# ── 4) DUPLICATE DIAGNOSTICS ─────────────────────────────────────────────────

cat("\n── Duplicate sim_id check ──\n")

dups_abc  <- abc_pred   %>% count(sim_id) %>% filter(n > 1)
dups_lstm <- lstm_preds %>% count(sim_id) %>% filter(n > 1)

if (nrow(dups_abc) > 0) {
  cat(sprintf("⚠ Duplicate sim_ids in abc_pred (%d affected):\n", nrow(dups_abc)))
  print(dups_abc)
} else {
  cat("✓ No duplicate sim_ids in abc_pred\n")
}

if (nrow(dups_lstm) > 0) {
  cat(sprintf("⚠ Duplicate sim_ids in lstm_preds (%d affected):\n", nrow(dups_lstm)))
  print(dups_lstm)
} else {
  cat("✓ No duplicate sim_ids in lstm_preds\n")
}

# ── 5) DEDUPLICATE ────────────────────────────────────────────────────────────

abc_pred   <- abc_pred   %>% distinct(sim_id, .keep_all = TRUE)
lstm_preds <- lstm_preds %>% distinct(sim_id, .keep_all = TRUE)

cat(sprintf("\nAfter deduplication:\n"))
cat(sprintf("  abc_pred rows:   %d\n", nrow(abc_pred)))
cat(sprintf("  lstm_preds rows: %d\n", nrow(lstm_preds)))

# ── 6) BUILD params_all ───────────────────────────────────────────────────────

params_all <- abc_pred %>%
  select(sim_id,
         true_n,
         prevalence_true        = true_preval,
         contact_rate_true      = true_crate,
         transmission_rate_true = true_ptran,
         recovery_rate_true     = true_recov,
         R0_true                = true_R0,
         R0_abc                 = abc_R0,
         prevalence_abc         = abc_preval,
         contact_rate_abc       = abc_crate,
         transmission_rate_abc  = abc_ptran,
         recovery_rate_abc      = abc_recov) %>%
  inner_join(lstm_preds, by = "sim_id") %>%
  mutate(recovery_rate_lstm = recovery_rate_true) %>%
  pivot_longer(
    cols          = -c(sim_id, true_n),
    names_to      = c(".value", "param_type"),
    names_pattern = "(.*)_(true|abc|lstm)"
  )

cat(sprintf("\nparams_all rows after pivot: %d\n", nrow(params_all)))
cat(sprintf("param_type levels: %s\n", paste(unique(params_all$param_type), collapse = ", ")))

write.csv(params_all, "~/sima/Calibrate_ABM/data-result/params_all.csv", row.names = FALSE)

# ── 7) FORWARD SIMULATIONS ────────────────────────────────────────────────────

cat("\nRunning forward simulations...\n")

all_sims <- params_all %>%
  mutate(run_id = row_number()) %>%
  pmap_dfr(function(run_id, sim_id, true_n, param_type,
                    prevalence, contact_rate,
                    transmission_rate, recovery_rate, R0) {
    
    calibrated_model <- ModelSIRCONN(
      name              = paste0("run_", run_id),
      n                 = true_n,
      prevalence        = prevalence,
      contact_rate      = contact_rate,
      transmission_rate = transmission_rate,
      recovery_rate     = recovery_rate
    )
    
    run(calibrated_model, ndays = 60)
    
    incidence_calib     <- plot_incidence(calibrated_model, plot = FALSE)
    infected_calibrated <- incidence_calib[, 1]

    
    tibble(
      time       = 1:length(infected_calibrated),
      I          = infected_calibrated,
      rep        = 1,
      param_type = param_type,
      sim_id     = sim_id
    )
  })

# ── 8) SUMMARISE ──────────────────────────────────────────────────────────────

summary_df <- all_sims %>%
  group_by(param_type, time) %>%
  summarise(
    median_I = median(I),
    lower95  = quantile(I, 0.025),
    upper95  = quantile(I, 0.975),
    .groups  = "drop"
  )

cat("\nPeak incidence per param_type:\n")
summary_df %>%
  group_by(param_type) %>%
  summarise(
    peak_median = max(median_I),
    peak_day    = time[which.max(median_I)],
    .groups     = "drop"
  ) %>%
  print()

# ── 9) PLOT ───────────────────────────────────────────────────────────────────

p <- ggplot(summary_df,
            aes(x = time, y = median_I,
                color = param_type, fill = param_type)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(
    x     = "Day",
    y     = "Active Infected (I)",
    title = "True vs. ABC vs. LSTM-Estimated SIRCONN Trajectories",
    color = "Parameter set",
    fill  = "Parameter set"
  ) +
  theme_minimal(base_size = 14)

print(p)

ggsave("~/sima/Calibrate_ABM/figures/comparison_plot.png",
       plot = p, width = 10, height = 8, dpi = 300)
cat("\nPlot saved as 'comparison_plot.png'\n")

# ── 10) SAVE ALL RESULTS ──────────────────────────────────────────────────────

results_summary <- list(
  abc_predictions     = abc_pred,
  lstm_predictions    = lstm_preds,
  combined_parameters = params_all,
  simulation_results  = all_sims,
  summary_statistics  = summary_df
)

saveRDS(results_summary,
        "~/sima/Calibrate_ABM/data-result/comparison_results.rds")
cat("Results saved as 'comparison_results.rds'\n")

# ── 11) BIAS CALCULATIONS ─────────────────────────────────────────────────────

cat("\nCalculating parameter biases...\n")

true_params <- params_all %>%
  filter(param_type == "true") %>%
  select(sim_id,
         prevalence_true        = prevalence,
         contact_rate_true      = contact_rate,
         transmission_rate_true = transmission_rate,
         recovery_rate_true     = recovery_rate,
         R0_true                = R0)

bias_df <- params_all %>%
  filter(param_type != "true") %>%
  left_join(true_params, by = "sim_id") %>%
  mutate(
    bias_prevalence        = prevalence        - prevalence_true,
    bias_contact_rate      = contact_rate      - contact_rate_true,
    bias_transmission_rate = transmission_rate - transmission_rate_true,
    bias_recovery_rate     = recovery_rate     - recovery_rate_true,
    bias_R0                = R0                - R0_true,
    rel_bias_contact_rate      = 100 * bias_contact_rate      / contact_rate_true,
    rel_bias_transmission_rate = 100 * bias_transmission_rate / transmission_rate_true,
    rel_bias_R0                = 100 * bias_R0                / R0_true
  )

bias_summary <- bias_df %>%
  group_by(param_type) %>%
  summarise(
    mean_bias_contact_rate          = mean(bias_contact_rate,      na.rm = TRUE),
    mean_bias_transmission_rate     = mean(bias_transmission_rate, na.rm = TRUE),
    mean_bias_R0                    = mean(bias_R0,                na.rm = TRUE),
    mean_rel_bias_contact_rate      = mean(rel_bias_contact_rate,      na.rm = TRUE),
    mean_rel_bias_transmission_rate = mean(rel_bias_transmission_rate, na.rm = TRUE),
    mean_rel_bias_R0                = mean(rel_bias_R0,               na.rm = TRUE),
    rmse_contact_rate               = sqrt(mean(bias_contact_rate^2,      na.rm = TRUE)),
    rmse_transmission_rate          = sqrt(mean(bias_transmission_rate^2, na.rm = TRUE)),
    rmse_R0                         = sqrt(mean(bias_R0^2,               na.rm = TRUE)),
    .groups = "drop"
  )

cat("\nBias summary:\n")
print(bias_summary)

write.csv(bias_df,      "~/sima/Calibrate_ABM/data-result/bias_df.csv",      row.names = FALSE)
write.csv(bias_summary, "~/sima/Calibrate_ABM/data-result/bias_summary.csv", row.names = FALSE)
cat("\nBias tables saved.\n")

# ── 12) BIAS PLOT ─────────────────────────────────────────────────────────────

bias_long <- bias_df %>%
  select(sim_id, param_type,
         contact_rate      = bias_contact_rate,
         transmission_rate = bias_transmission_rate,
         R0                = bias_R0) %>%
  pivot_longer(cols      = -c(sim_id, param_type),
               names_to  = "parameter",
               values_to = "bias")

p_bias <- ggplot(bias_long,
                 aes(x = param_type, y = bias, fill = param_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(
    x     = "Method",
    y     = "Bias (estimated − true)",
    title = "Parameter Bias: ABC vs. LSTM",
    fill  = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(p_bias)

ggsave("~/sima/Calibrate_ABM/figures/bias_plot.png",
       plot = p_bias, width = 10, height = 6, dpi = 300)
cat("Bias plot saved as 'bias_plot.png'\n")

# ── end script ────────────────────────────────────────────────────────────────
