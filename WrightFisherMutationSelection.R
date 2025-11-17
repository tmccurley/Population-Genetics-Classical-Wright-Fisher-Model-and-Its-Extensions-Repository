################################################################################
# Wright–Fisher Simulations in R
# (Neutral and Mutation Only — Selection Removed)
# Author: Thomas McCurley Project
################################################################################

# Setup ------------------------------------------------------------------------
rm(list = ls())
set.seed(20251030)
dir.create("results", showWarnings = FALSE)

library(ggplot2)

# Helper functions -------------------------------------------------------------

# Simulate neutral WF trajectories until absorption
simulate_wf_neutral <- function(N, p0, nTraj, max_gen = 10000) {
  p_mat <- matrix(NA, nrow = max_gen + 1, ncol = nTraj)
  p_mat[1, ] <- p0
  abs_time <- rep(NA, nTraj)
  fix_state <- rep(NA, nTraj)
  
  for (j in 1:nTraj) {
    p <- p0
    for (t in 1:max_gen) {
      if (p == 0 || p == 1) {
        abs_time[j] <- t - 1
        fix_state[j] <- p
        p_mat[(t + 1):(max_gen + 1), j] <- p
        break
      }
      x <- rbinom(1, 2 * N, p)
      p <- x / (2 * N)
      p_mat[t + 1, j] <- p
    }
  }
  list(p_mat = p_mat, abs_time = abs_time, fix_state = fix_state)
}

# Simulate WF with mutation
simulate_wf_mutation <- function(N, p0, mu, nu, nTraj, T) {
  p_mat <- matrix(NA, nrow = T + 1, ncol = nTraj)
  p_mat[1, ] <- p0
  for (j in 1:nTraj) {
    p <- p0
    for (t in 1:T) {
      q <- p * (1 - mu) + (1 - p) * nu
      x <- rbinom(1, 2 * N, q)
      p <- x / (2 * N)
      p_mat[t + 1, j] <- p
    }
  }
  p_mat
}

# Theoretical variance (neutral)
theoretical_variance <- function(p0, N, T) {
  t_vals <- 0:T
  p0 * (1 - p0) * (1 - (1 - 1 / (2 * N))^t_vals)
}

# Plot helper
save_plot <- function(plot, filename) {
  ggsave(filename, plot, width = 7, height = 4, dpi = 300)
}

# ------------------------------------------------------------------------------

# PARAMETERS -------------------------------------------------------------------
N <- 300
p0 <- 0.3
nTraj <- 500

# ========================
# 1) NEUTRAL MODEL
# ========================
neutral_res <- simulate_wf_neutral(N, p0, nTraj, max_gen = 10000)
p_mat <- neutral_res$p_mat
T <- nrow(p_mat) - 1

# Empirical mean/variance
mean_p <- rowMeans(p_mat, na.rm = TRUE)
var_p <- apply(p_mat, 1, var, na.rm = TRUE)
var_theory <- theoretical_variance(p0, N, T)

# Absorption info
absorbed <- !is.na(neutral_res$fix_state)
p1_est <- mean(neutral_res$fix_state[absorbed] == 1)
p0_est <- mean(neutral_res$fix_state[absorbed] == 0)
avg_time_fix1 <- mean(neutral_res$abs_time[neutral_res$fix_state == 1], na.rm = TRUE)
avg_time_fix0 <- mean(neutral_res$abs_time[neutral_res$fix_state == 0], na.rm = TRUE)

# Plot trajectories
traj_df <- data.frame(gen = rep(0:T, nTraj),
                      freq = as.vector(p_mat),
                      traj = rep(1:nTraj, each = T + 1))
neutral_traj_plot <- ggplot(traj_df[traj_df$traj <= 50, ], aes(x = gen, y = freq, group = traj)) +
  geom_line(alpha = 0.5, color = "black") +
  geom_hline(yintercept = p0, color = "red", linetype = "dashed") +
  labs(title = paste0("Neutral Wright–Fisher trajectories (N=", N, ")"),
       x = "Generation", y = "Allele frequency p_t")
save_plot(neutral_traj_plot, "results/neutral_trajectories.png")

# Variance plot
var_df <- data.frame(
  generation = 0:T,
  empirical = var_p,
  theoretical = var_theory
)
neutral_var_plot <- ggplot(var_df, aes(x = generation)) +
  geom_line(aes(y = empirical), color = "black", linewidth = 1) +
  geom_line(aes(y = theoretical), color = "blue", linetype = "dashed", linewidth = 1) +
  labs(title = "Neutral Model: Empirical vs Theoretical Variance",
       x = "Generation", y = "Var(p_t)")
save_plot(neutral_var_plot, "results/neutral_variance.png")

# Absorption barplot
abs_df <- data.frame(State = c("Loss (0)", "Fixation (1)"),
                     Proportion = c(p0_est, p1_est))
neutral_abs_plot <- ggplot(abs_df, aes(x = State, y = Proportion, fill = State)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("steelblue", "orange")) +
  labs(title = "Neutral Model: Absorption Outcomes", y = "Proportion") +
  theme(legend.position = "none")
save_plot(neutral_abs_plot, "results/neutral_absorption_props.png")

# Absorption time histogram
time_df <- data.frame(time = neutral_res$abs_time[absorbed])
neutral_time_plot <- ggplot(time_df, aes(x = time)) +
  geom_histogram(fill = "lightgreen", bins = 30, color = "white") +
  labs(title = "Neutral Model: Absorption Times", x = "Generation", y = "Count")
save_plot(neutral_time_plot, "results/neutral_absorption_times.png")

# ========================
# 2) MUTATION MODEL
# ========================
mu <- 0.002
nu <- 0.0005
T_mut <- 3000

p_mat_mut <- simulate_wf_mutation(N, p0, mu, nu, nTraj, T_mut)
mean_p_mut <- rowMeans(p_mat_mut)
var_p_mut <- apply(p_mat_mut, 1, var)

p_eq <- nu / (mu + nu)
alpha_param <- 2 * N * nu
beta_param <- 2 * N * mu

# Trajectories (mutation)
traj_df_mut <- data.frame(gen = rep(0:T_mut, nTraj),
                          freq = as.vector(p_mat_mut),
                          traj = rep(1:nTraj, each = T_mut + 1))
mut_traj_plot <- ggplot(traj_df_mut[traj_df_mut$traj <= 50, ], aes(x = gen, y = freq, group = traj)) +
  geom_line(alpha = 0.4, color = "gray30") +
  geom_hline(yintercept = p_eq, color = "green4", linetype = "dashed") +
  labs(title = "Wright–Fisher with Mutation", x = "Generation", y = "Allele frequency p_t")
save_plot(mut_traj_plot, "results/mut_trajectories.png")

# Final distribution
final_freqs <- p_mat_mut[T_mut + 1, ]
final_df <- data.frame(freq = final_freqs)
x_vals <- seq(0, 1, length.out = 500)
beta_density <- dbeta(x_vals, alpha_param, beta_param)
beta_df <- data.frame(x = x_vals, y = beta_density)

mut_hist_plot <- ggplot(final_df, aes(x = freq)) +
  geom_histogram(aes(y = ..density..), fill = "skyblue", color = "white", bins = 30, alpha = 0.6) +
  geom_line(data = beta_df, aes(x = x, y = y), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Mutation Model: Stationary Distribution", x = "Allele frequency", y = "Density")
save_plot(mut_hist_plot, "results/mut_final_hist_beta.png")

################################################################################
cat("✅ All simulations (neutral + mutation) complete.\nFiles saved in 'results/' folder.\n")
################################################################################
