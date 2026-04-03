# =============================================================================
# Interactive Session 2: Random Search
# =============================================================================
# Lotka-Volterra parameter search via uniform random sampling and ODE SSE loss.
# Run from repository root:
#   Rscript vignettes/LV_random_search.R
#
# Requires working directory such that R/ODE_SSE_loss.R and data/ are reachable
# (typically the repo root).

resolve_repo_root <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", ca[grep("^--file=", ca)])
    if (length(f)) {
        return(normalizePath(file.path(dirname(normalizePath(f[1])), "..")))
    }
    if (file.exists(file.path(getwd(), "R", "ODE_SSE_loss.R"))) {
        return(normalizePath(getwd()))
    }
    up <- normalizePath(file.path(getwd(), ".."))
    if (file.exists(file.path(up, "R", "ODE_SSE_loss.R"))) {
        return(up)
    }
    normalizePath(getwd())
}

root <- resolve_repo_root()
source(file.path(root, "R", "ODE_SSE_loss.R"))

# ----- 1. Uniform priors: alpha, beta, gamma, delta ~ Unif(prior_min[k], prior_max[k]) -----
prior_min <- c(alpha = 0, beta = 0, gamma = 0, delta = 0)
prior_max <- c(alpha = 2, beta = 2, gamma = 2, delta = 2)
stopifnot(all(prior_max > prior_min))

# ----- 2. Number of trials -----
N <- 10^4 * 3
stopifnot(N >= 1L)

# ----- Data & metadata -----
data_path <- file.path(root, "data", "lv_ssa_fake_counts.csv")
meta_path <- file.path(root, "data", "lv_ssa_metadata.csv")
stopifnot(file.exists(data_path), file.exists(meta_path))
obs  <- read.csv(data_path, check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)

set.seed(1)

# ----- 3. Random search -----
t0 <- proc.time()
trials <- data.frame(
    alpha  = numeric(N),
    beta   = numeric(N),
    gamma  = numeric(N),
    delta  = numeric(N),
    loss   = numeric(N)
)
nm <- c("alpha", "beta", "gamma", "delta")
for (i in seq_len(N)) {
    par <- stats::setNames(
        mapply(function(lo, hi) stats::runif(1L, min = lo, max = hi), prior_min, prior_max),
        nm
    )
    trials[i, nm] <- unlist(par)
    trials$loss[i] <- ODE_SSE_loss(as.list(trials[i, nm]), obs, metadata = meta)
}
runtime_sec <- as.numeric((proc.time() - t0)["elapsed"])

# ----- 4. Save all trials; report best -----
out_dir <- file.path(root, "vignettes", "random_search")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
trials_path <- file.path(out_dir, "lv_random_search_trials.csv")
write.csv(trials, file = trials_path, row.names = FALSE)

i_best <- which.min(trials$loss)
best <- trials[i_best, , drop = FALSE]
message(
    "Best trial (#", i_best, "/", N, "): ",
    "alpha=", round(best$alpha, 6),
    ", beta=", round(best$beta, 6),
    ", gamma=", round(best$gamma, 6),
    ", delta=", round(best$delta, 6),
    ", loss=", signif(best$loss, 8),
    ". All trials: ", trials_path
)

# ----- 5. Refit best parameters & plot (ggplot style ~ temp_scripts/plot_lv_ssa_fake_counts.R) -----
library(ggplot2)

dat <- obs[, c("time", "rabbit", "fox")]
dat <- dat[order(dat$time), , drop = FALSE]
i0 <- which.min(abs(dat$time - 0))
stopifnot(abs(dat$time[i0]) < 1e-8)
initial_state <- list(
    rabbit = as.numeric(dat$rabbit[i0]),
    fox    = as.numeric(dat$fox[i0])
)
best_list <- stats::setNames(as.numeric(best[1, nm]), nm)
best_params <- as.list(best_list)
best_params$N <- as.numeric(meta$N[1])
pred <- lv_populations_at_time(best_params, dat$time, initial_state)

# Observed: identical to temp_scripts/plot_lv_ssa_fake_counts.R (points only).
# Fit: lines only, same yellow (rabbit) / red (fox) as that script.
col_rabbit <- "#EFC000"
col_fox <- "#BC3C29"

obs_df <- rbind(
    data.frame(time = dat$time, species = "Rabbit", count = dat$rabbit),
    data.frame(time = dat$time, species = "Fox", count = dat$fox)
)
obs_df$species <- factor(obs_df$species, levels = c("Rabbit", "Fox"))

fit_df <- rbind(
    data.frame(time = pred$time, species = "Rabbit", count = pred$rabbit),
    data.frame(time = pred$time, species = "Fox", count = pred$fox)
)
fit_df$species <- factor(fit_df$species, levels = c("Rabbit", "Fox"))

p_fit <- ggplot() +
    geom_line(
        data = fit_df,
        aes(x = time, y = count, colour = species, group = species),
        linewidth = 1.2
    ) +
    geom_point(
        data = obs_df,
        aes(x = time, y = count, colour = species),
        size = 7,
        stroke = 0.3
    ) +
    scale_colour_manual(
        values = c(Rabbit = col_rabbit, Fox = col_fox),
        name = NULL
    ) +
    labs(x = "Time", y = "Count") +
    theme_bw(base_size = 20) +
    theme(
        plot.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(
            fill = "white",
            colour = "grey35",
            linewidth = 0.5
        ),
        legend.margin = margin(6, 8, 6, 8)
    )

plot_path <- file.path(out_dir, "lv_random_search_fit.png")
ggsave(
    filename = plot_path,
    plot = p_fit,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
message("Saved fit plot: ", plot_path)

# =============================================================================
# Plots: point estimate of best parameters
# =============================================================================
param_names <- c("alpha", "beta", "gamma", "delta")

long_pe <- data.frame(
    parameter = factor(param_names, levels = param_names),
    value     = as.numeric(best[1, param_names])
)

p_pe <- ggplot(long_pe, aes(x = value)) +
    geom_vline(aes(xintercept = value), colour = "#08519C", linewidth = 1.2) +
    geom_point(aes(y = 0), colour = "#08519C", size = 4) +
    geom_text(
        aes(y = 0, label = sprintf("%.3f", value)),
        colour = "#08519C", vjust = -1.5, size = 4.5, fontface = "bold"
    ) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    coord_cartesian(xlim = c(0, 2)) +
    labs(x = NULL, y = "Point estimate") +
    theme_bw(base_size = 20) +
    theme(
        axis.text.y  = element_blank(),
        panel.grid   = element_blank()
    )

ggsave(file.path(out_dir, "lv_random_search_point_estimate.png"),
       p_pe, width = 10, height = 6, dpi = 300, bg = "white")
message("Saved point estimate plot to ", out_dir)
