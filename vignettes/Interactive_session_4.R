# =============================================================================
# Interactive Session 4: ABC Rejection
# =============================================================================
# Lotka-Volterra: ABC Rejection using tau-leaping SSE as summary statistic.
#
# Uses abc_rejection() from R/abc.R (thin wrapper over EasyABC::ABC_rejection)
# and tau_leaping_SSE_loss() from R/tau_leaping_SSE_loss.R to compute the
# summary statistic. Tau-leaping is much faster than exact SSA.
#
# Prior: Uniform(0, 2) for alpha, beta, gamma, delta (same as LV_MCMC.R).
# Summary statistic: SSE between observed data and one tau-leaping trajectory.
# Target summary statistic: 0 (a perfect match).
#
# Run from repository root:
#   Rscript vignettes/LV_abc_rejection.R

# ---- Resolve repository root ------------------------------------------------
resolve_repo_root <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", ca[grep("^--file=", ca)])
    if (length(f)) {
        return(normalizePath(file.path(dirname(normalizePath(f[1])), "..")))
    }
    if (file.exists(file.path(getwd(), "R", "abc.R"))) {
        return(normalizePath(getwd()))
    }
    up <- normalizePath(file.path(getwd(), ".."))
    if (file.exists(file.path(up, "R", "abc.R"))) {
        return(up)
    }
    normalizePath(getwd())
}

root <- resolve_repo_root()
source(file.path(root, "R", "abc.R"))
source(file.path(root, "R", "tau_leaping_SSE_loss.R"))

# ---- Data --------------------------------------------------------------------
data_path <- file.path(root, "data", "lv_ssa_fake_counts.csv")
meta_path <- file.path(root, "data", "lv_ssa_metadata.csv")
stopifnot(file.exists(data_path), file.exists(meta_path))
obs  <- read.csv(data_path, check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)

# ---- Parameter labels --------------------------------------------------------
par_names <- c("alpha", "beta", "gamma", "delta")
parameters_labels <- data.frame(parameter = par_names, stringsAsFactors = FALSE)

# ---- Priors (EasyABC format): Uniform(0, 2) for each parameter --------------
prior_distributions <- list(
    c("unif", 0, 2),
    c("unif", 0, 2),
    c("unif", 0, 2),
    c("unif", 0, 2)
)

# ---- Model function (EasyABC simulator contract) -----------------------------
# Input : data.frame `parameters` — one row per particle, columns = par_names
# Output: data.frame whose first columns are the parameters, followed by
#         summary-statistic columns (here: a single column "SSE").
nSimulations <- 0

model <- function(parameters, parallel = TRUE) {
    n <- nrow(parameters)
    results <- vector("list", n)
    for (i in seq_len(n)) {
        nSimulations <<- nSimulations + 1L
        parms <- as.list(parameters[i, , drop = TRUE])
        names(parms) <- par_names
        parms$N <- as.numeric(meta$N[1])

        sse <- tryCatch(
            tau_leaping_SSE_loss(parms, obs, metadata = meta), # nolint: object_usage_linter.
            error = function(e) Inf
        )

        results[[i]] <- data.frame(
            alpha = as.numeric(parms$alpha),
            beta  = as.numeric(parms$beta),
            gamma = as.numeric(parms$gamma),
            delta = as.numeric(parms$delta),
            SSE   = sse
        )
    }
    do.call(rbind, results)
}

# ---- Observed summary statistic: SSE = 0 (perfect match) --------------------
statistics_target <- data.frame(SSE = 0)

# ---- ABC Rejection settings --------------------------------------------------
nParticles         <- 18000
tolerance_quantile <- 500 / 18000

# ---- Run ---------------------------------------------------------------------
set.seed(1)
t_run <- proc.time()

abc_result <- abc_rejection(
    statistics_target   = statistics_target,
    model               = model,
    parameters_labels   = parameters_labels,
    prior_distributions = prior_distributions,
    nParticles          = nParticles,
    tolerance_quantile  = tolerance_quantile,
    progress_bar        = TRUE
)

runtime_sec <- as.numeric((proc.time() - t_run)["elapsed"])
message("ABC rejection done. Runtime: ", signif(runtime_sec, 4), " s")
message("Total simulations: ", nSimulations)

# ---- Save results ------------------------------------------------------------
out_dir <- file.path(root, "vignettes", "ABC_rejection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

accepted_params <- abc_result[["Iteration_1"]]$parameters
write.csv(
    accepted_params,
    file.path(out_dir, "lv_abc_rejection_params.csv"),
    row.names = FALSE
)
write.csv(
    abc_result[["Iteration_1"]]$reference,
    file.path(out_dir, "lv_abc_rejection_reference.csv"),
    row.names = FALSE
)
message("Saved results to ", out_dir)

# ---- Plots -------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)

    long_df <- do.call(rbind, lapply(par_names, function(p) {
        data.frame(parameter = p, value = accepted_params[[p]],
                   stringsAsFactors = FALSE)
    }))

    true_vals <- data.frame(
        parameter  = par_names,
        true_value = c(1, 1, 1, 1)
    )

    p_hist <- ggplot(long_df, aes(x = value)) +
        geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, colour = "white") +
        geom_vline(data = true_vals, aes(xintercept = true_value),
                   colour = "red", linetype = "dashed", linewidth = 0.8) +
        facet_wrap(~parameter, ncol = 2) +
        xlim(0, 2) +
        labs(
            title = paste0("ABC Rejection posterior (n = ", nParticles,
                           ", tol = ", tolerance_quantile, ")"),
            x = "Parameter value", y = "Count"
        ) +
        theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "lv_abc_rejection_hist.png"),
           p_hist, width = 9, height = 6, dpi = 150)

    message("Saved plots to ", out_dir)
} else {
    message("Package ggplot2 not installed; skipping plots.")
}

# ---- Summary -----------------------------------------------------------------
cat("\n=== ABC Rejection Summary ===\n")
cat("nParticles:          ", nParticles, "\n")
cat("Tolerance quantile:  ", tolerance_quantile, "\n")
cat("Accepted particles:  ", nrow(accepted_params), "\n")
cat("Runtime (sec):       ", signif(runtime_sec, 4), "\n")
cat("\nPosterior means:\n")
print(colMeans(accepted_params))
cat("\nPosterior medians:\n")
print(apply(accepted_params, 2, median))
cat("\nTrue values: alpha = 1, beta = 1, gamma = 1, delta = 1\n")

# =============================================================================
# Plots: posterior predictive population fit + posterior density
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))
source(file.path(root, "R", "tau_leaping.R"))

posterior <- abc_result[["Iteration_1"]]$reference
N_meta <- as.numeric(meta$N[1])

dat <- obs[, c("time", "rabbit", "fox")]
dat <- dat[order(dat$time), , drop = FALSE]
t_obs <- dat$time
i0 <- which.min(abs(t_obs - 0))
stopifnot(abs(t_obs[i0]) < 1e-8)
initial_state <- list(
    rabbit = as.integer(dat$rabbit[i0]),
    fox    = as.integer(dat$fox[i0])
)

reaction_propensities <- function(state, p) {
    with(c(state, p), c(
        alpha * rabbit,
        (beta  / N) * rabbit * fox,
        (delta / N) * rabbit * fox,
        gamma * fox
    ))
}
reaction_stoichiometries <- list(
    list(rabbit =  1L, fox =  0L),
    list(rabbit = -1L, fox =  0L),
    list(rabbit =  0L, fox =  1L),
    list(rabbit =  0L, fox = -1L)
)

n_draw <- nrow(posterior)
n_t    <- length(t_obs)
rabbit_mat <- matrix(NA_real_, nrow = n_draw, ncol = n_t)
fox_mat    <- matrix(NA_real_, nrow = n_draw, ncol = n_t)

message("Running tau-leaping for ", n_draw, " posterior draws ...")
for (i in seq_len(n_draw)) {
    parms_i <- list(
        alpha = as.numeric(posterior$alpha[i]),
        beta  = as.numeric(posterior$beta[i]),
        gamma = as.numeric(posterior$gamma[i]),
        delta = as.numeric(posterior$delta[i]),
        N     = N_meta
    )
    raw <- tau_leaping(
        initial_state            = initial_state,
        parameters               = parms_i,
        reaction_propensities    = reaction_propensities,
        reaction_stoichiometries = reaction_stoichiometries,
        time_points              = t_obs,
        tau                      = 0.01
    )
    rabbit_mat[i, ] <- as.numeric(raw[, "rabbit"])
    fox_mat[i, ]    <- as.numeric(raw[, "fox"])
}

rabbit_mean <- colMeans(rabbit_mat)
fox_mean    <- colMeans(fox_mat)
rabbit_lo <- apply(rabbit_mat, 2L, quantile, probs = 0.025, names = FALSE)
rabbit_hi <- apply(rabbit_mat, 2L, quantile, probs = 0.975, names = FALSE)
fox_lo    <- apply(fox_mat, 2L, quantile, probs = 0.025, names = FALSE)
fox_hi    <- apply(fox_mat, 2L, quantile, probs = 0.975, names = FALSE)

summ <- data.frame(
    time = t_obs, rabbit_mean = rabbit_mean, rabbit_lo = rabbit_lo,
    rabbit_hi = rabbit_hi, fox_mean = fox_mean, fox_lo = fox_lo, fox_hi = fox_hi
)

col_rabbit <- "#EFC000"
col_fox    <- "#BC3C29"

summ_all <- rbind(
    data.frame(species = "Rabbit", time = summ$time, lo = summ$rabbit_lo, hi = summ$rabbit_hi),
    data.frame(species = "Fox", time = summ$time, lo = summ$fox_lo, hi = summ$fox_hi)
)
summ_all$species <- factor(summ_all$species, levels = c("Rabbit", "Fox"))

obs_df <- rbind(
    data.frame(time = dat$time, species = "Rabbit", count = dat$rabbit),
    data.frame(time = dat$time, species = "Fox",    count = dat$fox)
)
obs_df$species <- factor(obs_df$species, levels = c("Rabbit", "Fox"))

fit_df <- rbind(
    data.frame(time = summ$time, species = "Rabbit", count = summ$rabbit_mean),
    data.frame(time = summ$time, species = "Fox",    count = summ$fox_mean)
)
fit_df$species <- factor(fit_df$species, levels = c("Rabbit", "Fox"))

p_pop <- ggplot() +
    geom_ribbon(
        data = summ_all,
        aes(x = time, ymin = lo, ymax = hi, fill = species),
        alpha = 0.2
    ) +
    geom_line(
        data = fit_df,
        aes(x = time, y = count, colour = species, group = species),
        linewidth = 1.2
    ) +
    geom_point(
        data = obs_df,
        aes(x = time, y = count, colour = species),
        size = 7, stroke = 0.3
    ) +
    scale_fill_manual(values = c(Rabbit = col_rabbit, Fox = col_fox), guide = "none") +
    scale_colour_manual(values = c(Rabbit = col_rabbit, Fox = col_fox), name = NULL) +
    labs(x = "Time", y = "Count") +
    theme_bw(base_size = 20) +
    theme(
        plot.title = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "grey35", linewidth = 0.5),
        legend.margin = margin(6, 8, 6, 8)
    )

ggsave(file.path(out_dir, "lv_abc_rejection_population_fit_ci.png"),
       p_pop, width = 10, height = 6, dpi = 300, bg = "white")

# ---- Posterior density ----
long_post <- do.call(rbind, lapply(par_names, function(pn) {
    data.frame(parameter = pn, value = posterior[[pn]], stringsAsFactors = FALSE)
}))
long_post$parameter <- factor(long_post$parameter, levels = par_names)

p_den <- ggplot(long_post, aes(x = value)) +
    geom_density(fill = "#3182BD", alpha = 0.45, colour = "#08519C", linewidth = 0.6) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    coord_cartesian(xlim = c(0, 2)) +
    labs(x = NULL, y = "Density") +
    theme_bw(base_size = 20)

ggsave(file.path(out_dir, "lv_abc_rejection_density.png"),
       p_den, width = 10, height = 6, dpi = 300, bg = "white")
message("Saved population fit + density plots to ", out_dir)
