setwd( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/R"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
files_sources <- list.files(pattern = "\\.[rR]$") #<<<<<<<<<<<<<<<<<<<<<
sapply(files_sources, source) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
setwd( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/vignettes"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
library(abcsmcrf)
library(truncnorm)
library(ggplot2)
data <- read.csv( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/data/data.csv"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
set.seed(1)
N_particles <- rep(2500, 5)
parameter_ranges <- data.frame(
    alpha = c(0, 2),
    beta = c(0, 2),
    gamma = c(0, 2),
    delta = c(0, 2),
    row.names = c("min", "max")
)
initial_state <- list(
    rabbit = 700,
    fox = 300
)
reaction_propensities <- function(state, p) {
    with(c(state, p), c(
        alpha * rabbit,
        (beta / 1000) * rabbit * fox,
        (delta / 1000) * rabbit * fox,
        gamma * fox
    ))
}
reaction_stoichiometries <- list(
    list(rabbit = 1L, fox = 0L),
    list(rabbit = -1L, fox = 0L),
    list(rabbit = 0L, fox = 1L),
    list(rabbit = 0L, fox = -1L)
)
statistics_target <- data.frame(
    matrix(c(data$rabbit, data$fox), nrow = 1)
)
colnames(statistics_target) <- c(paste0("rabbit_", data$time), paste0("fox_", data$time))
model <- function(parameters, parallel = TRUE) {
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("tau_leaping", "initial_state", "reaction_propensities", "reaction_stoichiometries", "data"))
        stats <- parLapply(
            cl = cl, 1:nrow(parameters),
            function(i) {
                candidate_populations <- tau_leaping(
                    initial_state = initial_state,
                    parameters = parameters[i, ],
                    reaction_propensities = reaction_propensities,
                    reaction_stoichiometries = reaction_stoichiometries,
                    time_points = data$time,
                    tau = 0.01, max_pop = 5000L
                )
                output <- as.data.frame(matrix(c(candidate_populations[, "rabbit"], candidate_populations[, "fox"]), nrow = 1))
                return(output)
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            candidate_populations <- tau_leaping(
                initial_state = initial_state,
                parameters = parameters[i, ],
                reaction_propensities = reaction_propensities,
                reaction_stoichiometries = reaction_stoichiometries,
                time_points = data$time,
                tau = 0.01, max_pop = 5000L
            )
            stats <- rbind(stats, c(candidate_populations[, "rabbit"], candidate_populations[, "fox"]))
        }
    }
    stats <- as.data.frame(stats)
    colnames(stats) <- c(paste0("rabbit_", data$time), paste0("fox_", data$time))
    return(cbind(parameters, stats))
}
dprior <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "alpha")) {
        probs <- probs * dunif(parameters[["alpha"]], min = parameter_ranges["min", "alpha"], max = parameter_ranges["max", "alpha"])
    }
    if (parameter_id %in% c("all", "beta")) {
        probs <- probs * dunif(parameters[["beta"]], min = parameter_ranges["min", "beta"], max = parameter_ranges["max", "beta"])
    }
    if (parameter_id %in% c("all", "gamma")) {
        probs <- probs * dunif(parameters[["gamma"]], min = parameter_ranges["min", "gamma"], max = parameter_ranges["max", "gamma"])
    }
    if (parameter_id %in% c("all", "delta")) {
        probs <- probs * dunif(parameters[["delta"]], min = parameter_ranges["min", "delta"], max = parameter_ranges["max", "delta"])
    }
    return(probs)
}
rprior <- function(Nparameters) {
    data.frame(
        alpha = runif(Nparameters, min = parameter_ranges["min", "alpha"], max = parameter_ranges["max", "alpha"]),
        beta  = runif(Nparameters, min = parameter_ranges["min", "beta"], max = parameter_ranges["max", "beta"]),
        gamma = runif(Nparameters, min = parameter_ranges["min", "gamma"], max = parameter_ranges["max", "gamma"]),
        delta = runif(Nparameters, min = parameter_ranges["min", "delta"], max = parameter_ranges["max", "delta"])
    )
}
dperturb <- function(parameters, parameters_previous, parameters_previous_sampled, iteration, parameter_id = "all") {
    Beaumont_variances <- 2 * pmax(sapply(parameters_previous_sampled, var), 1e-10)
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "alpha")) {
        probs <- probs * dnorm(parameters[["alpha"]], mean = parameters_previous[["alpha"]], sd = sqrt(Beaumont_variances[["alpha"]]))
    }
    if (parameter_id %in% c("all", "beta")) {
        probs <- probs * dnorm(parameters[["beta"]], mean = parameters_previous[["beta"]], sd = sqrt(Beaumont_variances[["beta"]]))
    }
    if (parameter_id %in% c("all", "gamma")) {
        probs <- probs * dnorm(parameters[["gamma"]], mean = parameters_previous[["gamma"]], sd = sqrt(Beaumont_variances[["gamma"]]))
    }
    if (parameter_id %in% c("all", "delta")) {
        probs <- probs * dnorm(parameters[["delta"]], mean = parameters_previous[["delta"]], sd = sqrt(Beaumont_variances[["delta"]]))
    }
    return(probs)
}
rperturb <- function(parameters_unperturbed, parameters_previous_sampled, iteration) {
    Beaumont_variances <- 2 * pmax(sapply(parameters_previous_sampled, var), 1e-10)
    parameters_perturbed <- parameters_unperturbed
    parameters_perturbed[["alpha"]] <- rtruncnorm(
        n = nrow(parameters_perturbed),
        a = parameter_ranges["min", "alpha"], b = parameter_ranges["max", "alpha"],
        mean = parameters_perturbed[["alpha"]],
        sd = sqrt(Beaumont_variances[["alpha"]])
    )
    parameters_perturbed[["beta"]] <- rtruncnorm(
        n = nrow(parameters_perturbed),
        a = parameter_ranges["min", "beta"], b = parameter_ranges["max", "beta"],
        mean = parameters_perturbed[["beta"]],
        sd = sqrt(Beaumont_variances[["beta"]])
    )
    parameters_perturbed[["gamma"]] <- rtruncnorm(
        n = nrow(parameters_perturbed),
        a = parameter_ranges["min", "gamma"], b = parameter_ranges["max", "gamma"],
        mean = parameters_perturbed[["gamma"]],
        sd = sqrt(Beaumont_variances[["gamma"]])
    )
    parameters_perturbed[["delta"]] <- rtruncnorm(
        n = nrow(parameters_perturbed),
        a = parameter_ranges["min", "delta"], b = parameter_ranges["max", "delta"],
        mean = parameters_perturbed[["delta"]],
        sd = sqrt(Beaumont_variances[["delta"]])
    )
    return(parameters_perturbed)
}
abcsmcrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    rperturb = rperturb,
    dperturb = dperturb,
    nParticles = N_particles,
    parallel = TRUE
)
accepted_parameters <- abcsmcrf_results[[paste0("Iteration_", length(N_particles) + 1)]]$parameters
predicted_populations <- lapply(1:nrow(accepted_parameters), function(i) {
    predicted_populations <- as.data.frame(tau_leaping(
        initial_state = initial_state,
        parameters = accepted_parameters[i, ],
        reaction_propensities = reaction_propensities,
        reaction_stoichiometries = reaction_stoichiometries,
        time_points = seq(0, 14, by = 0.1),
        tau = 0.01, max_pop = 5000L
    ))
    predicted_populations$time <- seq(0, 14, by = 0.1)
    predicted_populations$iteration <- i
    return(predicted_populations)
})
predicted_populations <- do.call(rbind, predicted_populations)
predicted_populations_mean <- aggregate(cbind(rabbit, fox) ~ time, data = predicted_populations, FUN = mean)
predicted_populations_CI_low <- aggregate(cbind(rabbit, fox) ~ time, data = predicted_populations, FUN = function(x) quantile(x, 0.025))
predicted_populations_CI_high <- aggregate(cbind(rabbit, fox) ~ time, data = predicted_populations, FUN = function(x) quantile(x, 0.975))
#-----------------------------------------------------------------------Plot fitted & observed population dynamics
color_rabbit <- "#EFC000"
color_fox <- "#BC3C29"
obs_df <- rbind(
    data.frame(time = data$time, species = "Rabbit", count = data$rabbit),
    data.frame(time = data$time, species = "Fox", count = data$fox)
)
obs_df$species <- factor(obs_df$species, levels = c("Rabbit", "Fox"))
fit_ranges_df <- rbind(
    data.frame(species = "Rabbit", time = predicted_populations_mean$time, lo = predicted_populations_CI_low$rabbit, hi = predicted_populations_CI_high$rabbit),
    data.frame(species = "Fox", time = predicted_populations_mean$time, lo = predicted_populations_CI_low$fox, hi = predicted_populations_CI_high$fox)
)
fit_ranges_df$species <- factor(fit_ranges_df$species, levels = c("Rabbit", "Fox"))
fit_df <- rbind(
    data.frame(time = predicted_populations_mean$time, species = "Rabbit", count = predicted_populations_mean$rabbit),
    data.frame(time = predicted_populations_mean$time, species = "Fox", count = predicted_populations_mean$fox)
)
fit_df$species <- factor(fit_df$species, levels = c("Rabbit", "Fox"))
p_pop <- ggplot() +
    geom_ribbon(
        data = fit_ranges_df,
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
    scale_fill_manual(values = c(Rabbit = color_rabbit, Fox = color_fox), guide = "none") +
    scale_colour_manual(values = c(Rabbit = color_rabbit, Fox = color_fox), name = NULL) +
    labs(x = "Time", y = "Count") +
    theme_bw(base_size = 20) +
    theme(
        plot.title = element_blank(), panel.grid = element_blank(),
        legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "grey35", linewidth = 0.5),
        legend.margin = margin(6, 8, 6, 8)
    )
ggsave(
    filename = "Interactive_session_6_population_dynamics.png",
    plot = p_pop,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
#-----------------------------------------------------------------------Plot inferred parameter values
long_post <- do.call(rbind, lapply(c("alpha", "beta", "gamma", "delta"), function(p) {
    data.frame(parameter = p, value = accepted_parameters[[p]], stringsAsFactors = FALSE)
}))
long_post$parameter <- factor(long_post$parameter, levels = c("alpha", "beta", "gamma", "delta"))
p_params <- ggplot(long_post, aes(x = value)) +
    geom_density(fill = "#3182BD", alpha = 0.45, colour = "#08519C", linewidth = 0.6) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    coord_cartesian(xlim = c(0, 2)) +
    labs(x = NULL, y = "Density") +
    theme_bw(base_size = 20)
ggsave(
    filename = "Interactive_session_6_parameters.png",
    plot = p_params,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
