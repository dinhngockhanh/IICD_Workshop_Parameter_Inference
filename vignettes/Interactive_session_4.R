setwd( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/R"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
files_sources <- list.files(pattern = "\\.[rR]$") #<<<<<<<<<<<<<<<<<<<<<
sapply(files_sources, source) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
setwd( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/vignettes"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
library(ggplot2)
data <- read.csv( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/data/data.csv"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
set.seed(1)
N_particles <- 12500
N_accept <- 500
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
rprior <- function(N_parameters = 1) {
    data.frame(
        alpha = runif(N_parameters, min = parameter_ranges["min", "alpha"], max = parameter_ranges["max", "alpha"]),
        beta  = runif(N_parameters, min = parameter_ranges["min", "beta"], max = parameter_ranges["max", "beta"]),
        gamma = runif(N_parameters, min = parameter_ranges["min", "gamma"], max = parameter_ranges["max", "gamma"]),
        delta = runif(N_parameters, min = parameter_ranges["min", "delta"], max = parameter_ranges["max", "delta"])
    )
}
candidate_parameters <- rprior(N_parameters = N_particles)
candidate_errors <- rep(NA, N_particles)
pb <- txtProgressBar(
    min = 0, max = N_particles,
    style = 3, width = 50, char = "+"
)
for (i in 1:N_particles) {
    setTxtProgressBar(pb, i)
    candidate_populations <- tau_leaping(
        initial_state = initial_state,
        parameters = candidate_parameters[i, ],
        reaction_propensities = reaction_propensities,
        reaction_stoichiometries = reaction_stoichiometries,
        time_points = data$time,
        tau = 0.01, max_pop = 5000L
    )
    candidate_errors[i] <- sum((data$rabbit - candidate_populations[, "rabbit"])^2 + (data$fox - candidate_populations[, "fox"])^2)
}
accepted_indices <- order(candidate_errors)[1:N_accept]
accepted_parameters <- candidate_parameters[accepted_indices, ]
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
    filename = "Interactive_session_4_population_dynamics.png",
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
    filename = "Interactive_session_4_parameters.png",
    plot = p_params,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
