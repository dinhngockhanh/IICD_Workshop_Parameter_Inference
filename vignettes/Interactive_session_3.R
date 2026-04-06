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
N_particles <- 10000
N_burn_in <- 2000
Gaussian_noise_sd <- 200
Gaussian_perturbation_sd <- 0.5
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
dprior <- function(parameters) {
    prod(dunif(
        x = c(parameters$alpha, parameters$beta, parameters$gamma, parameters$delta),
        min = as.numeric(parameter_ranges["min", ]), max = as.numeric(parameter_ranges["max", ])
    ))
}
rprior <- function(N_parameters = 1) {
    data.frame(
        alpha = runif(N_parameters, min = parameter_ranges["min", "alpha"], max = parameter_ranges["max", "alpha"]),
        beta  = runif(N_parameters, min = parameter_ranges["min", "beta"], max = parameter_ranges["max", "beta"]),
        gamma = runif(N_parameters, min = parameter_ranges["min", "gamma"], max = parameter_ranges["max", "gamma"]),
        delta = runif(N_parameters, min = parameter_ranges["min", "delta"], max = parameter_ranges["max", "delta"])
    )
}
dperturb <- function(parameters, parameters_previous) {
    prod(dnorm(
        x = c(parameters$alpha, parameters$beta, parameters$gamma, parameters$delta),
        mean = c(parameters_previous$alpha, parameters_previous$beta, parameters_previous$gamma, parameters_previous$delta),
        sd = Gaussian_perturbation_sd
    ))
}
rperturb <- function(parameters_previous) {
    data.frame(
        alpha = rnorm(1, mean = parameters_previous$alpha, sd = Gaussian_perturbation_sd),
        beta = rnorm(1, mean = parameters_previous$beta, sd = Gaussian_perturbation_sd),
        gamma = rnorm(1, mean = parameters_previous$gamma, sd = Gaussian_perturbation_sd),
        delta = rnorm(1, mean = parameters_previous$delta, sd = Gaussian_perturbation_sd)
    )
}
loglikelihood <- function(parameters, initial_state, data) {
    predicted_populations <- Lotka_Volterra_ODE_solver(
        parameters = parameters,
        times = data$time,
        initial_state = initial_state
    )
    residual_sum <- sum((data$rabbit - predicted_populations$rabbit)^2 + (data$fox - predicted_populations$fox)^2)
    n <- 2 * nrow(data)
    loglikelihood <- -(n / 2) * log(2 * pi * Gaussian_noise_sd * Gaussian_noise_sd) - residual_sum / (2 * Gaussian_noise_sd * Gaussian_noise_sd)
    return(loglikelihood)
}
mcmc <- function() {
    mcmc_chain <- data.frame()
    current_particle <- rprior()
    mcmc_chain <- rbind(mcmc_chain, current_particle)
    pb <- txtProgressBar(
        min = 0, max = N_particles,
        style = 3, width = 50, char = "+"
    )
    for (i in 1:N_particles) {
        setTxtProgressBar(pb, i)
        proposed_particle <- rperturb(current_particle)
        acceptance_threshold <- runif(1)
        proposed_particle_prior <- dprior(proposed_particle)
        if (proposed_particle_prior == 0) {
            acceptance_probability <- 0
        } else {
            acceptance_probability <- min(1, exp(
                loglikelihood(proposed_particle, initial_state, data) - loglikelihood(current_particle, initial_state, data) +
                    log(dprior(proposed_particle)) - log(dprior(current_particle)) +
                    log(dperturb(current_particle, proposed_particle)) - log(dperturb(proposed_particle, current_particle))
            ))
        }
        if (acceptance_probability > acceptance_threshold) {
            current_particle <- proposed_particle
        }
        mcmc_chain <- rbind(mcmc_chain, current_particle)
    }
    mcmc_chain <- mcmc_chain[N_burn_in + 1:nrow(mcmc_chain), ]
    return(mcmc_chain)
}
mcmc_chain <- mcmc()
predicted_populations <- lapply(1:nrow(mcmc_chain), function(i) {
    predicted_populations <- Lotka_Volterra_ODE_solver(
        parameters = mcmc_chain[i, ],
        times = seq(0, 14, by = 0.1),
        initial_state = initial_state
    )
    predicted_populations$rabbit <- predicted_populations$rabbit + rnorm(nrow(predicted_populations), mean = 0, sd = Gaussian_noise_sd)
    predicted_populations$fox <- predicted_populations$fox + rnorm(nrow(predicted_populations), mean = 0, sd = Gaussian_noise_sd)
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
    filename = "Interactive_session_3_population_dynamics.png",
    plot = p_pop,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
#-----------------------------------------------------------------------Plot inferred parameter values
long_post <- do.call(rbind, lapply(c("alpha", "beta", "gamma", "delta"), function(p) {
    data.frame(parameter = p, value = mcmc_chain[[p]], stringsAsFactors = FALSE)
}))
long_post$parameter <- factor(long_post$parameter, levels = c("alpha", "beta", "gamma", "delta"))
p_params <- ggplot(long_post, aes(x = value)) +
    geom_density(fill = "#3182BD", alpha = 0.45, colour = "#08519C", linewidth = 0.6) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    coord_cartesian(xlim = c(0, 2)) +
    labs(x = NULL, y = "Density") +
    theme_bw(base_size = 20)
ggsave(
    filename = "Interactive_session_3_parameters.png",
    plot = p_params,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
